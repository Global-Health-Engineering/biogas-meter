import os
import git
import json
import ntpath
import argparse
import numpy as np
import pandas as pd
from datetime import datetime
import CoolProp as CP


class MFCperInterval(object):
    df: pd.DataFrame
    metaData: str
    
    def __init__(self, df, metaData):
        self.df = df
        with open(metaData) as fh:
            self.meta = json.load(fh)
        self.flds = tuple(v["FLD"] for v in self.meta["mfc"].values())
        self.HEOS = tuple(CP.AbstractState("HEOS", fld) for fld in self.flds)

    def get_rhomass(self, idx, p_bar, T_degC):
        """ return fluid density in g/L """
        self.HEOS[idx].update(CP.PT_INPUTS, p_bar*1e5, T_degC+273.15)
        return self.HEOS[idx].rhomass()

    def get_rhomolar(self, idx, p_bar, T_degC):
        """ return fluid density in mol/L """
        self.HEOS[idx].update(CP.PT_INPUTS, p_bar*1e5, T_degC+273.15)
        return self.HEOS[idx].rhomolar() / 1e3

    def get_all_info(self):
        # make empty dictionary
        d = {}

        # get number of measurements per stream
        for k in self.meta["mfc"].keys():
            d[f"{k}: ground truth count"] = len(self.df[f"{k}: Flow (ln/min)"].dropna())

        for col in self.df.columns:
            if "setpoint" not in col.lower():
                # get mean values per interval
                d[col] = self.df[col].dropna().mean()
                # get standard deviation per interval
                d[f"{col} STD"] = self.df[col].dropna().std()
                # get standard uncertainty per interval
                d[f"{col} STUNC"] = d[f"{col} STD"] / (d[f"{col.split(': ')[0]}: ground truth count"])**0.5
        
        # calculate total volumetric flow
        d["Total vol. flow (ln/min)"] = sum([self.df[f"{k}: Flow (ln/min)"].dropna().mean() for k in self.meta["mfc"].keys()])
        
        # calculate relative combined uncertainty of total volumetric flow
        if d["Total vol. flow (ln/min)"] != 0:
            d["Total vol. flow (ln/min) COMB-UNC"] = sum([d[f"{k}: Flow (ln/min) STUNC"]**2 for k in self.meta["mfc"].keys()])**0.5
        else:
            d["Total vol. flow (ln/min) COMB-UNC"] = 0
        
        for k in self.meta["mfc"].keys():
            # calculate volumetric fraction
            if "Total vol. flow (ln/min)" != 0:
                d[f"{k}: Vol. fraction"] = d[f"{k}: Flow (ln/min)"] / d["Total vol. flow (ln/min)"]
            else:
                d[f"{k}: Vol. fraction"] = 0
        
            # calculate relative combined uncertainty of volumatric fraction
            if d[f"{k}: Flow (ln/min)"] != 0:
                d[f"{k}: Vol. fraction COMB-UNC"] = (
                (d[f"{k}: Flow (ln/min) STUNC"] / d[f"{k}: Flow (ln/min)"])**2
                + (d["Total vol. flow (ln/min) COMB-UNC"] / d["Total vol. flow (ln/min)"])**2)**0.5
            else:
                d[f"{k}: Vol. fraction COMB-UNC"] = 0
        return d

    def get_mass_fractions(self, d):
        for k, v in self.meta["mfc"].items():
            rhomass = self.df.apply(lambda x: self.get_rhomass(idx = v["ID"],
                                                               p_bar = x[f"{k}: Pressure outlet (bar(a))"],
                                                               T_degC = x[f"{k}: Temperature (°C)"]), axis=1)
            self.df[f"{k}: Mass flow (g/min)"] = self.df[f"{k}: Flow (ln/min)"] * rhomass
            d[f"{k}: Mass flow (g/min)"] = self.df[f"{k}: Mass flow (g/min)"].mean()
        self.df["Total mass flow (g/min)"] = sum([self.df[f"{k}: Mass flow (g/min)"] for k in self.meta["mfc"].keys()])
        d["Total mass flow (g/min)"] = self.df["Total mass flow (g/min)"].mean()
        
        for k in self.meta["mfc"].keys():
            # calculate mass fraction
            self.df[f"{k}: Mass fraction"] = 0
            self.df.loc[df["Total mass flow (g/min)"] != 0, f"{k}: Mass fraction"] = self.df[f"{k}: Mass flow (g/min)"] / self.df["Total mass flow (g/min)"]
            d[f"{k}: Mass fraction"] = self.df[f"{k}: Mass fraction"].mean()
        return d

    def get_molar_fractions(self, d):
        for k, v in self.meta["mfc"].items():
            rhomolar = self.df.apply(lambda x: self.get_rhomolar(idx = v["ID"],
                                                                 p_bar = x[f"{k}: Pressure outlet (bar(a))"],
                                                                 T_degC = x[f"{k}: Temperature (°C)"]), axis=1)
            self.df[f"{k}: Molar flow (mol/min)"] = self.df[f"{k}: Flow (ln/min)"] * rhomolar
        self.df["Total molar flow (mol/min)"] = sum([self.df[f"{k}: Molar flow (mol/min)"] for k in self.meta["mfc"].keys()])
        d["Total molar flow (mol/min)"] = self.df["Total molar flow (mol/min)"].mean()
        
        for k in self.meta["mfc"].keys():
            # calculate molar fraction
            self.df[f"{k}: Mole fraction"] = 0
            self.df.loc[df["Total molar flow (mol/min)"] != 0, f"{k}: Mole fraction"] = self.df[f"{k}: Molar flow (mol/min)"] / self.df["Total molar flow (mol/min)"]
            d[f"{k}: Mole fraction"] = self.df[f"{k}: Mole fraction"].mean()
        return d

    def __call__(self, mass_props=True, molar_props=True):
        d = self.get_all_info()
        if mass_props:
            d = self.get_mass_fractions(d)
        if molar_props:
            d = self.get_molar_fractions(d)

        unique_flds = True if len(set(self.flds)) == len(self.flds) else False
        
        # create column mapper to swap mass flow controller S/N with fluid name
        name_map = {}
        for key, val in self.meta["mfc"].items():
            for col in d.keys():
                if key in col:
                    new_col_name = val["FLD"] if unique_flds else f'{val["FLD"]}_{val["ID"]}'
                    name_map[col] = f"{new_col_name} {col.lstrip(f'{key}:').strip(' ')}"
        
        # rename keys (swap mass flow controller S/N with fluid name)
        d_renamed = {}
        for k, v in d.items():
            try:
                d_renamed[name_map[k]] = v
            except KeyError:
                d_renamed[k] = v

        return d_renamed


class SBGperInterval(object):
    df: pd.DataFrame

    def __init__(self, df):
        self.df = df

    def __call__(self):
        d = {}
        
        # get point count
        d["SBG count"] = len(self.df["Flow (lph)"])
        # get mean pressure
        d["Pressure SBG [Pa]"] = self.df["Pressure (Pa)"].mean()
        # get standard deviation of pressure
        d["Pressure SBG [Pa] STD"] = self.df["Pressure (Pa)"].std()
        # get standard uncertainty of pressure
        d["Pressure SBG [Pa] STUNC"] = d["Pressure SBG [Pa] STD"] / len(self.df["Pressure (Pa)"])**0.5

        # convert from Standard Liters Per Hour to Normal Liters Per Minute
        # 1 SLPM = 1 NLPM * (273.15 K / 293.15 K) * (14.696 psi / 14.504 psi)
        # write it to temporary dataframe
        self.df["Flow [ln/min]"] = self.df["Flow (lph)"] / 60 * 293.15 / 273.15 * 14.504 / 14.696
        # get mean flow
        d["Flow SBG [ln/min]"] = self.df["Flow [ln/min]"].dropna().mean()
        # get standard deviation of flow
        d["Flow SBG [ln/min] STD"] = self.df["Flow [ln/min]"].dropna().std()
        # get standard uncertainty of flow
        d["Flow SBG [ln/min] STUNC"] = d["Flow SBG [ln/min] STD"] / len(self.df["Flow (lph)"].dropna())**0.5

        return d


class AnalyzeAll(object):
    mfc_dir: str
    sbg_dir: str
    metaData: str

    def __init__(self, mfc_dir, sbg_dir, metaData):
        self.mfc_dir = mfc_dir
        self.sbg_dir = sbg_dir
        self.metaData = metaData
        with open(metaData) as fh:
            self.meta = json.load(fh)

    def get_df(self, df_type, fname):
        assert df_type.lower() in ["mfc", "sbg"]
        d = {"mfc": {"path": self.mfc_dir,
                     "col": "TimeStamp"},
             "sbg": {"path": self.sbg_dir,
                     "col": "Local Time"}}
        path = d[df_type.lower()]["path"]
        col = d[df_type.lower()]["col"]
        # read data
        df = pd.read_csv(os.path.join(path, fname))
        # convert strings to datetime objects
        df[col] = pd.to_datetime(df[col])
        # remove time zone info
        df[col] = df[col].apply(lambda x: x.replace(tzinfo=None))
        # make "TimeStamp" column new index column
        df.set_index(col, drop=True, inplace=True)
        # drop index name
        df.index.name = None
        return df

    def get_mfc_df(self, fname):
        return self.get_df("mfc", fname)

    def get_sbg_df(self, fname):
        df = self.get_df("sbg", fname)
        # drop "UNIX Timestamp" column
        df.drop("UNIX Timestamp", axis=1, inplace=True)
        return df

    def get_interval_df(self, df, stamp):
        time_start = datetime.strptime(f"{stamp['date/yyyy-mm-dd']} {stamp['time_start/hh:mm:ss']}", "%Y-%m-%d %H:%M:%S")
        time_end = datetime.strptime(f"{stamp['date/yyyy-mm-dd']} {stamp['time_end/hh:mm:ss']}", "%Y-%m-%d %H:%M:%S")
        return df[(df.index >= time_start) & (df.index <= time_end)].copy()

    def __call__(self, mass_props=True, molar_props=True, derived_dir=None):
        date = None
        df = pd.DataFrame()
        for m in self.meta["measurements"]:
            for i, stamp in enumerate(m["stamps"]):
                if stamp["date/yyyy-mm-dd"] != date:
                    df_mfc = self.get_mfc_df(f"{stamp['date/yyyy-mm-dd']}.csv")
                    df_sbg = self.get_sbg_df(f"{stamp['date/yyyy-mm-dd']}.csv")
                df_mfc_i = self.get_interval_df(df_mfc, stamp)
                df_sbg_i = self.get_interval_df(df_sbg, stamp)

                mfcperinterval = MFCperInterval(df=df_mfc_i, metaData=self.metaData)
                sbgperinterval = SBGperInterval(df=df_sbg_i)

                d = {}
                d["date/yyyy-mm-dd"] = stamp["date/yyyy-mm-dd"]
                d["time_start/hh:mm:ss"] = stamp["time_start/hh:mm:ss"]
                d["time_end/hh:mm:ss"] = stamp["time_end/hh:mm:ss"]
                # get ground truth means, stds, uncertainties, etc. for given interval
                d.update(mfcperinterval(mass_props=mass_props, molar_props=molar_props))
                # get Smart Biogas means and stds
                d.update(sbgperinterval())
                _df = pd.DataFrame(d, index=[i])
                df = pd.concat([df, _df], axis=0)

        if derived_dir:
            fname = os.path.basename(self.metaData)
            df.to_csv(os.path.join(derived_dir, f"{fname.split('.')[0]}.csv"), sep=",", index=False)
        return df


def get_git_root(path):
    git_repo = git.Repo(path, search_parent_directories=True)
    return git_repo.working_dir


def main():
    metaData = os.path.join(get_git_root(os.getcwd()), "data", "metadata", "test_exp.json")

    mfc_dir = os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "mfc")
    sbg_dir = os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "smart_biogas")
    derived_dir = os.path.join(get_git_root(os.getcwd()), "data", "derived_data")

    an = AnalyzeAll(mfc_dir = mfc_dir,
                    sbg_dir = sbg_dir,
                    metaData = metaData)
    an(mass_props = False,
       molar_props = False,
       derived_dir = derived_dir)


if __name__ == "__main__":
    main()
