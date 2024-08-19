import os
import json
import ntpath
import numpy as np
import pandas as pd
import CoolProp as CP
from get_git_root import get_git_root

class MFCPerInterval(object):
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
            if col.split(': ')[0] in self.meta["mfc"].keys():
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

    def __call__(self, mass_props=False, molar_props=False):
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


class SBGPerInterval(object):
    df: pd.DataFrame

    def __init__(self, df):
        self.df = df

    def __call__(self):
        d = {}
        
        # get point count
        d["SBG count"] = len(self.df["Flow (lph)"])
        # get mean pressure
        d["Pressure SBG (Pa)"] = self.df["Pressure (Pa)"].mean()
        # get standard deviation of pressure
        d["Pressure SBG (Pa) STD"] = self.df["Pressure (Pa)"].std()
        # get standard uncertainty of pressure
        d["Pressure SBG (Pa) STUNC"] = d["Pressure SBG (Pa) STD"] / len(self.df["Pressure (Pa)"])**0.5

        # convert from Standard Liters Per Hour to Normal Liters Per Minute
        # 1 SLPM = 1 NLPM * (273.15 K / 293.15 K) * (14.696 psi / 14.504 psi)
        self.df["Flow (ln/min)"] = self.df["Flow (lph)"] / 60 * 293.15 / 273.15 * 14.504 / 14.696
        # get mean flow
        d["Flow SBG (ln/min)"] = self.df["Flow (ln/min)"].dropna().mean()
        # get standard deviation of flow
        d["Flow SBG (ln/min) STD"] = self.df["Flow (ln/min)"].dropna().std()
        # get standard uncertainty of flow
        d["Flow SBG (ln/min) STUNC"] = d["Flow SBG (ln/min) STD"] / len(self.df["Flow (lph)"].dropna())**0.5

        return d


class CleanSensirion(object):
    """ object converting and cleaning edf file to df/csv """
    edf_file: str

    def __init__(self, edf_file):
        self.edf_file = edf_file
        self.df = pd.read_csv(edf_file,
                              engine = "python",
                              encoding = "utf-8",
                              sep = r"\t",
                              skiprows = 10)
        self.df = self.clean_df(self.df)
        self.df = self.map_cols(self.df)

    def clean_df(self, df):
        df.replace({np.nan: None})
        df.drop("Epoch_UTC", axis=1, inplace=True)
        df.Local_Date_Time = pd.to_datetime(df.Local_Date_Time)
        df.Local_Date_Time = df.Local_Date_Time.apply(lambda x: x.replace(tzinfo=None))
        # df.Local_Date_Time = df.Local_Date_Time.dt.round("1s")
        df.index = df.Local_Date_Time
        df.drop("Local_Date_Time", axis=1, inplace=True)
        df.index.name = None
        return df

    def map_cols(self, df):
        mapper = {}
        for col in list(df.columns):
            if col.split("_")[0] == "F":
                mapper[col] = "Sensirion Vol. flow (ls/min)"
            elif col.split("_")[0] == "T":
                mapper[col] = "Sensirion T (deg C)"
        df.rename(mapper, axis=1, inplace=True)
        return df

    def get_df(self):
        return self.df

    def to_csv(self, derived_data_dir):
        self.df.to_csv(os.path.join(derived_data_dir,
                                    f"sensirion_{ntpath.basename(self.edf_file.strip('.edf'))}.csv"))


class SensirionPerInterval(object):
    df: pd.DataFrame

    def __init__(self, df):
        self.df = df

    def __call__(self):
        d = {}
        d["Sensirion count"] = len([j for j in self.df["Sensirion Vol. flow (ls/min)"] if j is not None])
        for col in self.df.columns:
            d[col] = self.df[col].mean()
            d[f"{col} STD"] = self.df[col].std()
            d[f"{col} STUNC"] = d[f"{col} STD"] / (d["Sensirion count"])**0.5
        return d


class HumidityPerInterval(object):
    df: pd.DataFrame

    def __init__(self, df):
        self.df = df

    def __call__(self):
        d = {}
        d["AHT count"] = len([j for j in self.df["AHT temp (deg C)"] if j is not None])
        for col in self.df.columns:
            d[col] = self.df[col].mean()
            d[f"{col} STD"] = self.df[col].std()
            d[f"{col} STUNC"] = d[f"{col} STD"] / (d["AHT count"])**0.5
        return d

def main():
    edf_file = os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "sensirion", "2024-08-19.edf")
    derived_data_dir = os.path.join(get_git_root(os.getcwd()), "data", "derived_data")
    s = CleanSensirion(edf_file)
    s.to_csv(derived_data_dir)


if __name__ == "__main__":
    main()
