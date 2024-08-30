import os
import json
import ntpath
import numpy as np
import pandas as pd
import CoolProp as CP
from datetime import datetime
from get_git_root import get_git_root


class MassFlowControllers(object):
    meta_data: str
    raw_data_dir: str
    
    def __init__(self,
                 meta_data,
                 raw_data_dir = os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "mfc")):
        self.meta_data = meta_data
        self.raw_data_dir = raw_data_dir
        with open(meta_data) as fh:
            self.meta = json.load(fh)
        self.flds = tuple(v for v in self.meta["mfc"]["sensors"].values())
        self.HEOS = {fld: CP.AbstractState("HEOS", fld) for fld in self.flds}

    def set_df(self, fname):
        self.df = pd.read_csv(os.path.join(self.raw_data_dir, fname))
        # convert strings to datetime objects
        self.df["TimeStamp"] = pd.to_datetime(self.df["TimeStamp"])
        # remove time zone info
        self.df["TimeStamp"] = self.df["TimeStamp"].apply(lambda x: x.replace(tzinfo=None))
        # make "TimeStamp" column new index column
        self.df.set_index("TimeStamp", drop=True, inplace=True)
        # drop index name
        self.df.index.name = None
        # standardize flow values
        self.set_standard_flow()

    def set_standard_flow(self):
        """ standardize all flow values in raw data """
        for col in self.df.columns:
            if "Flow" in col:
                assert col.split(' ')[-1] == "(ln/min)"
                newCol = f"{' '.join(col.split(' ')[:-1])} (ls/min)"
                self.df[newCol] = self.get_standard_flow(self.df[col])
                self.df.drop(col, axis=1, inplace=True)

    def get_standard_flow(self, V_lnpm):
        """ return flow at standard conditions (273.15 K and 1 bar(a)) """
        p_mfc = self.meta["mfc"]["normalizing conditions"]["p_bara"]
        T_mfc = self.meta["mfc"]["normalizing conditions"]["T_K"]
        return V_lnpm * p_mfc / 1 * 273.15 / T_mfc

    def get_rhomass(self, fld, p_bar, T_degC):
        """ return fluid density in g/L """
        try:
            self.HEOS[fld].update(CP.PT_INPUTS, p_bar*1e5, T_degC+273.15)
            return self.HEOS[fld].rhomass()
        except ValueError:
            return np.nan

    def get_rhomolar(self, fld, p_bar, T_degC):
        """ return fluid density in mol/L """
        try:
            self.HEOS[fld].update(CP.PT_INPUTS, p_bar*1e5, T_degC+273.15)
            return self.HEOS[fld].rhomolar() / 1e3
        except ValueError:
            return np.nan

    def get_mean_std_stunc(self, df, d_res=None):
        # make empty dictionary
        if d_res is None:
            d_res = {}

        # get number of measurements per stream
        for k in self.meta["mfc"]["sensors"].keys():
            d_res[f"{k}: ground truth count"] = len(df[f"{k}: Flow (ls/min)"].dropna())

        for col in df.columns:
            if col.split(': ')[0] in self.meta["mfc"]["sensors"].keys():
                if "setpoint" not in col.lower():
                    # get mean values per interval
                    d_res[col] = df[col].dropna().mean()
                    # get standard deviation per interval
                    d_res[f"{col} STD"] = df[col].dropna().std()
                    # get standard uncertainty per interval
                    d_res[f"{col} STUNC"] = d_res[f"{col} STD"] / (d_res[f"{col.split(': ')[0]}: ground truth count"])**0.5
        
        # calculate total volumetric flow
        d_res["Total vol. flow (ls/min)"] = sum([df[f"{k}: Flow (ls/min)"].dropna().mean() for k in self.meta["mfc"]["sensors"].keys()])
        
        # calculate relative combined uncertainty of total volumetric flow
        if d_res["Total vol. flow (ls/min)"] != 0:
            d_res["Total vol. flow (ls/min) COMB-UNC"] = sum([d_res[f"{k}: Flow (ls/min) STUNC"]**2 for k in self.meta["mfc"]["sensors"].keys()])**0.5
        else:
            d_res["Total vol. flow (ls/min) COMB-UNC"] = 0
        
        for k in self.meta["mfc"]["sensors"].keys():
            # calculate volumetric fraction
            if "Total vol. flow (ls/min)" != 0:
                d_res[f"{k}: Vol. fraction"] = d_res[f"{k}: Flow (ls/min)"] / d_res["Total vol. flow (ls/min)"]
            else:
                d_res[f"{k}: Vol. fraction"] = 0
        
            # calculate relative combined uncertainty of volumatric fraction
            if d_res[f"{k}: Flow (ls/min)"] != 0:
                d_res[f"{k}: Vol. fraction COMB-UNC"] = (
                (d_res[f"{k}: Flow (ls/min) STUNC"] / d_res[f"{k}: Flow (ls/min)"])**2
                + (d_res["Total vol. flow (ls/min) COMB-UNC"] / d_res["Total vol. flow (ls/min)"])**2)**0.5
            else:
                d_res[f"{k}: Vol. fraction COMB-UNC"] = 0
        return d_res

    def get_mass_fractions(self, df, d_res):
        for k, v in self.meta["mfc"]["sensors"].items():
            rhomass = df.apply(lambda x: self.get_rhomass(fld = v,
                                                          p_bar = x[f"{k}: Pressure outlet (bar(a))"],
                                                          T_degC = x[f"{k}: Temperature (°C)"]), axis=1)
            df[f"{k}: Mass flow (g/min)"] = df[f"{k}: Flow (ls/min)"] * rhomass
            d_res[f"{k}: Mass flow (g/min)"] = df[f"{k}: Mass flow (g/min)"].mean()
        df["Total mass flow (g/min)"] = sum([df[f"{k}: Mass flow (g/min)"] for k in self.meta["mfc"]["sensors"].keys()])
        d_res["Total mass flow (g/min)"] = df["Total mass flow (g/min)"].mean()
        
        for k in self.meta["mfc"]["sensors"].keys():
            # calculate mass fraction
            df[f"{k}: Mass fraction"] = 0
            df.loc[df["Total mass flow (g/min)"] != 0, f"{k}: Mass fraction"] = df[f"{k}: Mass flow (g/min)"] / df["Total mass flow (g/min)"]
            d_res[f"{k}: Mass fraction"] = df[f"{k}: Mass fraction"].mean()
        return d_res

    def get_molar_fractions(self, df, d_res):
        for k, v in self.meta["mfc"]["sensors"].items():
            rhomolar = df.apply(lambda x: self.get_rhomolar(fld = v,
                                                            p_bar = x[f"{k}: Pressure outlet (bar(a))"],
                                                            T_degC = x[f"{k}: Temperature (°C)"]), axis=1)
            df[f"{k}: Molar flow (mol/min)"] = df[f"{k}: Flow (ls/min)"] * rhomolar
        df["Total molar flow (mol/min)"] = sum([self.df[f"{k}: Molar flow (mol/min)"] for k in self.meta["mfc"]["sensors"].keys()])
        d_res["Total molar flow (mol/min)"] = df["Total molar flow (mol/min)"].mean()
        
        for k in self.meta["mfc"]["sensors"].keys():
            # calculate molar fraction
            df[f"{k}: Mole fraction"] = 0
            df.loc[df["Total molar flow (mol/min)"] != 0, f"{k}: Mole fraction"] = df[f"{k}: Molar flow (mol/min)"] / df["Total molar flow (mol/min)"]
            d_res[f"{k}: Mole fraction"] = df[f"{k}: Mole fraction"].mean()
        return d_res

    def get_interval_df(self, stamp):
        time_start = datetime.strptime(f"{stamp['date/yyyy-mm-dd']} {stamp['time_start/hh:mm:ss']}", "%Y-%m-%d %H:%M:%S")
        time_end = datetime.strptime(f"{stamp['date/yyyy-mm-dd']} {stamp['time_end/hh:mm:ss']}", "%Y-%m-%d %H:%M:%S")
        return self.df[(self.df.index >= time_start) & (self.df.index <= time_end)].copy()

    def __call__(self, mass_props=False, molar_props=False, derived_dir=None):
        date = None
        df_res = pd.DataFrame()
        for m in self.meta["measurements"]:
            for i, stamp in enumerate(m["stamps"]):
                if stamp["date/yyyy-mm-dd"] != date:
                    date = stamp["date/yyyy-mm-dd"]
                    self.set_df(f"{date}.csv")
                df_i = self.get_interval_df(stamp)

                d_res = {}
                d_res["date/yyyy-mm-dd"] = stamp["date/yyyy-mm-dd"]
                d_res["time_start/hh:mm:ss"] = stamp["time_start/hh:mm:ss"]
                d_res["time_end/hh:mm:ss"] = stamp["time_end/hh:mm:ss"]
                # get ground truth means, stds, uncertainties, etc. for given interval
                d_res = self.get_mean_std_stunc(df_i, d_res=d_res)
                if mass_props:
                    d_res = self.get_mass_fractions(df_i, d_res=d_res)
                if molar_props:
                    d_res = self.get_molar_fractions(df_i, d_res=d_res)
                
                _df = pd.DataFrame(d_res, index=[i])
                df_res = pd.concat([df_res, _df], axis=0)

        unique_flds = True if len(set(self.flds)) == len(self.flds) else False
        # create column mapper to swap mass flow controller S/N with fluid name
        name_map = {}
        for i, (k, v) in enumerate(self.meta["mfc"]["sensors"].items()):
            for col in df_res.columns:
                if k in col:
                    new_col_name = v if unique_flds else f'{v}_{i}'
                    name_map[col] = f"{new_col_name} {col.lstrip(f'{k}:').strip(' ')}"
        df_res.rename(columns=name_map, inplace=True)

        if derived_dir:
            fname = os.path.basename(self.meta_data)
            df_res.to_csv(os.path.join(derived_dir, f"{fname.split('.')[0]}_mfc.csv"), sep=",", index=False)
        return df_res


class SmartBiogas(object):
    meta_data: str
    raw_data_dir: str

    def __init__(self,
                 meta_data,
                 raw_data_dir = os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "sbg")):
        self.meta_data = meta_data
        with open(meta_data) as fh:
            self.meta = json.load(fh)
        self.raw_data_dir = raw_data_dir

    def set_df(self, fname):
        self.df = pd.read_csv(os.path.join(self.raw_data_dir, fname))
        # convert strings to datetime objects
        self.df["Local Time"] = pd.to_datetime(self.df["Local Time"])
        # remove time zone info
        self.df["Local Time"] = self.df["Local Time"].apply(lambda x: x.replace(tzinfo=None))
        # make "TimeStamp" column new index column
        self.df.set_index("Local Time", drop=True, inplace=True)
        # drop index name
        self.df.index.name = None
        # standardize flow values
        self.set_standard_flow()

    def set_standard_flow(self):
        """ standardize all flow values in raw data """
        self.df["Flow (ls/min)"] = self.get_standard_flow(self.df["Flow (lph)"])
        self.df.drop("Flow (lph)", axis=1, inplace=True)

    def get_standard_flow(self, V_lph):
        """ return flow at standard conditions (273.15 K and 1 bar(a)) """
        p_sbg = self.meta["sbg"]["normalizing conditions"]["p_bara"]
        T_sbg = self.meta["sbg"]["normalizing conditions"]["T_K"]
        return V_lph / 60 * p_sbg / 1 * 273.15 / T_sbg

    def get_interval_df(self, stamp):
        time_start = datetime.strptime(f"{stamp['date/yyyy-mm-dd']} {stamp['time_start/hh:mm:ss']}", "%Y-%m-%d %H:%M:%S")
        time_end = datetime.strptime(f"{stamp['date/yyyy-mm-dd']} {stamp['time_end/hh:mm:ss']}", "%Y-%m-%d %H:%M:%S")
        return self.df[(self.df.index >= time_start) & (self.df.index <= time_end)].copy()

    def get_mean_std_stunc(self, df, d_res=None):
        # make empty dictionary if d_res is None
        if d_res is None:
            d_res = {}
        # get point count
        d_res["SBG count"] = len(df["Flow (ls/min)"])
        # iterate over two properties of interest and get means, std, stunc
        for col in ["Flow (ls/min)", "Pressure (Pa)"]:
            ncol = f"{col.split(' ')[0]} SBG {col.split(' ')[1]}"
            # get mean property
            d_res[ncol] = df[col].mean()
            # get standard deviation of a property
            d_res[f"{ncol} STD"] = df[col].std()
            # get standard uncertainty of a property
            d_res[f"{ncol} STUNC"] = d_res[f"{ncol} STD"] / len(df[col])**0.5
        # get only mean temperature, as not enough measurements are provided to calculate std and stunc
        d_res["Temperature SBG (°C)"] = df["Temperature (°C)"].mean()
        return d_res

    def __call__(self, derived_dir=None):
        date = None
        df_res = pd.DataFrame()
        for m in self.meta["measurements"]:
            for i, stamp in enumerate(m["stamps"]):
                if stamp["date/yyyy-mm-dd"] != date:
                    date = stamp["date/yyyy-mm-dd"]
                    self.set_df(f"{date}.csv")
                df_i = self.get_interval_df(stamp)

                d_res = {}
                d_res["date/yyyy-mm-dd"] = stamp["date/yyyy-mm-dd"]
                d_res["time_start/hh:mm:ss"] = stamp["time_start/hh:mm:ss"]
                d_res["time_end/hh:mm:ss"] = stamp["time_end/hh:mm:ss"]
                # get means, stds, and uncertainty for given interval
                d_res = self.get_mean_std_stunc(df_i, d_res=d_res)
                _df = pd.DataFrame(d_res, index=[i])
                df_res = pd.concat([df_res, _df], axis=0)
        if derived_dir:
            fname = os.path.basename(self.meta_data)
            df_res.to_csv(os.path.join(derived_dir, f"{fname.split('.')[0]}_sbg.csv"), sep=",", index=False)
        return df_res


class CleanSensirionMeasurements(object):
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


class Sensirion(object):
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


class AHT20(object):
    meta_data: str
    raw_data_dir: str

    def __init__(self,
                 meta_data,
                 raw_data_dir = os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "rH-temp")):
        self.meta_data = meta_data
        with open(meta_data) as fh:
            self.meta = json.load(fh)
        self.raw_data_dir = raw_data_dir

    def set_df(self, fname):
        self.df = pd.read_csv(os.path.join(self.raw_data_dir, fname))
        # convert strings to datetime objects
        self.df["Datetime"] = pd.to_datetime(self.df["Datetime"])
        # make "Datetime" column new index column
        self.df.set_index("Datetime", drop=True, inplace=True)
        # drop index name
        self.df.index.name = None
        # rename columns
        name_map = {col: f"AHT {col}" for col in self.df.columns}
        self.df.rename(columns=name_map, inplace=True)

    def get_interval_df(self, stamp):
        time_start = datetime.strptime(f"{stamp['date/yyyy-mm-dd']} {stamp['time_start/hh:mm:ss']}", "%Y-%m-%d %H:%M:%S")
        time_end = datetime.strptime(f"{stamp['date/yyyy-mm-dd']} {stamp['time_end/hh:mm:ss']}", "%Y-%m-%d %H:%M:%S")
        return self.df[(self.df.index >= time_start) & (self.df.index <= time_end)].copy()

    def get_mean_std_stunc(self, df, d_res=None):
        if d_res is None:
            d_res = {}
        d_res["AHT count"] = len([i for i in df["AHT Temp (deg C)"] if i is not np.nan])
        for col in self.df.columns:
            d_res[col] = df[col].mean()
            d_res[f"{col} STD"] = df[col].std()
            d_res[f"{col} STUNC"] = d_res[f"{col} STD"] / (d_res["AHT count"])**0.5
        return d_res

    def __call__(self, derived_dir=None):
        date = None
        df_res = pd.DataFrame()
        for i, m in enumerate(self.meta["measurements"]):
            for j, stamp in enumerate(m["stamps"]):
                if stamp["date/yyyy-mm-dd"] != date:
                    date = stamp["date/yyyy-mm-dd"]
                    self.set_df(f"{date}.csv")
                df_j = self.get_interval_df(stamp)

                d_res = {}
                d_res["date/yyyy-mm-dd"] = stamp["date/yyyy-mm-dd"]
                d_res["time_start/hh:mm:ss"] = stamp["time_start/hh:mm:ss"]
                d_res["time_end/hh:mm:ss"] = stamp["time_end/hh:mm:ss"]
                # get means, stds, and uncertainty for given interval
                d_res = self.get_mean_std_stunc(df_j, d_res=d_res)
                _df = pd.DataFrame(d_res, index=[i+j])
                df_res = pd.concat([df_res, _df], axis=0)
        if derived_dir:
            fname = os.path.basename(self.meta_data)
            df_res.to_csv(os.path.join(derived_dir, f"{fname.split('.')[0]}_aht20.csv"), sep=",", index=False)
        return df_res


class CompositionSensor(object):
    meta_data: str
    raw_data_dir: str

    def __init__(self,
                 meta_data,
                 raw_data_dir = os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "composition")):
        self.meta_data = meta_data
        with open(meta_data) as fh:
            self.meta = json.load(fh)
        self.raw_data_dir = raw_data_dir

    def set_df(self, fname):
        self.df = pd.read_csv(os.path.join(self.raw_data_dir, fname),
                              index_col = "Time Stamp(s)",
                              skiprows = 14)
        self.df.index = pd.to_datetime(self.df.index)
        self.df.index.name = None
        self.df.drop("Records", axis=1, inplace=True)

    def get_interval_df(self, stamp):
        time_start = datetime.strptime(f"{stamp['date/yyyy-mm-dd']} {stamp['time_start/hh:mm:ss']}", "%Y-%m-%d %H:%M:%S")
        time_end = datetime.strptime(f"{stamp['date/yyyy-mm-dd']} {stamp['time_end/hh:mm:ss']}", "%Y-%m-%d %H:%M:%S")
        return self.df[(self.df.index >= time_start) & (self.df.index <= time_end)].copy()

    def get_mean_std_stunc(self, df, d_res=None):
        if d_res is None:
            d_res = {}
        d_res["DC voltage count"] = len([i for i in df["DC Voltage(VDC)"] if i is not np.nan])
        d_res["DC voltage (VDC)"] = df["DC Voltage(VDC)"].mean()
        d_res["DC voltage STD (VDC)"] = df["DC Voltage(VDC)"].std()
        d_res["DC voltage STUNC (VDC)"] = d_res["DC voltage STD (VDC)"] / (d_res["DC voltage count"])**0.5
        return d_res

    def __call__(self, derived_dir=None):
        date = None
        df_res = pd.DataFrame()
        for i, m in enumerate(self.meta["measurements"]):
            for j, stamp in enumerate(m["stamps"]):
                if stamp["date/yyyy-mm-dd"] != date:
                    date = stamp["date/yyyy-mm-dd"]
                    self.set_df(f"{date}.csv")
                df_j = self.get_interval_df(stamp)

                d_res = {}
                d_res["date/yyyy-mm-dd"] = stamp["date/yyyy-mm-dd"]
                d_res["time_start/hh:mm:ss"] = stamp["time_start/hh:mm:ss"]
                d_res["time_end/hh:mm:ss"] = stamp["time_end/hh:mm:ss"]
                # get means, stds, and uncertainty for given interval
                d_res = self.get_mean_std_stunc(df_j, d_res=d_res)
                _df = pd.DataFrame(d_res, index=[i+j])
                df_res = pd.concat([df_res, _df], axis=0)
        if derived_dir:
            fname = os.path.basename(self.meta_data)
            df_res.to_csv(os.path.join(derived_dir, f"{fname.split('.')[0]}_composition.csv"), sep=",", index=False)
        return df_res


def sbg_dry_run():
    meta_data = os.path.join(get_git_root(os.getcwd()), "data", "metadata", "sbg_dry_run.json")
    derived_data_dir = os.path.join(get_git_root(os.getcwd()), "data", "derived_data")
    mfc = MassFlowControllers(meta_data)
    sbg = SmartBiogas(meta_data)
    dfm = mfc()
    dfs = sbg()
    dfm.set_index(["date/yyyy-mm-dd", "time_start/hh:mm:ss", "time_end/hh:mm:ss"], inplace=True)
    dfs.set_index(["date/yyyy-mm-dd", "time_start/hh:mm:ss", "time_end/hh:mm:ss"], inplace=True)

    df = dfm.join(dfs, how="left")
    df["SBG err"] = df["Total vol. flow (ls/min)"] - df["Flow SBG (ls/min)"]
    df["SBG rel-err"] = (df["Total vol. flow (ls/min)"] - df["Flow SBG (ls/min)"]) / df["Total vol. flow (ls/min)"]
    df.to_csv(os.path.join(derived_data_dir, f"{os.path.basename(meta_data).split('.')[0]}.csv"), sep=",")


def sbg_humid_run():
    meta_data = os.path.join(get_git_root(os.getcwd()), "data", "metadata", "sbg_humid_run.json")
    aht_data_dir = os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "rH-temp")
    derived_data_dir = os.path.join(get_git_root(os.getcwd()), "data", "derived_data")
    mfc = MassFlowControllers(meta_data)
    sbg = SmartBiogas(meta_data)
    aht = AHT20(meta_data)
    dfm = mfc()
    dfs = sbg()
    dfh = aht()
    dfm.set_index(["date/yyyy-mm-dd", "time_start/hh:mm:ss", "time_end/hh:mm:ss"], inplace=True)
    dfs.set_index(["date/yyyy-mm-dd", "time_start/hh:mm:ss", "time_end/hh:mm:ss"], inplace=True)
    dfh.set_index(["date/yyyy-mm-dd", "time_start/hh:mm:ss", "time_end/hh:mm:ss"], inplace=True)

    df = dfm.join(dfs, how="left").join(dfh, how="left")
    df["SBG err"] = df["Total vol. flow (ls/min)"] - df["Flow SBG (ls/min)"]
    df["SBG rel-err"] = (df["Total vol. flow (ls/min)"] - df["Flow SBG (ls/min)"]) / df["Total vol. flow (ls/min)"]
    df.to_csv(os.path.join(derived_data_dir, f"{os.path.basename(meta_data).split('.')[0]}.csv"), sep=",")


def merge_dry_and_humid_sbg():
    derived_data_dir = os.path.join(get_git_root(os.getcwd()), "data", "derived_data")
    dry_res = "sbg_dry_run.csv"
    humid_res = "sbg_humid_run.csv"
    dfd = pd.read_csv(os.path.join(derived_data_dir, dry_res), sep=",")
    dfh = pd.read_csv(os.path.join(derived_data_dir, humid_res), sep=",")
    dfd.set_index(["date/yyyy-mm-dd", "time_start/hh:mm:ss", "time_end/hh:mm:ss"], inplace=True)
    dfh.set_index(["date/yyyy-mm-dd", "time_start/hh:mm:ss", "time_end/hh:mm:ss"], inplace=True)

    # populate humidity values with zeros for dry measurements with no humidity data acquisition
    dfd["AHT count"] = [0]*len(dfd.index)
    dfd["AHT Humidity (%rH)"] = [0]*len(dfd.index)
    dfres = pd.concat([dfd, dfh])
    dfres.to_csv(os.path.join(derived_data_dir, "sbg_dry_and_humid.csv"), sep=",")


def composition_sensor_calibration():
    meta_data = os.path.join(get_git_root(os.getcwd()), "data", "metadata", "all_sensors_humid_run.json")
    derived_data_dir = os.path.join(get_git_root(os.getcwd()), "data", "derived_data")
    mfc = MassFlowControllers(meta_data)
    aht = AHT20(meta_data)
    cs = CompositionSensor(meta_data)
    dfm = mfc()
    dfh = aht()
    dfc = cs()
    dfm.set_index(["date/yyyy-mm-dd", "time_start/hh:mm:ss", "time_end/hh:mm:ss"], inplace=True)
    dfh.set_index(["date/yyyy-mm-dd", "time_start/hh:mm:ss", "time_end/hh:mm:ss"], inplace=True)
    dfc.set_index(["date/yyyy-mm-dd", "time_start/hh:mm:ss", "time_end/hh:mm:ss"], inplace=True)
    df = dfm.join(dfh, how="left").join(dfc, how="left")
    df.to_csv(os.path.join(derived_data_dir, "composition_calibration.csv"), sep=",")


if __name__ == "__main__":
    # sensirion_test()
    # sbg_dry_run()
    # sbg_humid_run()
    # merge_dry_and_humid_sbg()
    composition_sensor_calibration()
