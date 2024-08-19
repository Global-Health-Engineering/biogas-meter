import os
import git
import json
import ntpath
import numpy as np
import pandas as pd
import CoolProp as CP
from datetime import datetime
from get_git_root import get_git_root


class MFCaveraging(object):
    meta_file: str
    raw_csv: str
    
    def __init__(self,
                 meta_file,
                 raw_csv):
        with open(meta_file) as fh:
            self.meta = json.load(fh)
        self.csv_name = ntpath.basename(raw_csv)
        self.df_raw = pd.read_csv(raw_csv, sep=",")
        self.df_raw = self.remove_cols_with_no_metadata(self.df_raw)
        self.df_raw = self.df_raw.replace({np.nan: None})
        self.flds = tuple(val["FLD"] for val in self.meta.values())
        self.HEOS = tuple(CP.AbstractState("HEOS", fld) for fld in self.flds)

    def remove_cols_with_no_metadata(self, df):
        """ remove data from dataframe for which no metadata is passed """
        for col in df.columns:
            if (col != "TimeStamp") and (col.split(":")[0] not in self.meta.keys()):
                df.drop(col, axis=1, inplace=True)
        return df

    def get_rhomass(self, idx, p_bar, T_degC):
        """ return fluid density in g/L """
        self.HEOS[idx].update(CP.PT_INPUTS, p_bar*1e5, T_degC+273.15)
        return self.HEOS[idx].rhomass()

    def get_rhomolar(self, idx, p_bar, T_degC):
        """ return fluid density in mol/L """
        self.HEOS[idx].update(CP.PT_INPUTS, p_bar*1e5, T_degC+273.15)
        return self.HEOS[idx].rhomolar() / 1e3

    def get_mean_measurements(self, interval):
        # make empty dataframe
        df = pd.DataFrame()
        
        # copy raw data
        _df = self.df_raw.copy()

        # convert strings to datetime objects
        _df["FullTimeStamp"] = pd.to_datetime(_df["TimeStamp"])
        
        # round datetime to interval in new column
        _df["TimeStamp"] = _df["FullTimeStamp"].dt.round(interval)
        
        # remove time zone info in new datetime column
        _df["TimeStamp"] = _df["TimeStamp"].apply(lambda x: x.replace(tzinfo=None))
        
        # drop full time information
        _df.drop("FullTimeStamp", axis=1, inplace=True)

        for i, group in _df.groupby("TimeStamp"):
            d = {}

            # number of measurements per stream
            for key in self.meta.keys():
                d[f"{key}: ground truth count"] = len([j for j in group[f"{key}: Flow (ln/min)"] if i is not None])

            for col in group.columns:
                if col != "TimeStamp":
                    # mean values per interval
                    d[col] = group[col].mean()
                    
                    if "setpoint" not in col.lower():
                        # standard deviation (per interval)
                        d[f"{col} STD"] = group[col].std()
        
                        # standard uncertainty (per interval)
                        d[f"{col} STUNC"] = d[f"{col} STD"] / (d[f"{col.split(': ')[0]}: ground truth count"])**0.5
        
            df = pd.concat([df, pd.DataFrame(d, index=[i])])
        
        # total volumetric flow
        df["Total vol. flow (ln/min)"] = sum([df[f"{key}: Flow (ln/min)"] for key in self.meta.keys()])
        
        # relative combined uncertainty of total volumetric flow
        df["Total vol. flow (ln/min) COMB-UNC"] = 0
        df.loc[df["Total vol. flow (ln/min)"] != 0, "Total vol. flow (ln/min) COMB-UNC"] = (
            sum([df[f"{key}: Flow (ln/min) STUNC"]**2 for key in self.meta.keys()])**0.5)
        
        for key in self.meta.keys():
            # calculate volumetric fraction
            df[f"{key}: Vol. fraction"] = 0
            df.loc[df["Total vol. flow (ln/min)"] != 0, f"{key}: Vol. fraction"] = df[f"{key}: Flow (ln/min)"] / df["Total vol. flow (ln/min)"]
        
            # calculate relative combined uncertainty of volumatric fraction
            df[f"{key}: Vol. fraction COMB-UNC"] = 0
            df.loc[df[f"{key}: Flow (ln/min)"] != 0, f"{key}: Vol. fraction COMB-UNC"] = (
                (df[f"{key}: Flow (ln/min) STUNC"] / df[f"{key}: Flow (ln/min)"])**2
                + (df["Total vol. flow (ln/min) COMB-UNC"] / df["Total vol. flow (ln/min)"])**2)**0.5
        
        return df

    def get_per_sec_mean_measurements(self):
        return self.get_mean_measurements("1s")

    def get_per_min_mean_measurements(self):
        return self.get_mean_measurements("1min")

    def get_mass_fractions(self, df):
        for key, value in self.meta.items():
            rhomass = df.apply(lambda x: self.get_rhomass(idx = value["ID"],
                                                          p_bar = x[f"{key}: Pressure outlet (bar(a))"],
                                                          T_degC = x[f"{key}: Temperature (°C)"]), axis=1)
            df[f"{key}: Mass flow (g/min)"] = df[f"{key}: Flow (ln/min)"] * rhomass
        df["Total mass flow (g/min)"] = sum([df[f"{key}: Mass flow (g/min)"] for key in self.meta.keys()])
        for key in self.meta.keys():
            
            # calculate mass fraction
            df[f"{key}: Mass fraction"] = 0
            df.loc[df["Total mass flow (g/min)"] != 0, f"{key}: Mass fraction"] = df[f"{key}: Mass flow (g/min)"] / df["Total mass flow (g/min)"]
        
        return df

    def get_molar_fractions(self, df):
        for key, value in self.meta.items():
            rhomolar = df.apply(lambda x: self.get_rhomolar(idx = value["ID"],
                                                            p_bar = x[f"{key}: Pressure outlet (bar(a))"],
                                                            T_degC = x[f"{key}: Temperature (°C)"]), axis=1)
            df[f"{key}: Molar flow (mol/min)"] = df[f"{key}: Flow (ln/min)"] * rhomolar
        df["Total molar flow (mol/min)"] = sum([df[f"{key}: Molar flow (mol/min)"] for key in self.meta.keys()])
        
        for key in self.meta.keys():
        
            # calculate molar fraction
            df[f"{key}: Mole fraction"] = 0
            df.loc[df["Total molar flow (mol/min)"] != 0, f"{key}: Mole fraction"] = df[f"{key}: Molar flow (mol/min)"] / df["Total molar flow (mol/min)"]
        
        return df

    def __call__(self,
                 sec_avg = False,
                 mass_props = True,
                 molar_props = True,
                 derived_data_dir = None):
        df = self.get_per_sec_mean_measurements() if sec_avg else self.get_per_min_mean_measurements()
        if mass_props:
            df = self.get_mass_fractions(df)
        if molar_props:
            df = self.get_molar_fractions(df)

        unique_flds = True if len(set(self.flds)) == len(self.flds) else False
        # create column mapper to swap mass flow controller S/N with fluid name
        name_map = {}
        for key, val in self.meta.items():
            for col in df.columns:
                if key in col:
                    new_col_name = val["FLD"] if unique_flds else f'{val["FLD"]}_{val["ID"]}'
                    name_map[col] = f"{new_col_name} {col.lstrip(f'{key}:').strip(' ')}"
        # rename columns (swap mass flow controller S/N with fluid name)
        df = df.rename(columns=name_map)

        if derived_data_dir:
            df.to_csv(os.path.join(derived_data_dir, f"mfc_{self.csv_name}"), sep=",")

        return df


class SensirionAveraging(object):
    edf_file: str

    def __init__(self, edf_file):
        self.edf_file = edf_file
        self.raw_df = pd.read_csv(edf_file,
                                  engine = "python",
                                  encoding = "utf-8",
                                  sep = r"\t",
                                  skiprows = 10)
        self.raw_df = self.set_col_names(self.raw_df)

    def set_col_names(self, df):
        mapper = {}
        for col in list(df.columns):
            if col.split("_")[0] == "F":
                mapper[col] = "Vol. flow (ls/min)"
            elif col.split("_")[0] == "T":
                mapper[col] = "T (deg C)"
        df.rename(mapper, axis=1, inplace=True)
        # convert from Standard Liters Per Minute to Normal Liters Per Minute
        # 1 SLPM = 1 NLPM * (273.15 K / 293.15 K) * (14.696 psi / 14.504 psi)
        df["Vol. flow (ln/min)"] = df["Vol. flow (ls/min)"] * 293.15 / 273.15 * 14.504 / 14.696
        df.drop("Vol. flow (ls/min)", axis=1, inplace=True)
        return df

    def get_mean_measurements(self, interval):
        # make empty dataframe
        df = pd.DataFrame()

        # make a copy of raw data with renamed columns
        _df = self.raw_df.copy()

        # swap nan with None
        _df.replace({np.nan: None})

        # drop UTC time info
        _df.drop("Epoch_UTC", axis=1, inplace=True)

        # convert strings to datetime objects
        _df.Local_Date_Time = pd.to_datetime(_df.Local_Date_Time)

        # drop time zone information
        _df.Local_Date_Time = _df.Local_Date_Time.apply(lambda x: x.replace(tzinfo=None))

        # round to inverval
        _df.Local_Date_Time = _df.Local_Date_Time.dt.round(interval)

        # set index to cleaned time
        _df.index = _df.Local_Date_Time

        # drop full time information
        _df.drop("Local_Date_Time", axis=1, inplace=True)

        # remove index name
        _df.index.name = None

        # caculate means and statistics per interval
        for i, group in _df.groupby(_df.index):
            d = {}
            d["count"] = len([j for j in group["Vol. flow (ln/min)"] if j is not None])
            for col in _df.columns:
                d[col] = group[col].mean()
                d[f"{col} STD"] = group[col].std()
                d[f"{col} STUNC"] = d[f"{col} STD"] / (d["count"])**0.5
            df = pd.concat([df, pd.DataFrame(d, index=[i])])

        return df

    def get_per_sec_mean_measurements(self):
        return self.get_mean_measurements("1s")

    def get_per_min_mean_measurements(self):
        return self.get_mean_measurements("1min")

    def __call__(self,
                 sec_avg = False,
                 derived_data_dir = None):
        df = self.get_per_sec_mean_measurements() if sec_avg else self.get_per_min_mean_measurements()

        if derived_data_dir:
            df.to_csv(os.path.join(derived_data_dir,
                                   f"sensirion_{ntpath.basename(self.edf_file).rstrip('.edf')}.csv"))


def main():
    fname = "2024-08-15"
    metadata_name = "mfc_air_meta.json"
    raw_csv = os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "mfc", f"{fname}.csv")
    meta = os.path.join(get_git_root(os.getcwd()), "data", "metadata", metadata_name)
    derived_data_dir = os.path.join(get_git_root(os.getcwd()), "data", "derived_data")

    sec_avg = True
    mass_props = False
    molar_props = False

    mfc = MFCaveraging(meta_file = meta,
                       raw_csv = raw_csv)
    mfc(sec_avg = sec_avg,
        mass_props = mass_props,
        molar_props = molar_props,
        derived_data_dir = derived_data_dir)

    raw_data = os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "sensirion", f"{fname}.edf")
    
    sa = SensirionAveraging(raw_data)
    sa(sec_avg = sec_avg,
       derived_data_dir = derived_data_dir)


if __name__ == "__main__":
    main()
