import os
import git
import json
import ntpath
import argparse
import numpy as np
import pandas as pd
from datetime import datetime
from get_git_root import get_git_root
from sensors_all import MFCPerInterval
from sensors_all import SensirionPerInterval
from sensors_all import HumidityPerInterval


class Analyze_MFC_Sen_AHT(object):
    mfc_dir: str
    sen_dir: str
    aht_dir: str
    metaData: str

    def __init__(self, mfc_dir, sen_dir, aht_dir, metaData):
        self.mfc_dir = mfc_dir
        self.sen_dir = sen_dir
        self.aht_dir = aht_dir
        self.metaData = metaData
        with open(metaData) as fh:
            self.meta = json.load(fh)

    def get_mfc_df(self, fname):
        df = pd.read_csv(os.path.join(self.mfc_dir, fname))
        # convert strings to datetime objects
        df["TimeStamp"] = pd.to_datetime(df["TimeStamp"])
        # remove time zone info
        df["TimeStamp"] = df["TimeStamp"].apply(lambda x: x.replace(tzinfo=None))
        # make "TimeStamp" column new index column
        df.set_index("TimeStamp", drop=True, inplace=True)
        # drop index name
        df.index.name = None
        return df

    def get_sen_df(self, edf_fname):
        df = pd.read_csv(os.path.join(self.sen_dir, edf_fname),
                         engine = "python",
                         encoding = "utf-8",
                         sep = r"\t",
                         skiprows = 10)
        # replace nan with None
        df.replace({np.nan: None})
        # drop UTC timestamps
        df.drop("Epoch_UTC", axis=1, inplace=True)
        # convert timestamps to datetime objects
        df.Local_Date_Time = pd.to_datetime(df.Local_Date_Time)
        # remove time zone information
        df.Local_Date_Time = df.Local_Date_Time.apply(lambda x: x.replace(tzinfo=None))
        # set time to index
        df.index = df.Local_Date_Time
        # drop old time column
        df.drop("Local_Date_Time", axis=1, inplace=True)
        # remove index name
        df.index.name = None
        # rename columns
        mapper = {}
        for col in list(df.columns):
            if col.split("_")[0] == "F":
                mapper[col] = "Sensirion Vol. flow (ls/min)"
            elif col.split("_")[0] == "T":
                mapper[col] = "Sensirion T (deg C)"
        df.rename(mapper, axis=1, inplace=True)
        return df

    def get_aht20_df(self, fname):
        df = pd.read_csv(os.path.join(self.aht_dir, fname),
                         index_col="Datetime")
        # remove index name
        df.index.name = None
        # convert timestamps to datetime objects
        df.index = pd.to_datetime(df.index)
        # rename columns
        mapper = {"Temp (deg C)": "AHT temp (deg C)",
                  "Humidity (rH)": "AHT humidity (%rH)"}
        df.rename(mapper, axis=1, inplace=True)
        return df

    def get_interval_df(self, df, stamp):
        time_start = datetime.strptime(f"{stamp['date/yyyy-mm-dd']} {stamp['time_start/hh:mm:ss']}", "%Y-%m-%d %H:%M:%S")
        time_end = datetime.strptime(f"{stamp['date/yyyy-mm-dd']} {stamp['time_end/hh:mm:ss']}", "%Y-%m-%d %H:%M:%S")
        return df[(df.index >= time_start) & (df.index <= time_end)].copy()

    def convert_lspm_to_lnpm(self, V_lspm, T_degC, p_bara):
        # convert from Standard Liters Per Minute to Normal Liters Per Minute
        return V_lspm * 273.15 / (273.15 + T_degC) * p_bara / 1.01325

    def __call__(self, mass_props=False, molar_props=False, derived_dir=None):
        date = None
        df = pd.DataFrame()
        for m in self.meta["measurements"]:
            for i, stamp in enumerate(m["stamps"]):
                if stamp["date/yyyy-mm-dd"] != date:
                    date = stamp["date/yyyy-mm-dd"]
                    df_mfc = self.get_mfc_df(f"{stamp['date/yyyy-mm-dd']}.csv")
                    df_sen = self.get_sen_df(f"{stamp['date/yyyy-mm-dd']}.edf")
                    df_aht = self.get_aht20_df(f"{stamp['date/yyyy-mm-dd']}.csv")
                df_mfc_i = self.get_interval_df(df_mfc, stamp)
                df_sen_i = self.get_interval_df(df_sen, stamp)
                df_aht_i = self.get_interval_df(df_aht, stamp)
                    
                mfc = MFCPerInterval(df=df_mfc_i, metaData=self.metaData)
                sen = SensirionPerInterval(df=df_sen_i)
                aht = HumidityPerInterval(df=df_aht_i)

                d = {}
                d["date/yyyy-mm-dd"] = stamp["date/yyyy-mm-dd"]
                d["time_start/hh:mm:ss"] = stamp["time_start/hh:mm:ss"]
                d["time_end/hh:mm:ss"] = stamp["time_end/hh:mm:ss"]
                # get ground truth means, stds, uncertainties, etc. for given interval
                d.update(mfc(mass_props=mass_props, molar_props=molar_props))
                d.update(sen())
                d.update(aht())

                _df = pd.DataFrame(d, index=[i])
                df = pd.concat([df, _df], axis=0)

        df["Sensirion Vol. flow (ln/min)"] = self.convert_lspm_to_lnpm(V_lspm = df["Sensirion Vol. flow (ls/min)"],
                                                                       T_degC = df["Sensirion T (deg C)"],
                                                                       p_bara = df["Air Pressure outlet (bar(a))"])

        if derived_dir:
            fname = os.path.basename(self.metaData)
            df.to_csv(os.path.join(derived_dir, f"{fname.split('.')[0]}.csv"), sep=",", index=False)
        return df


def main():
    metaData = os.path.join(get_git_root(os.getcwd()), "data", "metadata", "mfc_sensirion_aht20.json")
    mfc_dir = os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "mfc")
    sen_dir = os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "sensirion")
    aht_dir = os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "rH-temp")
    derived_dir = os.path.join(get_git_root(os.getcwd()), "data", "derived_data")

    amsa = Analyze_MFC_Sen_AHT(mfc_dir = mfc_dir,
                               sen_dir = sen_dir,
                               aht_dir = aht_dir,
                               metaData = metaData)

    amsa(derived_dir = derived_dir)


if __name__ == "__main__":
    main()
