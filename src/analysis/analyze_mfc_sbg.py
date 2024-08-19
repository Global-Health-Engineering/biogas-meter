import os
import git
import json
import ntpath
import argparse
import numpy as np
import pandas as pd
from datetime import datetime
from sensors_all import MFCPerInterval
from sensors_all import SBGPerInterval


class Analyze_MFC_SBG(object):
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
        assert os.path.isfile(os.path.join(path, fname))
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
                    date = stamp["date/yyyy-mm-dd"]
                    df_mfc = self.get_mfc_df(f"{stamp['date/yyyy-mm-dd']}.csv")
                    df_sbg = self.get_sbg_df(f"{stamp['date/yyyy-mm-dd']}.csv")
                df_mfc_i = self.get_interval_df(df_mfc, stamp)
                df_sbg_i = self.get_interval_df(df_sbg, stamp)

                mfcperinterval = MFCPerInterval(df=df_mfc_i, metaData=self.metaData)
                sbgperinterval = SBGPerInterval(df=df_sbg_i)

                d = {}
                d["date/yyyy-mm-dd"] = stamp["date/yyyy-mm-dd"]
                d["time_start/hh:mm:ss"] = stamp["time_start/hh:mm:ss"]
                d["time_end/hh:mm:ss"] = stamp["time_end/hh:mm:ss"]
                # get ground truth means, stds, uncertainties, etc. for given interval
                d.update(mfcperinterval(mass_props=mass_props, molar_props=molar_props))
                # get Smart Biogas means and stds
                d.update(sbgperinterval())
                # calculate relative error of Smart Biogas measurements
                d["SBG rel-err"] = (d["Total vol. flow (ln/min)"] - d["Flow SBG (ln/min)"]) / d["Total vol. flow (ln/min)"]

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
    metaData = os.path.join(get_git_root(os.getcwd()), "data", "metadata", "full_sgb_dry_run.json")
    mfc_dir = os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "mfc")
    sbg_dir = os.path.join(get_git_root(os.getcwd()), "data", "raw_data", "sbg")
    derived_dir = os.path.join(get_git_root(os.getcwd()), "data", "derived_data")

    an = Analyze_MFC_SBG(mfc_dir = mfc_dir,
                         sbg_dir = sbg_dir,
                         metaData = metaData)
    an(mass_props = False,
       molar_props = False,
       derived_dir = derived_dir)


if __name__ == "__main__":
    main()
