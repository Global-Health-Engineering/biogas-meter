import os
import json
import pandas as pd
from datetime import datetime


def combine_data(ground_truth_dir, smart_biogas_dir, mfc_meta, time_stamps):
    # get fluid names
    with open(mfc_meta) as fh:
        d_mfc = json.load(fh)
    flds = [value["FLD"] for value in d_mfc.values()]

    with open(time_stamps) as fh:
        d_time = json.load(fh)

    # make empty dataframe to populate with values
    df = pd.DataFrame()

    # iterate over measurement series
    for m in d_time["measurements"]:
        # read ground truth data
        df_gt = pd.read_csv(os.path.join(ground_truth_dir, f"{m['date/yyyy-mm-dd']}.csv"), index_col=[0])
        
        # convert indexes to datetime objects
        df_gt.index = pd.to_datetime(df_gt.index)
        
        # drop set point values
        for fld in flds:
            df_gt.drop([f"{fld} Flow setpoint (ln/min)"], axis=1, inplace=True)

        # read Smart Biogas data
        df_sb = pd.read_csv(os.path.join(smart_biogas_dir, f"{m['date/yyyy-mm-dd']}.csv"))
        
        # convert time column to datetime objects
        df_sb["Local Time"] = pd.to_datetime(df_sb["Local Time"])
        
        # remove time zone info from "Local Time" column in Smart Biogas dataframe
        df_sb["Local Time"] = df_sb["Local Time"].apply(lambda x: x.replace(tzinfo=None))
        
        # make "Local Time" column a new index
        df_sb.set_index("Local Time", drop=True, inplace=True)
        
        # drop nex index name
        df_sb.index.name = None
        
        # drop "UNIX Timestamp" column from Smart Biogas dataframe
        df_sb.drop("UNIX Timestamp", axis=1, inplace=True)
        
        # iterate over measurement timestamps
        for stamp in m["stamps"]:
            # convert timestamps to datetime objects
            time_start = datetime.strptime(f"{m['date/yyyy-mm-dd']} {stamp['time_start/hh:mm:ss']}", "%Y-%m-%d %H:%M:%S")
            time_end = datetime.strptime(f"{m['date/yyyy-mm-dd']} {stamp['time_end/hh:mm:ss']}", "%Y-%m-%d %H:%M:%S")
            
            # get ground truth data for given time range and copy it to temporary dataframe
            _df = df_gt[(df_gt.index >= time_start) & (df_gt.index <= time_end)].copy()

            # get data from Smart Biogas meter for given time range
            _df_sb = df_sb[(df_sb.index >= time_start) & (df_sb.index <= time_end)].copy()
            
            # write Smart Biogas measured pressure to tempoprary dataframe
            _df["Pressure SBG (Pa)"] = _df_sb["Pressure (Pa)"]
            
            # convert from Standard Liters Per Hour to Normal Liters Per Minute
            # 1 SLPM = 1 NLPM * (273.15 K / 293.15 K) * (14.696 psi / 14.504 psi)
            # write it to temporary dataframe
            _df["Flow SBG (ln/min)"] = _df_sb["Flow (lph)"] / 60 * 293.15 / 273.15 * 14.504 / 14.696

            # concatenate temporary dataframe to the main dataframe
            df = pd.concat([df, _df], axis=0)
    return df


def main():
    ground_truth_dir = "/Users/jtkaczuk/codes/biogas-meter/data/derived_data/ground_truth"
    smart_biogas_dir = "/Users/jtkaczuk/codes/biogas-meter/data/raw_data/smart_biogas"
    mfc_meta = "/Users/jtkaczuk/codes/biogas-meter/data/metadata/mfc_meta.json"
    exp_metadata = "/Users/jtkaczuk/codes/biogas-meter/data/metadata/exp1.json"
    df = combine_data(ground_truth_dir = ground_truth_dir,
                      smart_biogas_dir = smart_biogas_dir,
                      mfc_meta = mfc_meta,
                      time_stamps = exp_metadata)
    print(df)


if __name__ == "__main__":
    main()
