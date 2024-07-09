import os
import json
import ntpath
import argparse
import numpy as np
import pandas as pd
from datetime import datetime
import CoolProp as CP


class Bronhorst_processing(object):
    meta_file: str
    raw_csv: str
    k : float
    
    def __init__(self,
                 meta_file,
                 raw_csv,
                 k = 1.96):
        self.csv_name = ntpath.basename(raw_csv)
        self.df_raw = pd.read_csv(raw_csv, sep=",")
        with open(meta_file) as fh:
            self.meta = json.load(fh)
        self.flds = tuple(val["FLD"] for val in self.meta.values())
        self.HEOS = tuple(CP.AbstractState("HEOS", fld) for fld in self.flds)
        self.set_coverage_factor(k)

    def set_coverage_factor(self, k):
        """ setter for coverage factor for expanded uncertainty calculations """
        self.k = k

    def get_rhomass(self, idx, p_bar, T_degC):
        """ return fluid density in g/L """
        self.HEOS[idx].update(CP.PT_INPUTS, p_bar*1e5, T_degC+273.15)
        return self.HEOS[idx].rhomass()

    def get_rhomolar(self, idx, p_bar, T_degC):
        """ return fluid density in mol/L """
        self.HEOS[idx].update(CP.PT_INPUTS, p_bar*1e5, T_degC+273.15)
        return self.HEOS[idx].rhomolar() / 1e3

    def get_per_second_avg_measurements(self):
        # make empty dataframe
        df = pd.DataFrame()
        # copy raw data
        _df = self.df_raw.copy()

        # convert strings to datetime objects
        _df["FullTimeStamp"] = pd.to_datetime(_df["TimeStamp"])
        # round datetime to 1 sec in new column
        _df["TimeStamp"] = _df["FullTimeStamp"].dt.round("1s")
        # remove time zone info in new datetime column
        _df["TimeStamp"] = _df["TimeStamp"].apply(lambda x: x.replace(tzinfo=None))
        # drop full time information
        _df.drop("FullTimeStamp", axis=1, inplace=True)

        for i, group in _df.groupby("TimeStamp"):
            d = {}
            for col in group.columns:
                if col != "TimeStamp":
                    d[col] = group[col].mean()
                    if "setpoint" not in col.lower():
                        # calculate standard deviation (per second)
                        d[f"{col} STD"] = group[col].std()
                        # calculate standard uncertainty (per second)
                        d[f"{col} STUNC"] = d[f"{col} STD"] / len(group)**0.5
            df = pd.concat([df, pd.DataFrame(d, index=[i])])
        # calculate total volumetric flow
        df["Total vol. flow (ln/min)"] = sum([df[f"{key}: Flow (ln/min)"] for key in self.meta.keys()])
        # calculate relative combined uncertainty of total volumetric flow
        df["Total vol. flow (ln/min) COMB-UNC"] = 0
        df.loc[df["Total vol. flow (ln/min)"] != 0, "Total vol. flow (ln/min) COMB-UNC"] = (
            sum([df[f"{key}: Flow (ln/min) STUNC"]**2 for key in self.meta.keys()])**0.5)
        for key in self.meta.keys():
            # calculate volumetric fraction
            df[f"{key}: Vol. fraction"] = 0
            df.loc[df["Total vol. flow (ln/min)"] != 0, f"{key}: Vol. fraction"] = df[f"{key}: Flow (ln/min)"] / df["Total vol. flow (ln/min)"]
            # calculate relative combined uncertainty of volumatric fraction
            df[f"{key}: Vol. fraction EXP-COMB-UNC"] = 0
            df.loc[df[f"{key}: Flow (ln/min)"] != 0, f"{key}: Vol. fraction EXP-COMB-UNC"] = self.k * (
                (df[f"{key}: Flow (ln/min) STUNC"] / df[f"{key}: Flow (ln/min)"])**2
                + (df["Total vol. flow (ln/min) COMB-UNC"] / df["Total vol. flow (ln/min)"])**2)**0.5
        return df

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

    def __call__(self, mass_props=True, molar_props=True, derived_data_dir=None):
        df = self.get_per_second_avg_measurements()
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
                    name_map[col] = f"{new_col_name} {col.lstrip(f'{key}:')}"
        # rename columns (swap mass flow controller S/N with fluid name)
        df = df.rename(columns=name_map)

        if derived_data_dir:
            df.to_csv(os.path.join(derived_data_dir, self.csv_name), sep=",")

        return df


def test():
    bp = Bronhorst_processing(
        meta_file="data/metadata/test.json",
        raw_csv="data/raw_data/test.csv")
    
    bp(mass_props=True,
       molar_props=True,
       derived_data_dir="data/derived_data")


def main():
    info = """
    Calculate flow information from raw Bronkhorst mass flow controller data.
    Arguments:
        'raw_csv' - input measurement file
        'meta_file' - file with metadata on measured gases in flow controllers
        'derived_data_dir' - directory where the processed file of the same name as 'raw_csv' will be saved
        'mass_props' - flag, will caculate mass-based measurements when set to True
        'molar_props' - flag, will caculate mole-based measurements when set to True
    """
    parser = argparse.ArgumentParser(prog="process_raw_bronkhorst_measurements", description=info)
    parser.add_argument("raw_csv")
    parser.add_argument("meta_file")
    parser.add_argument("derived_data_dir")
    parser.add_argument('--mass_props', action=argparse.BooleanOptionalAction)
    parser.add_argument('--molar_props', action=argparse.BooleanOptionalAction)
    args = parser.parse_args()

    bp = Bronhorst_processing(meta_file = args.meta_file,
                              raw_csv = args.raw_csv)
    bp(mass_props = args.mass_props,
       molar_props = args.molar_props,
       derived_data_dir = args.derived_data_dir)

    print(f"Done. File {os.path.join(args.derived_data_dir, ntpath.basename(args.raw_csv))} created.")


if __name__ == "__main__":
    main()
    # test()
