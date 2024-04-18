# src

Run `bronkhorst_processing.py` from the command line to postprocess raw data from Bronkhorst gas flow controllers by calling the following from the main repo directory.

Usage:
`
    process_raw_bronkhorst_measurements [-h] [--mass_props | --no-mass_props] [--molar_props | --no-molar_props] raw_csv meta_file derived_data_dir
`

Example:
`
    python src/bronkhorst_processing.py --mass_props --molar_props data/raw_data/test.csv data/metadata/test.json data/derived_data
`

Main function in `src/bronkhorst_processing.py` calculates the values below, averaged per second, for all gases in the system:  
- mean volumetric flow in ln/min
- standard deviation of volumetric flow of all fluids in ln/min
- one standard uncertainty of volumetric flow of all fluids ln/min
- mean temperature in degree C
- standard deviation of temperature in degree C
- one standard uncertainty of temperature in degree C
- mean pressure in bar(a)
- standard deviation of pressure in bar(a)
- one standard uncertainty of pressure in bar(a)
- mass flow rate for all flow streams and total mass flow in g/min
- mass fractions for all flow streams (dimensionless)
- molar flow rate for all flow streams and total molar flow in mol/min
- mole fractions for all flow streams (dimensionless)

Arguments:  
    `raw_csv` - input measurement file
    `meta_file` - file with metadata on measured gases in flow controllers
    `derived_data_dir` - directory where the processed file of the same name as 'raw_csv' will be saved
    `mass_props` - flag, will caculate mass-based measurements when set to True
    `molar_props` - flag, will caculate mole-based measurements when set to True

