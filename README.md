KPP Auto-reduction test bed
Derived from KPP-Compressor box model
Citation: Lin, H., et al. ... 
M.S. Long - Feb 10, 2022 -- initial release version

Updated May 2023 to include new input cases from GEOS-CF

# How to run
After building with `make`, the box model can be run using `./gckpp.exe`

Experiment 1 (the default)  runs a single initial condition in the format of the text files in the samples folder.  The local path and 
file name containing the initial condition should be specified in `filelist_exp1.txt`, for example `samples/LosAngeles_L1_20180702_1900.txt`.
The output of this is the species concentrations after 15 minutes, located in the file `final_concentrations_exp1.csv`





