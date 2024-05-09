# This script uses ChemicalCase to load all the files in the samples subdirectory, and save
# all the filenames that have cos(SZA) between cose(82) and cos(98) to a new
# file that is a list of all the twilight_sample files.
import os
import numpy as np
from ChemicalCase import ChemicalCase

# get the list of files in the samples directory
files = os.listdir('samples')
twilight_files = []
for file in files:
    # load the file
    print('Processing', file)
    case = ChemicalCase('samples/' + file)
    # get the cosine of the solar zenith angle
    cos_sza = case.cosine_of_solar_zenith_angle
    # check if it is between cos(82) and cos(98)
    if cos_sza > np.cos(98*np.pi/180) and cos_sza < np.cos(82*np.pi/180):
        # print('cos(SZA) =', cos_sza)
        twilight_files.append(file)

# write the list of twilight files to a new text file
with open('twilight_files.txt', 'w') as f:
    for file in twilight_files:
        f.write(file + '\n')
