# build_dirs.py
#
# Build directory structure for a new location to be
# analyzed using stochastic propagation modeling library
#
# Philip Blom (pblom@lanl.gov)

import os

from stochprop import eofs

# Set location and g2s file info
loc_dir = "US_RM"

prof_dir = "ncpag2s.r5ahhb/"
prof_prefix = "g2stxt_"
year_lims = [2010, 2016]

# build directories
if not os.path.isdir(loc_dir + "/"):
    print("Creating location directory...")
    os.system("mkdir " + loc_dir)

sub_dirs = ["coeffs", "eofs", "profs", "propagation", "samples"]
for dir in sub_dirs:
    if not os.path.isdir(loc_dir + "/" + dir):
        print('\t', "Creating " + dir + " directory...")
        os.system("mkdir " + loc_dir + "/" + dir)

for M in range(1, 13):
    if not os.path.isdir(loc_dir + "/profs/" + "{:02d}".format(M)):
        os.system("mkdir " + loc_dir + "/profs/" + "{:02d}".format(M))

seasons = ["fall", "spring", "summer", "winter"]
for season in seasons:
    if not os.path.isdir(loc_dir + "/samples/" + season):
        os.system("mkdir " + loc_dir + "/samples/" + season)

# move .dat profile files over
print('\n' + "Moving profiles into profs directory...")
for year in range(year_lims[0], year_lims[1] + 1):
    print('\t' + "Moving " + "{:04d}".format(year) + " profiles...")
    for M in range(1, 13):
        os.system("mv " + prof_dir + prof_prefix + "{:04d}".format(year) + "{:02d}".format(M) + "* " + loc_dir + "/profs/" + "{:02d}".format(M))

# Run QC
print('')
for M in range(1, 13):
    eofs.profs_check(loc_dir + "/profs/{:02d}/".format(M), prof_prefix + "*")
