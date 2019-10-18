# build_dirs.py
#
# Build directory structure for a new location to be
# analyzed using stochastic propagation modeling library
#
# Philip Blom (pblom@lanl.gov)

import os

from stochprop import eofs

if __name__ == '__main__':
    # ######################### #
    #   Define parameters for   #
    #      organizing data      #
    # ######################### #
    prof_dir = "dir/of/g2s/"

    prof_prefix = "g2stxt_"
    year_lims = [2010, 2016]
    
    # ######################### #
    #    Build directories to   #
    #      organize results     #
    # ######################### #
    sub_dirs = ["coeffs", "eofs", "profs", "propagation", "samples"]
    seasons = ["fall", "spring", "summer", "winter"]
    
    print("Creating directories...")
    for dir in sub_dirs:
        if not os.path.isdir(dir):
            os.system("mkdir " + dir)
    
    for M in range(1, 13):
        if not os.path.isdir("profs/" + "{:02d}".format(M)):
            os.system("mkdir profs/" + "{:02d}".format(M))

    for season in seasons:
        if not os.path.isdir("samples/" + season):
            os.system("mkdir samples/" + season)
    
    print('\n' + "Moving profiles into profs directory...")
    for year in range(year_lims[0], year_lims[1] + 1):
        print('\t' + "Moving " + "{:04d}".format(year) + " profiles...")
        for M in range(1, 13):
            os.system("mv " + prof_dir + prof_prefix + "{:04d}".format(year) + "{:02d}".format(M) + "* " + loc_dir + "/profs/" + "{:02d}".format(M))

    # Run QC
    for M in range(1, 13):
        eofs.profs_check("profs/{:02d}/".format(M), prof_prefix + "*")
