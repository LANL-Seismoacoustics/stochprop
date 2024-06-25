# propagation.py
#
# Run infraga-accel-3d or NCPAprop's modess to generate predictions to
# use for constructing stochastic propagation models
#
# Philip Blom (pblom@lanl.gov)

import os
import fnmatch
import warnings
import pickle
import itertools
import subprocess
import glob
import pkg_resources
import tempfile

import configparser as cnfg
import numpy as np

from importlib.util import find_spec

from datetime import datetime
from pyproj import Geod

from scipy.integrate import simps
from scipy.interpolate import interp1d, RectBivariateSpline
from scipy.stats import norm, gaussian_kde
from scipy.signal import savgol_filter
from scipy.special import gamma

import matplotlib.cm as cm
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colorbar import ColorbarBase

from mpl_toolkits.axes_grid1 import make_axes_locatable

import cartopy.crs as crs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

sph_proj = Geod(ellps='sphere')

map_proj = crs.PlateCarree()
resol = '100m'  # use data at this scale (not working at the moment)


# ############################## #
#  Plotting the Effective Sound  #
#   Speed Ratio for Seasonality  #
# ############################## #

def _plot_ess_ratio(A, z0, datetimes, grnd_index, results_path, include_ns, show_title, output_header):

    plt.rcParams.update({'font.size': 16})

    print('\n' + "Computing effective sound speed ratio for each day-of-year...")
    f1, ax1 = plt.subplots(2, figsize=(12, 6), sharex=True)
    ax1[0].set_xlim(0, 52)
    ax1[1].set_ylim(0, 100)
    ax1[0].set_xticks(range(0, 52, 8))
    ax1[0].set_ylabel("Peak ESS Ratio")
    ax1[1].set_xlabel("Week of Year")
    ax1[1].set_ylabel("Altitude [km]")
    if show_title:
        ax1[0].set_title("Effective Sound Speed (ESS) Ratio Analysis \n Eastward (blue), Westward (red)")

    if include_ns:
        f2, ax2 = plt.subplots(2, figsize=(12, 6), sharex=True)
        ax2[0].set_xlim(0, 52)
        ax2[1].set_ylim(0, 100)
        ax2[0].set_xticks(range(0, 52, 8))
        ax2[0].set_ylabel("Peak ESS Ratio")
        ax2[1].set_xlabel("Week of Year")
        ax2[1].set_ylabel("Altitude [km]")
        if show_title:
            ax2[0].set_title("Effective Sound Speed (ESS) Ratio Analysis \n Northward (purple), Southward (orange)")

    eff_sndspd_pk = np.empty((4, 365))
    for j, yday in enumerate(["{:03d}".format(m + 1) for m in range(365)]):

        yday_mask = [abs(j - dt_n.astype(datetime).timetuple().tm_yday) < 3 for dt_n in datetimes]
        eff_sndspd_ratio = np.empty((4, len(datetimes[yday_mask]), len(z0)))
        for n, An in enumerate(A[yday_mask]):
            u = An[1 * len(z0):2 * len(z0)]
            v = An[2 * len(z0):3 * len(z0)]
            d = An[3 * len(z0):4 * len(z0)]
            p = An[4 * len(z0):5 * len(z0)]

            c_eff = np.sqrt(0.14 * p / d) + np.sin(np.radians(90.0)) * u + np.cos(np.radians(90.0)) * v
            eff_sndspd_ratio[0][n] = c_eff / c_eff[grnd_index]
            
            c_eff = np.sqrt(0.14 * p / d) + np.sin(np.radians(-90.0)) * u + np.cos(np.radians(-90.0)) * v
            eff_sndspd_ratio[1][n] = c_eff / c_eff[grnd_index]

            if include_ns:
                c_eff = np.sqrt(0.14 * p / d) + np.sin(np.radians(0.0)) * u + np.cos(np.radians(0.0)) * v
                eff_sndspd_ratio[2][n] = c_eff / c_eff[grnd_index]
                
                c_eff = np.sqrt(0.14 * p / d) + np.sin(np.radians(180.0)) * u + np.cos(np.radians(180.0)) * v
                eff_sndspd_ratio[3][n] = c_eff / c_eff[grnd_index]

        plot_mask = np.logical_and(z0 < 100.0, np.mean(eff_sndspd_ratio[0], axis=0) > 1.0)
        if np.sum(plot_mask) > 1:
            im1 = ax1[1].scatter([float(yday) / 7.0] * len(z0[plot_mask]), z0[plot_mask], c=np.mean(eff_sndspd_ratio[0], axis=0)[plot_mask], s=1.0, cmap=cm.seismic_r, vmin=0.9, vmax=1.1)

        plot_mask = np.logical_and(z0 < 100.0, np.mean(eff_sndspd_ratio[1], axis=0) > 1.0)
        if np.sum(plot_mask) > 1:
            ax1[1].scatter([float(yday) / 7.0] * len(z0[plot_mask]), z0[plot_mask], c=np.mean(eff_sndspd_ratio[1], axis=0)[plot_mask], s=1.0, cmap=cm.seismic, vmin=0.9, vmax=1.1)

        for k in range(2):
            eff_sndspd_pk[k][j] = np.max(np.mean(eff_sndspd_ratio[k], axis=0)[np.logical_and(35.0 <= z0, z0 <= 70.0)])    

        if include_ns:
            plot_mask = np.logical_and(z0 < 100.0, np.mean(eff_sndspd_ratio[2], axis=0) > 1.0)
            if np.sum(plot_mask) > 1:
                ax2[1].scatter([float(yday) / 7.0] * len(z0[plot_mask]), z0[plot_mask], c=np.mean(eff_sndspd_ratio[2], axis=0)[plot_mask], s=1.0, cmap=cm.PuOr, vmin=0.9, vmax=1.1)

            plot_mask = np.logical_and(z0 < 100.0, np.mean(eff_sndspd_ratio[3], axis=0) > 1.0)
            if np.sum(plot_mask) > 1:
                im2 = ax2[1].scatter([float(yday) / 7.0] * len(z0[plot_mask]), z0[plot_mask], c=np.mean(eff_sndspd_ratio[3], axis=0)[plot_mask], s=1.0, cmap=cm.PuOr_r, vmin=0.9, vmax=1.1)

            for k in range(2, 4):
                eff_sndspd_pk[k][j] = np.max(np.mean(eff_sndspd_ratio[k], axis=0)[np.logical_and(35.0 <= z0, z0 <= 70.0)])    

    cbar = plt.colorbar(im1, ax=ax1.ravel().tolist(), location="bottom", shrink=0.8, aspect=40)
    cbar.set_ticks(np.arange(0.9, 1.1, 0.1))
    cbar.set_ticklabels(['1.1 (west)', '1.0', '1.1 (east)'])
    cbar.set_label("Effective Sound Speed (ESS) Ratio")

    if include_ns:
        cbar = plt.colorbar(im2, ax=ax2.ravel().tolist(), location="bottom", shrink=0.8, aspect=40)
        cbar.set_ticks(np.arange(0.9, 1.1, 0.1))
        cbar.set_ticklabels(['1.1 (south)', '1.0', '1.1 (north)'])

    # Output summary of ess_ratio unity crossings
    print('\n' + "Eastward waveguide changes...")
    for j in range(364):
        if (eff_sndspd_pk[0][j] - 1.0) * (eff_sndspd_pk[0][j + 1] - 1) <= 0.0:
            if eff_sndspd_pk[0][j] > eff_sndspd_pk[0][j + 1]:
                print('\t' + "Waveguide dissipates:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")")
            else:
                print('\t' + "Waveguide forms:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")")

    print('\n' + "Westward waveguide changes...")
    for j in range(364):
        if (eff_sndspd_pk[1][j] - 1.0) * (eff_sndspd_pk[1][j + 1] - 1) <= 0.0:
            if eff_sndspd_pk[1][j] > eff_sndspd_pk[1][j + 1]:
                print('\t' + "Waveguide dissipates:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")")
            else:
                print('\t' + "Waveguide forms:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")")
    if include_ns:
        print('\n' + "Northward waveguide changes...")
        for j in range(364):
            if (eff_sndspd_pk[2][j] - 1.0) * (eff_sndspd_pk[2][j + 1] - 1) <= 0.0:
                if eff_sndspd_pk[2][j] > eff_sndspd_pk[2][j + 1]:
                    print('\t' + "Waveguide dissipates:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")")
                else:
                    print('\t' + "Waveguide forms:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")")

        print('\n' + "Southward waveguide changes...")
        for j in range(364):
            if (eff_sndspd_pk[3][j] - 1.0) * (eff_sndspd_pk[3][j + 1] - 1) <= 0.0:
                if eff_sndspd_pk[3][j] > eff_sndspd_pk[3][j + 1]:
                    print('\t' + "Waveguide dissipates:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")")
                else:
                    print('\t' + "Waveguide forms:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")")
    print('')

    ax1[0].plot(np.arange(365.0) / 7.0, eff_sndspd_pk[0], '-b', linewidth=2.5)
    ax1[0].plot(np.arange(365.0) / 7.0, eff_sndspd_pk[1], '--r', linewidth=2.5)
    ax1[0].axhline(1.0, color='k', linestyle='dashed')

    if include_ns:
        ax2[0].plot(np.arange(365.0) / 7.0, eff_sndspd_pk[2], color='purple', linewidth=2.5)
        ax2[0].plot(np.arange(365.0) / 7.0, eff_sndspd_pk[3], color='orange', linewidth=2.5)
        ax2[0].axhline(1.0, color='k', linestyle='dashed')

    if results_path is not None:
        output_file = open(results_path + ".summary.txt", 'w')
        print(output_header, file=output_file)

        print('\n' + "Eastward waveguide changes...", file=output_file)
        for j in range(364):
            if (eff_sndspd_pk[0][j] - 1.0) * (eff_sndspd_pk[0][j + 1] - 1) <= 0.0:
                if eff_sndspd_pk[0][j] > eff_sndspd_pk[0][j + 1]:
                    print('\t' + "Waveguide dissipates:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")", file=output_file)
                else:
                    print('\t' + "Waveguide forms:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")", file=output_file)

        print('\n' + "Westward waveguide changes...", file=output_file)
        for j in range(364):
            if (eff_sndspd_pk[1][j] - 1.0) * (eff_sndspd_pk[1][j + 1] - 1) <= 0.0:
                if eff_sndspd_pk[1][j] > eff_sndspd_pk[1][j + 1]:
                    print('\t' + "Waveguide dissipates:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")", file=output_file)
                else:
                    print('\t' + "Waveguide forms:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")", file=output_file)
        if include_ns:
            print('\n' + "Northward waveguide changes...", file=output_file)
            for j in range(364):
                if (eff_sndspd_pk[2][j] - 1.0) * (eff_sndspd_pk[2][j + 1] - 1) <= 0.0:
                    if eff_sndspd_pk[2][j] > eff_sndspd_pk[2][j + 1]:
                        print('\t' + "Waveguide dissipates:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")", file=output_file)
                    else:
                        print('\t' + "Waveguide forms:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")", file=output_file)

            print('\n' + "Southward waveguide changes...", file=output_file)
            for j in range(364):
                if (eff_sndspd_pk[3][j] - 1.0) * (eff_sndspd_pk[3][j + 1] - 1) <= 0.0:
                    if eff_sndspd_pk[3][j] > eff_sndspd_pk[3][j + 1]:
                        print('\t' + "Waveguide dissipates:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")", file=output_file)
                    else:
                        print('\t' + "Waveguide forms:", datetime.strptime('20' + "{:03d}".format(j), '%y%j').date().strftime('%B %d'), " (yday: " + str(j) + ", week: " + str(int(np.round(j/7))) + ")", file=output_file)

        output_file.close()

        f1.savefig(results_path + ".ess-ratio.png", dpi=300)
        if include_ns:
            f2.savefig(results_path + ".ess-ratio.NS.png", dpi=300)

    plt.show()

# ############################## #
#  Running infraga and NCPAprop  #
#     with multiple profiles     #
# ############################## #
def run_infraga(profs_path, results_file, pattern="*.met", cpu_cnt=None, geom="3d", bounces=25, inclinations=[1.0, 60.0, 1.0], azimuths=[-180.0, 180.0, 3.0], freq=0.1, z_grnd=0.0, rng_max=1000.0, src_loc=[0.0, 0.0, 0.0], infraga_path="", local_temp_dir=None, prof_format="zTuvdp", reverse_winds=False, topo_file=None, verbose=False):
    """
        Run the infraga -prop algorithm to compute path geometry
        statistics for BISL using a suite of specifications
        and combining results into single file

        Parameters
        ----------
        profs_path: string
            Path to atmospheric specification files
        results_file: string
            Path and name of file where results will be written
        pattern: string
            Pattern identifying atmospheric specification within profs_path location
        cpu_cnt: int
            Number of threads to use in OpenMPI implementation.  None runs non-OpenMPI version of infraga
        geom: string
            Defines geometry of the infraga simulations (3d" or "sph")
        bounces: int
            Maximum number of ground reflections to consider in ray tracing
        inclinations: iterable object
            Iterable of starting, ending, and step for ray launch inclination
        azimuths: iterable object
            Iterable of starting, ending, and step for ray launch azimuths
        freq: float
            Frequency to use for Sutherland Bass absorption calculation
        z_grnd: float
            Elevation of the ground surface relative to sea level
        rng_max: float
            Maximum propagation range for propagation paths
        src_loc: iterable object
            The horizontal (latitude and longitude) and altitude of the source
        infraga_path: string
            Location of infraGA executables
        local_temp_dir: string
            Path of local directory for storing individual infraGA results
        """

    if os.path.isfile(results_file):
        print(results_file + " already exists  --->  Skipping infraGA/GeoAc ray tracing runs...")
    else:
        if not ((geom == "3d") or (geom == "sph")):
            msg = "Incompatible geometry option for infraga: {}.  Options are 3d' and 'sph'".format(geom)
            warnings.warn(msg)
        else:
            if profs_path[-1] == "/":
                profs_path = profs_path[:-1]

            with tempfile.TemporaryDirectory(prefix='infraga_') as tmpdirname:
                if local_temp_dir is not None:
                    print("Writing individual infraGA results into local directory: " + local_temp_dir + '\n')
                    if not os.path.isdir(local_temp_dir):
                        os.mkdir(local_temp_dir)
                    tmpdirname = local_temp_dir
                else:
                    print('Created temp directory:', tmpdirname + '\n')
                    
                dir_files = np.sort(os.listdir(profs_path))
                for file_name in dir_files:
                    if fnmatch.fnmatch(file_name, pattern):
                        file_id = os.path.splitext(file_name)[0]

                        if os.path.isfile(tmpdirname + "/" + file_id + ".arrivals.dat"):
                            print("InfraGA/GeoAc arrivals for " + file_name + " already computed...")
                        else:
                            print("Generating ray paths for " + file_name)
                            if cpu_cnt:
                                command = "mpirun -np " + str(cpu_cnt) + " " + infraga_path + " infraga-accel-" + geom + " -prop "
                            else:
                                command = infraga_path + " infraga-" + geom + " -prop "

                            command = command + profs_path + "/" + file_name
                            command = command + " incl_min=" + str(inclinations[0]) + " incl_max=" + str(inclinations[1]) + " incl_step=" + str(inclinations[2])
                            command = command + " az_min=" + str(azimuths[0]) + " az_max=" + str(azimuths[1]) + " az_step=" + str(azimuths[2])
                            if geom == "sph":
                                command = command + " src_lat=" + str(src_loc[0]) + " src_lon=" + str(src_loc[1])
                            command = command + " src_alt=" + str(src_loc[2]) + " freq=" + str(freq) + " max_rng=" + str(rng_max) + " prof_format=" + str(prof_format)
                            command = command + " calc_amp=False" + " bounces=" + str(bounces) + " write_rays=false"
                            
                            if reverse_winds:
                                command = command + " reverse_winds=True"

                            if topo_file is not None:
                                command = command + " topo_file=" + topo_file 
                            else:
                                command = command + " z_grnd=" + str(z_grnd)

                            command = command + " output_id=" + tmpdirname + "/" + file_id

                            if not verbose:
                                command = command + " > /dev/null"

                            print(command)
                            subprocess.call(command, shell=True)

                command = "cat " + tmpdirname + "/*.arrivals.dat > " + results_file
                subprocess.call(command, shell=True)
                print("")



def run_ncpaprop(ncpaprop_method, profs_path, results_path, pattern="*.met", azimuths=[0.0, 359.0, 3.0], freq=0.1, z_grnd=0.0, rng_max=1000.0, rng_resol=1.0, src_loc=[40.0, -110.0, 0.0], ncpaprop_path="", use_topo=False, local_temp_dir=None, keep_lossless=False, reverse_winds=False, cpu_cnt=1):
    """
        Run one of the NCPAprop methods to compute transmission
        loss values for a suite of atmospheric specifications at
        a set of frequency values
        
        Note: the methods here use the ncpaprop-release version that includes 
        an option for --filetag that writes output into a specific location 
        and enables simultaneous calculations via subprocess.popen()

        Example ncpaprop function calls (usable from infraGA/GeoAc examples directory):

            Effective Sound Speed Modal (modess)
                modess --multiprop --atmosfile ToyAtmo.met --zground_km 1.0 --freq 0.5 --azimuth_start 0 --azimuth_end 360 --azimuth_step 3 --maxrange_km 1000 --Nrng_steps 1000 --filetag test
            
            ePade Parabolic Equation (epape)
                epape --multiprop --starter self --atmosfile ToyAtmo.met --groundheight_km 1.0 --freq 0.5 --azimuth_start 0 --azimuth_end 360 --azimuth_step 3 --maxrange_km 1000 --Nrng_steps 1000  --filetag test
                epape --singleprop --topo --starter self --atmosfile ToyAtmo.met --topofile topo/line.dat --freq 0.5 --azimuth 0 --maxrange_km 1000 --Nrng_steps 1000 --sourceheight_km 0.0 --filetag test2

        Note on defining flat ground elevation:
            For some reason, most parameters are the same between modess and pape (e.g., freq, azimuth_start, Nrng_steps), but the parameter to define the ground elevation is difference (--zground_km in modess vs. --groundheight_km in epape).  I might email Claus about this...

        Notes on including terrain:
            Including terrain requires individual azimuth runs to match up azimuths with correct terrain lines (can't use Nx2D)

            There is a utility function in InfraGA/GeoAc that extracts terrain lines from etopo1 that can be used here.
                infraga utils extract-terrain --geom line --lat1 40.0 --lon1 -102.5 --azimuth -90.0 --range 750.0 --output-file line_topo.dat
                        
            Terrain line file requires header info (included in InfraGA util output):
                # 'infraga extract-terrain --geom line' summary:
                # lat: 37.114912493
                # lon: -116.069089728
                # azimuth: 59.0
                # range: 1000.0
                #% r, km
                #% z, km
                0.0  0.8320394823
                1.8518518518517908	0.8875136241398117
                3.703703703703796	0.9092146274618726
                5.555555555556105	0.9170049733706201
                ...

        Parameters
        ----------
        ncpaprop_method: string
            Method in NCPAprop to use (currently supports 'modess' or 'epape')
        profs_path: string
            Path to atmospheric specification files
        results_file: string
            Path and name of file where results will be written
        pattern: string
            Pattern identifying atmospheric specification within profs_path location
        azimuths: iterable object
            Iterable of starting, ending, and step for propagation azimuths
        freq: float
            Frequency for simulation
        z_grnd: float
            Elevation of the ground surface relative to sea level
        rng_max: float
            Maximum propagation range for propagation paths
        rng_resol: float
            Resolution for range ouputs (defaults to 1 km)
        src_loc: iterable object
            The horizontal (latitude and longitude) and altitude of the source
        ncpaprop_path: string
            Path to NCPAprop binaries (if not on path)
        use_topo: bool
            Flag to turn on terrain inclusion for simulations
        local_temp_dir: string
            Path of local directory for storing individual NCPAprop results
        keep_lossless: boolean
            Flag to keep the lossless (no absorption) results
        cpu_cnt : integer
            Number of CPUs to use in subprocess.popen loop for simultaneous calculations
        """

    if profs_path[-1] != "/":
        profs_path = profs_path + "/"

    if ncpaprop_method == "modess":
        output_suffix = ".nm"
    else:
        output_suffix = ".pe"

    if os.path.isfile(results_path + output_suffix):
        print(results_path + output_suffix + " already exists  --->  Skipping NCPAprop simulations...")
    else:
        with tempfile.TemporaryDirectory(prefix='infraga_') as tmpdirname:
            if local_temp_dir is not None:
                print("Writing individual NCPAprop results into local directory: " + local_temp_dir + '\n')
                if not os.path.isdir(local_temp_dir):
                    os.mkdir(local_temp_dir)
                tmpdirname = local_temp_dir
            else:
                print('Created temp directory:', tmpdirname + '\n')
                
            if use_topo:
                reverse_winds = True
                output_suffix = ".pe"
                command_prefix = ncpaprop_path + " ePape --singleprop --topo --starter self"

                # Add methods here to cycle through azimuths using appropriate terrain lines

                topo_extract = "infraga utils extract-terrain --geom nx2d --show-terrain False"
                topo_extract = topo_extract + " --lat1 " + str(src_loc[0]) + " --lon1 " + str(src_loc[1])
                topo_extract = topo_extract + " --nx2d-resol " + str(azimuths[-1]) + " --range " + str(rng_max + 1.0)
                topo_extract = topo_extract + " --output-file " + tmpdirname + "/terrain"
                subprocess.call(topo_extract, shell=True)

                dir_files = np.sort(os.listdir(profs_path))
                for file_name in dir_files:
                    if fnmatch.fnmatch(file_name, pattern):
                        file_id = os.path.splitext(file_name)[0]

                        print('\n' + "Running simulations using --atmosfile " + profs_path + file_name)
                        command_list = []
                        for az_val in np.arange(*azimuths):
                            command = command_prefix + " --atmosfile " + profs_path + file_name
                            command = command + " --topofile " + tmpdirname + "/terrain_%06.2fdeg.line.dat" % az_val
                            command = command + " --freq " + str(freq) + " --azimuth " + str(az_val)
                            command = command + " --maxrange_km " + str(rng_max) + " --Nrng_steps " + str(int(rng_max / rng_resol))
                            command = command + " --filetag " + tmpdirname + "/" + file_id + "_%.3fHz" % freq + "_%06.2fdeg" % az_val
                            command = command + " > /dev/null"
                    
                            command_list = command_list + [command]

                        for j in range(0, len(command_list), cpu_cnt):              
                            if cpu_cnt==1 or j + 1 == len(command_list):
                                print('\t' + "Running NCPAprop " + ncpaprop_method + " for azimuth " + str(j * azimuths[2]))
                            elif cpu_cnt==2 or j + 2 == len(command_list):
                                print('\t' + "Running NCPAprop " + ncpaprop_method + " for azimuths " + str(j * azimuths[2]) + ", " + str((j + 1) * azimuths[2]))
                            else:
                                print('\t' + "Running NCPAprop " + ncpaprop_method + " for azimuths " + str(j * azimuths[2]) + " - " + str(min(j + cpu_cnt, len(command_list)) * azimuths[2]))

                            procs_list = [subprocess.Popen(cmd, shell=True) for cmd in command_list[j:j + cpu_cnt]]
                            for proc in procs_list:
                                proc.communicate()
                                proc.wait()
                
                        print('\t' + "Merging azimuth results...")
                        command = "cat " + tmpdirname + "/" + file_id + "_%.3f" % freq + "Hz*" + output_suffix + " > " + tmpdirname + "/" + file_id + "_%.3fHz.tloss_2d" % freq + output_suffix
                        subprocess.call(command, shell=True)

                print('\t' + "Combining transmission loss predictions...")
                command = "cat " + tmpdirname + "/*_%.3f" % freq + "Hz.tloss_2d" + output_suffix + " > " + results_path + output_suffix
                print('\t' + command)
                subprocess.call(command, shell=True)                        

            else:
                if ncpaprop_method == "modess":
                    command_prefix = ncpaprop_path + " modess --multiprop --zground_km " + str(z_grnd)
                else:
                    command_prefix = ncpaprop_path + " epape --multiprop --starter self --groundheight_km " + str(z_grnd)

                dir_files = np.sort(os.listdir(profs_path))
                print("Running NCPAprop " + ncpaprop_method + " for atmospheric specifications in " + profs_path)
                command_list = []
                for file_name in dir_files:
                    if fnmatch.fnmatch(file_name, pattern):
                        file_id = os.path.splitext(file_name)[0]
                        if os.path.isfile(tmpdirname + "/" + file_id + "_%.3fHz" % freq + ".Nby2D_tloss_1d" + output_suffix):
                            print("NCPAprop results for " + file_name + " at " + str(freq) + " Hz already computed...")
                        else:
                            command = command_prefix + " --atmosfile " + profs_path + file_name
                            command = command + " --freq " + str(freq)
                            command = command + " --maxrange_km " + str(rng_max) + " --Nrng_steps " + str(int(rng_max / rng_resol))
                            command = command + " --azimuth_start " + str(azimuths[0]) + " --azimuth_end " + str(azimuths[1]) + " --azimuth_step " + str(azimuths[2])
                            command = command + " --filetag " + tmpdirname + "/" + file_id + "_%.3fHz" % freq
                            command = command + " > /dev/null"

                            command_list = command_list + [command]

                for j in range(0, len(command_list), cpu_cnt):              
                    if cpu_cnt==1 or j + 1 == len(command_list):
                        print('\t' + "Running NCPAprop " + ncpaprop_method + " for sample " + str(j + 1) + "of " + str(len(command_list)))
                    elif cpu_cnt==2 or j + 2 == len(command_list):
                        print('\t' + "Running NCPAprop " + ncpaprop_method + " for samples " + str(j + 1) + ", " + str(j + 2) + " of " + str(len(command_list)))
                    else:
                        print('\t' + "Running NCPAprop " + ncpaprop_method + " for samples " + str(j + 1) + " - " + str(min(j + cpu_cnt, len(command_list))) + " of " + str(len(command_list)))

                    procs_list = [subprocess.Popen(cmd, shell=True) for cmd in command_list[j:j + cpu_cnt]]
                    for proc in procs_list:
                        proc.communicate()
                        proc.wait()

                print('\t' + "Combining transmission loss predictions...")
                command = "cat " + tmpdirname + "/*_%.3f" % freq + "Hz*" + output_suffix + " > " + results_path + output_suffix

                print('\t' + command)
                subprocess.call(command, shell=True)

                if keep_lossless:
                    command = "cat " + tmpdirname + "/*_%.3f" % freq + "Hz*.lossless" + output_suffix + " > " + results_path + ".lossless" + output_suffix
                    print('\t' + command)
                    subprocess.call(command, shell=True)




# ############################ #
#          Stochastic          #
#      Propagation Models      #
# ############################ #
def find_azimuth_bin(az, bin_cnt=16):
    """
        Identify the azimuth bin index given some specified number of bins

        Parameters
        ----------
        az: float
            Azimuth in degrees
        bin_cnt: int
            Number of bins used in analysis

        Returns
        -------
        index: int
            Index of azimuth bin
    """

    # wrap angles to range (-180, 180) degrees
    reduced = np.degrees(np.arctan2(np.sin(np.radians(az)), np.cos(np.radians(az))))

    # Identify the index of the nearest bin
    bins = np.arange(-180.0, 180.0, 360.0 / (bin_cnt * 2.0))
    result = np.asarray(np.digitize(reduced, bins) / 2)
    result[result >= bin_cnt] = 0

    return result.astype(int)


class PathGeometryModel(object):
    """
        Propagation path geometry statistics computed using ray tracing
        analysis on a suite of specifications includes celerity-range and
        azimuth deviation/scatter statistics
    """

    _az_bin_cnt = 16
    _rng_max = 1000.0

    _default_az_dev_std = 4.0
    _min_az_dev_std = 1.0

    _tropo_strat_bnd = 1.0 / 0.32
    _strat_therm_bnd = 1.0 / 0.26
    _bnd_overlap = 0.05

    _wts0 = np.array([0.0539, 0.0899, 0.8562])
    _std0 = np.array([0.066, 0.08, 0.33])
    _mns0 = np.array([1.0 / 0.33, 1.0 / 0.29, 1.0 / 0.26])

    _win_min_pts = 25
    _rcel_std_min = 0.05
    _rcel_wt_min = 1.0e-2

    def __init__(self):
        self.rngs = np.array([])

        self._rcel_wts = []
        self._rcel_mns = []
        self._rcel_std = []

        self.az_dev_mns = []
        self.az_dev_std = []

    def eval_rcel_gmm(self, rng, rcel, az):
        """
            Evaluate reciprocal celerity Gaussian Mixture Model (GMM)
            at specified range, reciprocal celerity, and azimuth

            Parameters
            ----------
            rng: float
                Range from source
            rcel: float
                Reciprocal celerity (travel time divided by propagation range)
            az: float
                Propagation azimuth (relative to North)

            Returns
            -------
            pdf: float
                Probability of observing an infrasonic arrival with specified celerity at specified range and azimuth

        """
        rng_eval = np.array(rng)
        rng_eval[rng_eval > self._rng_max] = self._rng_max

        if len(np.atleast_1d(rng)) == 1:
            n_az = find_azimuth_bin(az, self._az_bin_cnt)
            fit_rcel_wts = np.array([f(rng_eval) for f in self._rcel_wts[n_az]])
            fit_rcel_mns = np.array([f(rng_eval) for f in self._rcel_mns[n_az]])
            fit_rcel_std = np.array([f(rng_eval) for f in self._rcel_std[n_az]])
            result = np.sum(fit_rcel_wts / fit_rcel_std * norm.pdf((rcel - fit_rcel_mns) / fit_rcel_std))
        else:
            mn = np.empty((len(rng_eval), 3))
            vr = np.empty((len(rng_eval), 3))
            wt = np.empty((len(rng_eval), 3))

            az_indices = find_azimuth_bin(az, self._az_bin_cnt)
            for n_az in range(self._az_bin_cnt):
                mask = tuple([az_indices == n_az])
                if np.any(mask):
                    mn[mask] = np.array([f(rng_eval[mask]) for f in self._rcel_mns[n_az]]).T
                    vr[mask] = np.array([f(rng_eval[mask]) for f in self._rcel_std[n_az]]).T
                    wt[mask] = np.array([f(rng_eval[mask]) for f in self._rcel_wts[n_az]]).T

            result = np.sum(wt / vr * norm.pdf((np.array([rcel] * 3).T - mn) / vr), axis=1)
        return result

    def eval_az_dev_mn(self, rng, az):
        """
            Evaluate the mean back azimuth deviation at a given range
            and propagation azimuth

            Parameters
            ----------
            rng: float
                Range from source
            az: float
                Propagation azimuth (relative to North)

            Returns
            -------
            bias: float
                Predicted bias in the arrival back azimuth at specified arrival range and azimuth

        """
        if len(np.atleast_1d(rng)) == 1:
            return self.az_dev_mns[find_azimuth_bin(az, self._az_bin_cnt)](min(rng, self._rng_max))
        else:
            rng_eval = np.array(rng)
            rng_eval[rng_eval > self._rng_max] = self._rng_max
            az_indices = find_azimuth_bin(az, self._az_bin_cnt)

            mn = np.empty_like(rng_eval)
            for n_az in range(self._az_bin_cnt):
                mask = az_indices == n_az
                if np.any(mask):
                    mn[mask] = self.az_dev_mns[n_az](rng_eval[mask])
            return mn

    def eval_az_dev_std(self, rng, az):
        """
            Evaluate the standard deviation of the back azimuth at a given range
            and propagation azimuth


            Parameters
            ----------
            rng: float
                Range from source
            az: float
                Propagation azimuth (relative to North)

            Returns
            -------
            stdev: float
                Standard deviation of arrival back azimuths at specified range and azimuth

        """
        if len(np.atleast_1d(rng)) == 1:
            return self.az_dev_std[find_azimuth_bin(az, self._az_bin_cnt)](min(rng, self._rng_max))
        else:
            rng_eval = np.array(rng)
            rng_eval[rng_eval > self._rng_max] = self._rng_max
            az_indices = find_azimuth_bin(az, self._az_bin_cnt)

            vr = np.empty_like(rng_eval)
            for n_az in range(self._az_bin_cnt):
                mask = az_indices == n_az
                if np.any(mask):
                    vr[mask] = self.az_dev_std[n_az](rng_eval[mask])
            return vr

    def build(self, arrivals_file, output_file, show_fits=False, rng_width=50.0, rng_spacing=10.0, rng_max=1000.0, geom="3d", src_loc=[0.0, 0.0, 0.0], min_turning_ht=0.0, az_bin_cnt=16, az_bin_wdth=30.0, station_centered=False):
        """
            Construct propagation statistics from a ray tracing arrival file (concatenated from
            multiple runs most likely) and output a path geometry model


            Parameters
            ----------
            arrivals_file: string
                Path to file containing infraGA/GeoAc arrival information
            output_file: string
                Path to file where results will be saved
            show_fits: boolean
                Option ot visualize model construction (for QC purposes)
            rng_width: float
                Range bin width in kilometers
            rng_spacing: float
                Spacing between range bins in kilometers
            rng_max: float
                Maximum range in model
            geom: string
                Geometry used in infraGA/GeoAc simulation.  Options are "3d" and "sph"
            src_loc: iterable
                [x, y, z] or [lat, lon, elev] location of the source used in infraGA/GeoAc simulations.  Note: '3d' simulations assume source at origin.
            min_turning_ht: float
                Minimum turning height used to filter out boundary layer paths if not of interest
            az_bin_cnt: int
                Number of azimuth bins to use in analysis
            az_bin_width: float
                Azimuth bin width in degrees for analysis
            station_centered: bool
                Flag to build station-centered model via back projection ray tracing (needed for inclusion of terrain)
                
        """
        if os.path.isfile(output_file):
            print(output_file + " already exists  --->  Skipping path geometry model construction...")
        else:
            if not ((geom == "3d") or (geom == "sph")):
                msg = "Incompatible geometry option for infraga: {}.  Options are 3d' and 'sph'".format(geom)
                warnings.warn(msg)
            else:
                print('Building celerity and azimuth statistics from file:', arrivals_file)

                self._az_bin_cnt = az_bin_cnt

                # define range bins and parameter arrays
                self._rng_max = rng_max
                rng_bins = np.arange(0.0, self._rng_max, rng_spacing)
                rng_cnt = len(rng_bins)

                az_dev_mns = np.empty((self._az_bin_cnt, rng_cnt))
                az_dev_std = np.empty((self._az_bin_cnt, rng_cnt))

                rcel_wts = np.empty((self._az_bin_cnt, rng_cnt, 3))
                rcel_mns = np.empty((self._az_bin_cnt, rng_cnt, 3))
                rcel_std = np.empty((self._az_bin_cnt, rng_cnt, 3))

                # load infraGA/GeoAc predictions
                if geom == "3d":
                    theta, phi, n, x, y, t, cel, z_max, incl, back_az, amp_geo, amp_atmo = np.loadtxt(arrivals_file, unpack=True)

                    rngs = np.sqrt(x**2 + y**2)
                    if station_centered:
                        az = 90.0 - np.degrees(np.arctan2(-y, -x))
                        az_dev = (90.0 - np.degrees(np.arctan2(y, x))) - phi 
                    else:
                        az = 90.0 - np.degrees(np.arctan2(y, x))
                        az_dev = (90.0 - np.degrees(np.arctan2(-y, -x))) - back_az
 
                else:
                    theta, phi, n, lats, lons, t, cel, z_max, incl, back_az, amp_geo, amp_atmo = np.loadtxt(arrivals_file, unpack=True)

                    if station_centered:
                        az, temp_az, rngs = sph_proj.inv(lons, lats, np.array([src_loc[1]] * len(lons)), np.array([src_loc[0]] * len(lats)))
                        rngs = rngs / 1000.0
                        az_dev = np.array(temp_az) - phi
                    else:
                        az, temp_az, rngs = sph_proj.inv(np.array([src_loc[1]] * len(lons)), np.array([src_loc[0]] * len(lats)), lons, lats)
                        rngs = rngs / 1000.0
                        az_dev = np.array(temp_az) - back_az

                rcel = 1.0 / cel
 
                # wrap angles to +/- 180 degrees
                phi[phi > 180.0] -= 360.0
                phi[phi < -180.0] += 360.0

                az[az > 180.0] -= 360.0
                az[az < -180.0] += 360.0

                az_dev[az_dev > 180.0] -= 360.0
                az_dev[az_dev < -180.0] += 360.0

                print("phi lims:", min(phi), max(phi))
                print("az lims:", min(az), max(az))
                print("az_dev lims:", min(az_dev), max(az_dev))

                # Cycle through azimuth bins creating fit
                az_wts = np.empty(self._az_bin_cnt)
                for n_az in range(self._az_bin_cnt):
                    if n_az == 0:
                        center = -180.0
                        az_mask = np.logical_or(az > 180.0 - az_bin_wdth / 2.0, az < -180.0 + az_bin_wdth / 2.0)
                    else:
                        center = -180 + (360.0 / self._az_bin_cnt) * n_az
                        az_mask = np.logical_and(center - az_bin_wdth / 2.0 <= az, az <= center + az_bin_wdth / 2.0)

                    # combine azimuth and turning height masks
                    az_mask = np.logical_and(az_mask, min_turning_ht < (z_max - src_loc[2]))
                    az_wts[n_az] = float(len(rngs[az_mask])) / float(len(rngs))

                    # display work if show_fits is enabled
                    if show_fits:
                        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 11))

                        ax1.set_xlim([0.0, self._rng_max])
                        ax1.set_ylim([0.2, 0.4])

                        ax3.set_xlim([0.0, self._rng_max])
                        ax3.set_ylim([-12.0, 12.0])

                        ax1.set_xlabel('Range [km]')
                        ax1.set_ylabel('Celerity [km/s]')

                        ax2.set_xlabel('Celerity [km/s]')
                        ax2.set_ylabel('Probability')

                        ax3.set_xlabel('Range [km]')
                        ax3.set_ylabel('Azimuth Deviation [degrees]')

                        ax4.set_xlabel('Azimuth Deviation [degrees]')
                        ax4.set_ylabel('Probability')

                        plt.suptitle('Path Geometry Statistics (' + str(center) + ')', fontsize=18)
                        plt.show(block=False)

                        ax1.plot(rngs[az_mask], 1.0 / rcel[az_mask], 'k.', markersize=2.0)
                        ax1.set_title('Celerity-range scatter')
                        plt.pause(0.01)

                        ax3.plot(rngs[az_mask], az_dev[az_mask], 'k.', markersize=2.0)
                        ax3.set_title('Back-azimuth deviation')
                        plt.pause(0.01)

                    # Compute fits (1D interpolations of parameters)
                    rng_wts = np.empty(rng_cnt)
                    for nr in range(rng_cnt):
                        if show_fits:
                            window1 = ax1.axvspan(rng_bins[nr], rng_bins[nr] + rng_width, facecolor='b', alpha=0.5)
                            window2 = ax3.axvspan(rng_bins[nr], rng_bins[nr] + rng_width, facecolor='b', alpha=0.5)
                            plt.pause(0.1)

                        # define combined azimuth and range mask
                        rng_mask = np.logical_and(rng_bins[nr] <= rngs, rngs <= rng_bins[nr] + rng_width)
                        rng_az_mask = np.logical_and(rng_mask, az_mask)
                        rng_wts[nr] = float(len(rngs[rng_az_mask])) / len(rngs[az_mask])

                        # set azimuth deviation fit
                        if az_dev[rng_az_mask].shape[0] > self._win_min_pts:
                            az_dev_mns[n_az][nr] = np.mean(az_dev[rng_az_mask], dtype=np.float64)
                            az_dev_std[n_az][nr] = max(np.std(az_dev[rng_az_mask], dtype=np.float64), self._min_az_dev_std)
                        else:
                            if nr == 0:
                                az_dev_mns[n_az][nr] = 0.0
                                az_dev_std[n_az][nr] = self._default_az_dev_std
                            else:
                                az_dev_mns[n_az][nr] = az_dev_mns[n_az][nr - 1] * 0.75
                                az_dev_std[n_az][nr] = self._default_az_dev_std + (az_dev_std[n_az][nr - 1] - self._default_az_dev_std) * 0.5

                        # set tropospheric contribution to reciprocal celerity fit
                        combo_mask = np.logical_and(rng_az_mask, rcel < self._tropo_strat_bnd + self._bnd_overlap)

                        wt1 = float(len(rngs[combo_mask]))
                        if rcel[combo_mask].shape[0] > self._win_min_pts:
                            rcel_mns[n_az][nr][0] = np.mean(rcel[combo_mask], dtype=np.float64)
                            rcel_std[n_az][nr][0] = max(np.std(rcel[combo_mask], dtype=np.float64) * 1.25, self._rcel_std_min)
                        else:
                            if nr == 0:
                                rcel_mns[n_az][nr][0] = self._mns0[0]
                                rcel_std[n_az][nr][0] = self._std0[0]
                            else:
                                rcel_mns[n_az][nr][0] = self._mns0[0] + (rcel_mns[n_az][nr - 1][0] - self._mns0[0]) * 0.25
                                rcel_std[n_az][nr][0] = self._std0[0] + (rcel_std[n_az][nr - 1][0] - self._std0[0]) * 0.25

                        # set stratospheric contribution to reciprocal celerity fit
                        strat_mask = np.logical_and(self._tropo_strat_bnd - self._bnd_overlap <= rcel, rcel <= self._strat_therm_bnd + self._bnd_overlap)
                        combo_mask = np.logical_and(rng_az_mask, strat_mask)

                        wt2 = float(len(rngs[combo_mask]))
                        if rcel[combo_mask].shape[0] > self._win_min_pts:
                            rcel_mns[n_az][nr][1] = np.mean(rcel[combo_mask], dtype=np.float64)
                            rcel_std[n_az][nr][1] = max(np.std(rcel[combo_mask], dtype=np.float64) * 1.25, self._rcel_std_min)
                        else:
                            if nr == 0:
                                rcel_mns[n_az][nr][1] = self._mns0[1]
                                rcel_std[n_az][nr][1] = self._std0[1]
                            else:
                                rcel_mns[n_az][nr][1] = self._mns0[1] + (rcel_mns[n_az][nr - 1][1] - self._mns0[1]) * 0.25
                                rcel_std[n_az][nr][1] = self._std0[1] + (rcel_std[n_az][nr - 1][1] - self._std0[1]) * 0.25

                        # set thermospheric contribution to reciprocal celerity fit
                        combo_mask = np.logical_and(rng_az_mask, self._strat_therm_bnd - self._bnd_overlap < rcel)

                        wt3 = float(len(rngs[combo_mask]))
                        if rcel[combo_mask].shape[0] > self._win_min_pts:
                            rcel_mns[n_az][nr][2] = np.mean(rcel[combo_mask], dtype=np.float64)
                            rcel_std[n_az][nr][2] = max(np.std(rcel[combo_mask], dtype=np.float64) * 1.25, self._rcel_std_min)
                        else:
                            if nr == 0:
                                rcel_mns[n_az][nr][2] = self._mns0[2]
                                rcel_std[n_az][nr][2] = self._std0[2]
                            else:
                                rcel_mns[n_az][nr][2] = self._mns0[2] + (rcel_mns[n_az][nr - 1][2] - self._mns0[2]) * 0.25
                                rcel_std[n_az][nr][2] = self._std0[2] + (rcel_std[n_az][nr - 1][2] - self._std0[2]) * 0.25

                        # set weights of reciprocal celerity distribution
                        if len(rngs[rng_az_mask]) > self._win_min_pts:
                            rcel_wts[n_az][nr][0] = wt1 / (wt1 + wt2 + wt3)
                            rcel_wts[n_az][nr][1] = wt2 / (wt1 + wt2 + wt3)
                            rcel_wts[n_az][nr][2] = wt3 / (wt1 + wt2 + wt3)
                        else:
                            rcel_wts[n_az][nr] = [0.0, 0.0, 0.0]

                        if show_fits:
                            ax4.cla()
                            ax4.set_xlim([-12.5, 12.5])
                            ax4.plot(np.linspace(-12.5, 12.5, 1000), 1.0 / az_dev_std[n_az][nr] * norm.pdf((np.linspace(-12.5, 12.5, 1000) - az_dev_mns[n_az][nr]) / az_dev_std[n_az][nr]), linewidth=4.0, color='Blue')
                            plt.pause(0.01)

                            cel_vals = np.linspace(0.2, 0.4, 1000)
                            cel_dist1 = (rcel_wts[n_az][nr][0] / len(rngs)) / rcel_std[n_az][nr][0] * norm.pdf((1.0 / cel_vals - rcel_mns[n_az][nr][0]) / rcel_std[n_az][nr][0])
                            cel_dist2 = (rcel_wts[n_az][nr][1] / len(rngs)) / rcel_std[n_az][nr][1] * norm.pdf((1.0 / cel_vals - rcel_mns[n_az][nr][1]) / rcel_std[n_az][nr][1])
                            cel_dist3 = (rcel_wts[n_az][nr][2] / len(rngs)) / rcel_std[n_az][nr][2] * norm.pdf((1.0 / cel_vals - rcel_mns[n_az][nr][2]) / rcel_std[n_az][nr][2])

                            ax2.cla()
                            ax2.set_xlim([0.2, 0.4])
                            ax2.plot(cel_vals, cel_dist1, linewidth=2.0, color='Green')
                            ax2.plot(cel_vals, cel_dist2, linewidth=2.0, color='Green')
                            ax2.plot(cel_vals, cel_dist3, linewidth=2.0, color='Green')
                            ax2.plot(cel_vals, cel_dist1 + cel_dist2 + cel_dist3, linewidth=4.0, color='Blue')
                            plt.pause(0.01)

                            window1.remove()
                            window2.remove()

                    rng_wts[rng_wts > 4.0 / len(rng_wts)] = 4.0 / len(rng_wts)
                    print("rng_wts:", rng_wts)

                    for nr in range(rng_cnt):
                        rcel_wts[n_az][nr] *= rng_wts[nr] / np.sum(rng_wts)

                    plt.close('all')

                
                # Normalize weights by total arrivals at all azimuths
                print("az_wts:", az_wts)

                for n_az in range(self._az_bin_cnt):
                    rcel_wts[n_az] *= az_wts[n_az] / np.sum(az_wts)

                priors = [0] * 6
                priors[0] = rng_bins
                priors[1] = az_dev_mns
                priors[2] = az_dev_std
                priors[3] = rcel_mns
                priors[4] = rcel_std
                priors[5] = rcel_wts

                pickle.dump(priors, open(output_file, "wb"))

    def load(self, model_file, smooth=False):
        """
        Load a path geometry model file for use

        Parameters
        ----------

        model_file: string
            Path to PGM file constructed using stochprop.propagation.PathGeometryModel.build()
        smooth: boolean
            Option to use scipy.signal.savgol_filter to smooth discrete GMM parameters along range

        """

        fit_params = pickle.load(open(model_file, "rb"), encoding='latin1')
        self._az_bin_cnt = len(fit_params[1])
        self._rng_max = max(fit_params[0])

        print("Loading path geometry model from " + model_file)
        print('\t' + "Azimuth bin count: " + str(self._az_bin_cnt))
        print('\t' + "Maximum range: " + str(self._rng_max))

        self.az_dev_mns = [0] * self._az_bin_cnt
        self.az_dev_std = [0] * self._az_bin_cnt

        self._rcel_wts = [0] * self._az_bin_cnt
        self._rcel_mns = [0] * self._az_bin_cnt
        self._rcel_std = [0] * self._az_bin_cnt

        for n_az in range(self._az_bin_cnt):
            self._rcel_mns[n_az] = [0] * 3
            self._rcel_std[n_az] = [0] * 3
            self._rcel_wts[n_az] = [0] * 3

        if smooth:
            for n_az in range(self._az_bin_cnt):
                self.az_dev_mns[n_az] = interp1d(fit_params[0], savgol_filter(fit_params[1][n_az], 5, 3), kind='cubic')
                self.az_dev_std[n_az] = interp1d(fit_params[0], savgol_filter(fit_params[2][n_az], 5, 3), kind='cubic')

                self._rcel_mns[n_az] = [0] * 3
                self._rcel_std[n_az] = [0] * 3
                self._rcel_wts[n_az] = [0] * 3

                for j in range(3):
                    self._rcel_mns[n_az][j] = interp1d(fit_params[0], savgol_filter(fit_params[3][n_az][:, j], 5, 3), kind='cubic')
                    self._rcel_std[n_az][j] = interp1d(fit_params[0], savgol_filter(fit_params[4][n_az][:, j], 5, 3), kind='cubic')
                    self._rcel_wts[n_az][j] = interp1d(fit_params[0], savgol_filter(fit_params[5][n_az][:, j], 5, 3), kind='cubic')
        else:
            for n_az in range(self._az_bin_cnt):
                self.az_dev_mns[n_az] = interp1d(fit_params[0], fit_params[1][n_az], kind='cubic')
                self.az_dev_std[n_az] = interp1d(fit_params[0], fit_params[2][n_az], kind='cubic')

                self._rcel_mns[n_az] = [0] * 3
                self._rcel_std[n_az] = [0] * 3
                self._rcel_wts[n_az] = [0] * 3

                for j in range(3):
                    self._rcel_mns[n_az][j] = interp1d(fit_params[0], fit_params[3][n_az][:, j], kind='cubic')
                    self._rcel_std[n_az][j] = interp1d(fit_params[0], fit_params[4][n_az][:, j], kind='cubic')
                    self._rcel_wts[n_az][j] = interp1d(fit_params[0], fit_params[5][n_az][:, j], kind='cubic')

    def display(self, file_id=None, subtitle=None, show_colorbar=True, hold_fig=False, cmap_max=None):
        """
        Display the propagation geometry statistics

        Parameters
        ----------
        file_id: string
            File prefix to save visualization
        subtitle: string
            Subtitle used in figures

        """

        resol = 100
        rngs = np.linspace(0.0, self._rng_max, resol)
        bias = np.empty([resol])
        width = np.empty([resol])

        bias_color = 'Blue'
        var_color = 'LightBlue'

        compass_file = find_spec('stochprop').submodule_search_locations[0] + '/resources/compass.png'

        f1, ax = plt.subplots(3, 3, figsize=(12, 9))

        for n1, n2 in itertools.product(list(range(3)), repeat=2):
            if n1 != 1 or n2 != 1:
                ax[n1, n2].set_ylim([-10, 10])
                ax[n1, n2].set_xlim([0, self._rng_max])
                # ax[n1, n2].set_xticks(np.arange(0, self._rng_max, 5))
            if n2 != 0:
                ax[n1, n2].set_yticklabels([])
            if n1 != 2:
                ax[n1, n2].set_xticklabels([])

        img = mpimg.imread(compass_file)
        ax[1, 1].axis('off')
        ax[1, 1].imshow(img)

        ax[2, 1].set_xlabel('Range [km]')
        ax[1, 0].set_ylabel('Azimuth Deviation [deg]')

        if subtitle:
            plt.suptitle(("Azimuth Deviation Statistics" + '\n' + subtitle), fontsize=22)
        else:
            plt.suptitle("Azimuth Deviation Statistics", fontsize=22)
        plt.show(block=False)

        def plot_az_stats(axis_id, azimuth):
            bias = self.eval_az_dev_mn(rngs, [azimuth] * len(rngs))
            width = self.eval_az_dev_std(rngs, [azimuth] * len(rngs))
            axis_id.fill_between(rngs, bias - 2.0 * width, bias + 2.0 * width, facecolor=var_color)
            axis_id.plot(rngs, bias, linewidth=2.0, color=bias_color)
            axis_id.plot(rngs, [0] * len(rngs), 'k:')
            plt.pause(0.1)

        plot_az_stats(ax[0, 0], -45.0)
        plot_az_stats(ax[0, 1], 0.0)
        plot_az_stats(ax[0, 2], 45.0)
        plot_az_stats(ax[1, 0], -90.0)
        plot_az_stats(ax[1, 2], 90.0)
        plot_az_stats(ax[2, 0], -135.0)
        plot_az_stats(ax[2, 1], 180.0)
        plot_az_stats(ax[2, 2], 135.0)

        if file_id:
            plt.savefig(file_id + "_az-dev.png", bbox_inches='tight', dpi=200)

        # Plot celerity-range statistics
        cels = np.linspace(0.2, 0.4, resol)
        rngs = np.linspace(0, self._rng_max, resol)
        R, V = np.meshgrid(rngs, cels)
        R = R.flatten()
        V = V.flatten()

        pdf = np.empty([resol, resol])

        palette = cm.nipy_spectral_r

        if cmap_max is not None:
            pdf_max = cmap_max
        else:
            pdf_max = 0.0
            for az in np.arange(-180.0, 180.0, 45.0):
                pdf_max = max(pdf_max, max(self.eval_rcel_gmm(R, 1.0 / V, [az] * len(R))))

        f2, ax = plt.subplots(3, 3, figsize=(12, 9))

        for n1, n2 in itertools.product(list(range(3)), repeat=2):
            if n1 != 1 or n2 != 1:
                ax[n1, n2].set_ylim([0.2, 0.4])
                ax[n1, n2].set_xlim([0, self._rng_max])
                # ax[n1, n2].set_xticks([0, 250, 500, 750, 1000])
            if n2 != 0:
                ax[n1, n2].set_yticklabels([])
            if n1 != 2:
                ax[n1, n2].set_xticklabels([])

        img = mpimg.imread(compass_file)
        ax[1, 1].axis('off')
        ax[1, 1].imshow(img)

        ax[2, 1].set_xlabel('Range [km]')
        ax[1, 0].set_ylabel('Celerity [km/s]')

        if subtitle:
            plt.suptitle(("Celerity-Range Statistics" + '\n' + subtitle), fontsize=22)
        else:
            plt.suptitle("Celerity-Range Statistics", fontsize=22)
        plt.show(block=False)

        pdf = self.eval_rcel_gmm(R, 1.0 / V, [-45.0] * len(R))
        rngcel_plot = ax[0, 0].scatter(R, V, c=pdf, cmap=palette, marker=".", alpha=1.0, edgecolor='none', vmin=0.0, vmax=pdf_max)

        if show_colorbar:            
            cbar = f1.colorbar(rngcel_plot, ax=[ax[0, 2], ax[1, 2], ax[2, 2]], label="Probability", aspect=40)
            cbar.ax.set_yticklabels([]) 

        plt.pause(0.1)

        def plot_celerity_stats(axis_id, azimuth):
            pdf = self.eval_rcel_gmm(R, 1.0 / V, [azimuth] * len(R))
            axis_id.scatter(R, V, c=pdf, cmap=palette, marker=".", alpha=1.0, edgecolor='none', vmin=0.0, vmax=pdf_max)
            plt.pause(0.1)

        plot_celerity_stats(ax[0, 1], 0.0)
        plot_celerity_stats(ax[0, 2], 45.0)
        plot_celerity_stats(ax[1, 0], -90.0)
        plot_celerity_stats(ax[1, 2], 90.0)
        plot_celerity_stats(ax[2, 0], -135.0)
        plot_celerity_stats(ax[2, 1], 180.0)
        plot_celerity_stats(ax[2, 2], 135.0)

        if file_id:
            plt.savefig(file_id + "_cel-rng.png", bbox_inches='tight', dpi=200)

        if hold_fig:
            plt.show()

class TLossModel(object):
    _az_bin_cnt = 16

    def __init__(self):
        self.rng_vals = [0]
        self.tloss_vals = [0]

    def build(self, tloss_file, output_file, show_fits=False, use_coh=False, az_bin_cnt=16, az_bin_wdth=30.0, rng_lims=[1.0, 1000.0], rng_cnt=100, rng_smpls="linear"):
        """
            Construct propagation statistics from a NCPAprop modess or pape file (concatenated from
            multiple runs most likely) and output a transmission loss model

            Parameters
            ----------
            tloss_file: string
                Path to file containing NCPAprop transmission loss information
            output_file: string
                Path to file where results will be saved
            show_fits: boolean
                Option ot visualize model construction (for QC purposes)
            use_coh: boolean
                Option to use coherent transmission loss
            az_bin_cnt: int
                Number of azimuth bins to use in analysis
            az_bin_width: float
                Azimuth bin width in degrees for analysis
        """

        if os.path.isfile(output_file):
            print(output_file + " already exists  --->  Skipping transmission loss model construction...")
        else:
            print('Builing transmission loss models from file:', tloss_file)

            self._az_bin_cnt = az_bin_cnt

            # read in data, convert tloss to dB relative to 1 km, and wrap azimuths to [-180.0:180.0]
            print('\t' + "Reading in data...")
            if tloss_file[-3:] == ".nm":
                rngs, az, tloss_re, tloss_im, tloss_coh = np.loadtxt(tloss_file, unpack=True)
                if use_coh:
                    tloss = 10.0 * np.log10(tloss_coh)
                else:
                    tloss = 10.0 * np.log10(np.sqrt(tloss_re**2 + tloss_im**2))
            else:
                rngs, az, tloss_re, tloss_im = np.loadtxt(tloss_file, unpack=True)
                tloss = 10.0 * np.log10(np.sqrt(tloss_re**2 + tloss_im**2))

            # output_rngs = np.sort(np.unique(rngs)[::2])
            if rng_smpls == "linear":
                output_rngs = np.linspace(rng_lims[0], rng_lims[1], rng_cnt)
            else:
                output_rngs = np.logspace(np.log10(rng_lims[0]), np.log10(rng_lims[1]), rng_cnt)

            az[az > 180.0] -= 360.0
            az[az < -180.0] += 360.0
 
            tloss[np.isneginf(tloss)] = min(tloss[np.isfinite(tloss)])
            tloss[np.isposinf(tloss)] = max(tloss[np.isfinite(tloss)])

            tloss_vals = np.linspace(-150.0, 0.0, len(output_rngs))
            pdf_vals = np.empty((self._az_bin_cnt, len(output_rngs), len(tloss_vals)))

            for az_index in range(self._az_bin_cnt):
                center = -180 + 360.0 * (az_index / self._az_bin_cnt)
                if az_index == 0:
                    az_mask = np.logical_or(az >= 180.0 - az_bin_wdth / 2.0, az <= -180.0 + az_bin_wdth / 2.0)
                else:
                    az_mask = np.logical_and(center - az_bin_wdth / 2.0 <= az, az <= center + az_bin_wdth / 2.0)

                if show_fits:
                    f, ((ax1, ax2)) = plt.subplots(2, 1, figsize=(7.5, 10))

                    ax1.set_xlabel('Range [km]')
                    ax1.set_ylabel('Transmission Loss [dB]')
                    ax1.set_xlim([0.0, 1000.0])
                    ax1.set_ylim([min(tloss) - 5.0, max(tloss) + 5.0])

                    ax2.set_xlabel('Range [km]')
                    ax2.set_ylabel('Transmission Loss [dB]')
                    ax2.set_xlim([0.0, 1000.0])
                    ax2.set_ylim([min(tloss) - 5.0, max(tloss) + 5.0])

                    plt.suptitle("Stochastic Transmission Loss Model \n Azimuth: " + str(center), fontsize=18)
                    plt.show(block=False)

                    ax1.plot(rngs[az_mask][::11], tloss[az_mask][::11], 'ko', markersize=1)
                    plt.pause(0.001)

                print('\t' + "Building statistics for azimuth bin centered at " + str(center))

                norm_mask = np.logical_and(az_mask, rngs == min(rngs))
                if use_coh:
                    tloss_norm = 10.0 * np.log10(np.mean(tloss_coh[norm_mask]) * min(rngs))
                else:
                    tloss_norm = 10.0 * np.log10(np.mean(np.sqrt(tloss_re**2 + tloss_im**2)[norm_mask] * min(rngs)))

                # Define tloss pdf at each range point from KDE
                for nr, rng_val in enumerate(output_rngs):
                    nearest_rng = rngs[np.argmin(abs(rngs - rng_val))]
                    masked_tloss = tloss[np.logical_and(az_mask, rngs==nearest_rng)] - tloss_norm

                    if np.std(masked_tloss) < 0.05:
                        pdf_vals[az_index][nr] = norm.pdf(tloss_vals, loc=np.mean(masked_tloss), scale=0.05)
                    else:
                        masked_tloss_mean = np.mean(masked_tloss)
                        masked_tloss_stdev = np.std(masked_tloss)
                        outlier_mask = abs(masked_tloss - masked_tloss_mean) < 4.0 * masked_tloss_stdev

                        kernel = gaussian_kde(masked_tloss[outlier_mask])
                        pdf_vals[az_index][nr] = kernel.evaluate(tloss_vals)

                    TL1 = 10.0 * np.log(rng_val**(-1.0 / 4.0))
                    TL2 = 10.0 * np.log(rng_val**(-1.0 / 2.0))
                    env_stdev = max(abs(TL1 - TL2) / 10.0, 1.0)
                    envelope = 1.0 / (1.0 + np.exp((tloss_vals - (TL1 + 2.0 * env_stdev)) / env_stdev))
                    pdf_vals[az_index][nr] = pdf_vals[az_index][nr] * envelope

                    if show_fits:
                        ax2.scatter([rng_val] * len(tloss_vals), tloss_vals, c=pdf_vals[az_index][nr], cmap=cm.nipy_spectral_r, marker='o', s=[12.5] * len(tloss_vals), alpha=0.5, edgecolor='none')
                        plt.pause(0.001)

                if show_fits:
                    plt.close()

            priors = [0] * 3
            priors[0] = output_rngs
            priors[1] = tloss_vals
            priors[2] = pdf_vals

            pickle.dump(priors, open(output_file, "wb"))
            print(' ')

    def load(self, model_file, verbose=True):
        """
        Load a transmission loss file for use

        Parameters
        ----------

        model_file: string
            Path to TLoss file constructed using stochprop.propagation.TLossModel.build()

        """
        
        fit_params = pickle.load(open(model_file, "rb"), encoding='latin1')
        self.rng_vals = fit_params[0]
        self.tloss_vals = fit_params[1]
        self._az_bin_cnt = len(fit_params[2])

        if verbose:
            print("Loading transmission loss model from " + model_file)
            print('\t' + "Azimuth bin cound: " + str(self._az_bin_cnt))
            print('\t' + "Maximum range: " + str(max(self.rng_vals)) + '\n')

        self.pdf_vals = [0] * self._az_bin_cnt
        self.pdf_fits = [0] * self._az_bin_cnt

        for az_index in range(self._az_bin_cnt):
            self.pdf_vals[az_index] = fit_params[2][az_index]
            self.pdf_fits[az_index] = RectBivariateSpline(self.rng_vals, self.tloss_vals, self.pdf_vals[az_index])

    def eval(self, rng, tloss, az):
        """
            Evaluate TLoss model at specified range, transmission loss, and azimuth

            Parameters
            ----------
            rng: float
                Range from source
            tloss: float
                Transmission loss
            az: float
                Propagation azimuth (relative to North)

            Returns
            -------
            pdf: float
                Probability of observing an infrasonic arrival with specified transmission loss at specified range and azimuth


        """
        az_index = find_azimuth_bin(az, self._az_bin_cnt)

        in_rng_bnds = np.logical_and(self.rng_vals[0] <= rng, rng <= self.rng_vals[-1])
        in_tloss_bnds = np.logical_and(self.tloss_vals[0] <= tloss, tloss <= self.tloss_vals[-1])
        in_bnds = np.logical_and(in_rng_bnds, in_tloss_bnds)

        if len(np.atleast_1d(rng)) == 1:
            if in_bnds:
                result = self.pdf_fits[az_index].ev(rng, tloss)
            else:
                result = 0.0

        else:
            result = np.zeros_like(rng)
            for n_az in range(self._az_bin_cnt):
                mask = np.logical_and(in_bnds, az_index == n_az)
                if np.any(mask):
                    result[mask] = self.pdf_fits[n_az].ev(np.array(rng)[mask], np.array(tloss)[mask])

        return result

    def display(self, file_id=None, title="Transmission Loss Statistics", show_colorbar=True, hold_fig=False, show_ref_tloss=False):
        """
        Display the transmission loss statistics

        Parameters
        ----------
        file_id: string
            File prefix to save visualization
        subtitle: string
            Subtitle used in figures

        """
        scale_max = 0.1
        compass_file = find_spec('stochprop').submodule_search_locations[0] + '/resources/compass.png'

        resol = 100
        rngs = np.linspace(0.0, 1000.0, resol)
        tloss = np.linspace(-60.0, 0.0, resol)

        R, TL = np.meshgrid(rngs, tloss)
        R = R.flatten()
        TL = TL.flatten()

        palette = cm.nipy_spectral_r
        f1, ax = plt.subplots(3, 3, figsize=(14, 9))

        for n1, n2 in itertools.product(range(3), repeat=2):
            if n1 != 1 or n2 != 1:
                ax[n1, n2].set_ylim([tloss[0], tloss[-1]])
                ax[n1, n2].set_xlim([rngs[0], rngs[-1]])
                ax[n1, n2].set_xticks([0, 250, 500, 750, 1000])
            if n2 != 0:
                ax[n1, n2].set_yticklabels([])
            if n1 != 2:
                ax[n1, n2].set_xticklabels([])

        img = mpimg.imread(compass_file)
        ax[1, 1].axis('off')
        ax[1, 1].imshow(img)

        ax[2, 1].set_xlabel('Range [km]')
        ax[1, 0].set_ylabel('Transmission Loss [dB]')

        if title:
            plt.suptitle(title, fontsize=22)
        plt.show(block=False)

        pdf = self.eval(R, TL, np.array([-45.0] * len(R)))
        tloss_plot = ax[0, 0].scatter(R, TL, c=pdf, cmap=palette, marker='o', s=[12.5] * len(R), alpha=0.5, edgecolor='none', vmin=0.0, vmax=scale_max)
        if show_ref_tloss:
            ax[0, 0].plot(rngs[rngs > 0], 10.0 * np.log10(1.0 / np.sqrt(rngs[rngs > 0])), '--k')
            ax[0, 0].plot(rngs[rngs > 0], 10.0 * np.log10(1.0 / rngs[rngs > 0]), '-k')

        if show_colorbar:
            cbar = f1.colorbar(tloss_plot, ax=[ax[0, 2], ax[1, 2], ax[2, 2]], label="Probability", aspect=40)
            cbar.ax.set_yticklabels([]) 

        plt.pause(0.01)

        def plot_tloss_stats(axis_id, azimuth):
            pdf = self.eval(R, TL, np.array([azimuth] * len(R)))
            axis_id.scatter(R, TL, c=pdf, cmap=palette, marker='o', s=[12.5] * len(R), alpha=0.5, edgecolor='none', vmin=0.0, vmax=scale_max)
            if show_ref_tloss:
                axis_id.plot(rngs[rngs > 0], 10.0 * np.log10(1.0 / np.sqrt(rngs[rngs > 0])), '--k')
                axis_id.plot(rngs[rngs > 0], 10.0 * np.log10(1.0 / rngs[rngs > 0]), '-k')
            plt.pause(0.01)

        plot_tloss_stats(ax[0, 1], 0.0)
        plot_tloss_stats(ax[0, 2], 45.0)
        plot_tloss_stats(ax[1, 0], -90.0)
        plot_tloss_stats(ax[1, 2], 90.0)
        plot_tloss_stats(ax[2, 0], -135.0)
        plot_tloss_stats(ax[2, 1], 180.0)
        plot_tloss_stats(ax[2, 2], 135.0)

        if file_id:
            plt.savefig(file_id + "_tloss.png", bbox_inches='tight', dpi=200)

        if hold_fig:
            plt.show()


#########################
##   IMS Noise Model   ## 
#########################
def ims_noise_model():

    ims_ns = np.loadtxt(find_spec('stochprop').submodule_search_locations[0] + '/resources/IMSNOISE_MIN_MED_MAX.txt')

    ns_mean = interp1d(ims_ns[:, 0], 10.0 * ims_ns[:, 2], bounds_error=False, fill_value="extrapolate")
    ns_sigma = interp1d(ims_ns[:, 0], 10.0 * (ims_ns[:, 3] - ims_ns[:, 1]) / 4.0, bounds_error=False, fill_value="extrapolate")

    def ims_ns_pdf(P, f, duration=10.0, array_dim=4):
        eff_ns_mean = 0.5 * (ns_mean(f) + 10.0 * np.log10(duration / array_dim))
        return norm.pdf(P, loc=eff_ns_mean, scale=ns_sigma(f))

    def ims_ns_cdf(P, f, duration=10.0, array_dim=4):
        eff_ns_mean = 0.5 * (ns_mean(f) + 10.0 * np.log10(duration / array_dim))
        return norm.cdf(P, loc=eff_ns_mean, scale=ns_sigma(f))

    return ims_ns_pdf, ims_ns_cdf



#########################
##   Kinney & Graham   ## 
##   Blastwave Model   ##
#########################
def kg_op(W, r, p_amb=101.325, T_amb=288.15, exp_type="chemical"):
    """
        Kinney & Graham scaling law peak overpressure model
                
        Parameters
        ----------
        W: float
            Explosive yield of the source [kg eq. TNT]
        r: float
            Propagation distance [km]
        p_amb: float
            Ambient atmospheric pressure [kPa]
        T_amb: float
            Ambient atmospheric temperature [deg K]
        exp_type: string
            Type of explosion modeled, options are "chemical" or "nuclear"
        
        Returns
        -------
        p0: float
            Peak overpressure [Pa]    
    """
    
    fd = (p_amb / 101.325)**(1.0 / 3.0) * (T_amb / 288.15)**(1.0 / 3.0)
    sc_rng = fd / W**(1.0 / 3.0) * r * 1000.0
    
    if exp_type=="chemical":
        term1 = 1.0 + (sc_rng / 4.5)**2
        term2 = np.sqrt(1.0 + (sc_rng / 0.048)**2)
        term3 = np.sqrt(1.0 + (sc_rng / 0.32)**2)
        term4 = np.sqrt(1.0 + (sc_rng / 1.35)**2)
        
        result = 808.0 * term1 / (term2 * term3 * term4)
    else:
        term1 = (1.0 + sc_rng / 800.0)
        term2 = np.sqrt(1.0 + (sc_rng / 87.0)**2)
        
        result = 3.2e6 / sc_rng**3 * term1 * term2

    return p_amb * 1.0e3 * result


def kg_ppd(W, r, p_amb=101.325, T_amb=288.15, exp_type="chemical"):
    """
        Kinney & Graham scaling law positive phase duration model
                
        Parameters
        ----------
        W : float
            Explosive yield of the source [kg eq. TNT]
        r : float
            Propagation distance [km]
        p_amb : float
            Ambient atmospheric pressure [kPa]
        T_amb : float
            Ambient atmospheric temperature [deg K]
        exp_type : string
            Type of explosion modeled, options are "chemical" or "nuclear"
            
        Returns
        -------
        t0 : float
            Positive phase duration [s]
    """

    fd = (p_amb / 101.325)**(1.0 / 3.0) * (T_amb / 288.15)**(1.0 / 3.0)
    sc_rng = fd / W**(1.0 / 3.0) * r * 1000.0
    
    if exp_type=="chemical":
        term1 = 1.0 + (sc_rng / 0.54)**10
        term2 = 1.0 + (sc_rng / 0.02)**3
        term3 = 1.0 + (sc_rng / 0.74)**6
        term4 = np.sqrt(1.0 + (sc_rng / 6.9)**2)
        
        result = 980.0 * term1 / (term2 * term3 * term4)
    else:
        term1 = np.sqrt(1.0 + (sc_rng / 100.0)**3)
        term2 = np.sqrt(1.0 + (sc_rng / 40.0))
        term3 = (1.0 + (sc_rng / 285.0)**5)**(1.0 / 6.0)
        term4 = (1.0 + (sc_rng / 50000.0))**(1.0 / 6.0)
        
        result = 180.0 * term1 / (term2 * term3 * term4)

    return result * W**(1.0 / 3.0) / 1e3


def blastwave_spectrum(f, W, r, p_amb=101.325, T_amb=288.15, exp_type="chemical", shaping_param=0.0):
    """ 
        Fourier transform amplitude for the acoustic
            blastwave in sasm.acoustic.blastwave().
        
        Note: shaping_param = 0 produces the Friedlander
        blastwave model.
        
        Note: the peak of the spectrum occurs at
        f_0 = \frac{1}{2 \pi t_0} \frac{1}{\sqrt{\alpha + 1}
        and t0 corresponding to a given peak frequency is
        t_0 = \frac{1}{2 \pi f_0} \frac{1}{\sqrt{\alpha + 1}
        
        Parameters
        ----------
        f: float
            Frequency [Hz]
        W: float
            Explosive yield of the source [kg eq. TNT]
        r: float
            Propagation distance [km]
        p_amb: float
            Ambient atmospheric pressure [kPa]
        T_amb: float
            Ambient atmospheric temperature [deg K]
        exp_type: string
            Type of explosion modeled, options are "chemical" or "nuclear"
        shaping_param: float
            Shaping parameter for waveform that controls high frequency trend (set to 0 for original Friedlander blastwave)
        
        Returns
        -------
        P : float
            Spectral value at frequency f
    """

    p0 = kg_op(W, r, p_amb, T_amb, exp_type)
    t0 = kg_ppd(W, r, p_amb, T_amb, exp_type)

    omega = 2.0 * np.pi * f
    x0 = (1.0 + shaping_param) - np.sqrt(1.0 + shaping_param)
    norm = x0**shaping_param * (1.0 - x0 / (1.0 + shaping_param)) * np.exp(-x0)
    
    return p0 / norm * t0 * gamma(shaping_param + 1.0) * (omega * t0) / (1.0 + (omega**2 * t0**2))**(shaping_param / 2.0 + 1.0)


# ############################# #
#      Detection Capbility      #
#     & Network Performance     #
# ############################# #
def plot_detection_stats(tlms, yld_vals, array_dim, output_path=None, show_fig=True):

    plt.rcParams.update({'font.size': 14})
    _, ims_ns_cdf = ims_noise_model()

    P_vals = np.linspace(-50.0, max(10.0 * np.log10(blastwave_spectrum(np.array(tlms[0]), 2.0 * max(yld_vals), 1.0))) + 10.0, 200)

    fig_size = (3 * len(tlms[0]) + 1, 3 * len(yld_vals))
    fig, axes = plt.subplots(len(yld_vals), len(tlms[0]) + 1, figsize=fig_size, subplot_kw={'projection': 'polar'})
    for m in range(len(yld_vals)):
        print('\n' + "Computing detection statistics for surface explosion with yield " + str(yld_vals[m] * 1.0e-3) + " tons eq TNT...")
        for n in range(len(tlms[0])):
            print('\t' + "Analyzing TLM at " + str(tlms[0][n]) + " Hz...")
            tlm_freq = tlms[0][n]
            tlm = tlms[1][n]
            src0 = 10.0 * np.log10(blastwave_spectrum(tlm_freq, 2.0 * yld_vals[m], 1.0))
   
            rng_vals = np.linspace(5.0, tlm.rng_vals[-1], 100)
            dr = abs(rng_vals[1] - rng_vals[0])

            P_grid, rng_grid = np.meshgrid(P_vals, rng_vals)

            duration = rng_grid * (0.34 - 0.27) / (0.34 * 0.27)
            duration[duration < 5.0] = 5.0

            if len(yld_vals) > 1 and len(tlms[0]) > 1:
                ax_j = axes[m][n]
            elif len(yld_vals) > 1:
                ax_j = axes[m]
            else:
                ax_j = axes[n]

            ax_j.set_theta_zero_location("N")
            ax_j.set_theta_direction(-1)

            ax_j.set_xticks(np.linspace(0, 2 * np.pi, 4, endpoint=False))
            ax_j.set_yticks(np.linspace(0, np.round(tlm.rng_vals[-1], -1), 5))

            if n == 0 and m == 0:
                ax_j.set_xticklabels(['N', 'E', 'S', 'W'])
            else:
                ax_j.set_xticklabels([])
                ax_j.set_yticklabels([])

            [i.set_color('black') for i in ax_j.get_yticklabels()]
            if m == 0:
                ax_j.annotate(str(tlm_freq) + " Hz", (0.075, 0.925), xycoords='axes fraction', horizontalalignment='center')
            if n == 0:
                ax_j.annotate(str(yld_vals[m] * 1.0e-3) + " tons", (0.025, 0.075), xycoords='axes fraction', horizontalalignment='center')

            for az_index in range(tlm._az_bin_cnt):
                center = -180 + 360.0 * (az_index / tlm._az_bin_cnt)
                arrival = tlm.eval(rng_grid, P_grid - src0, np.ones_like(rng_grid) * center)
                rng_prob = simps(arrival * ims_ns_cdf(P_grid, tlm_freq, duration=duration, array_dim=array_dim), P_vals)

                im = ax_j.bar(np.radians(np.array([center] * len(rng_vals))), width=np.radians(360.0 / tlm._az_bin_cnt), height=dr * 2.0, bottom=rng_vals, color=cm.hot_r(rng_prob), alpha=0.9)
                        
    print('\n' + "Plotting detection statistics...", '\n')

    if len(yld_vals) > 1:
        for m in range(len(yld_vals)):
            axes[m][-1].axis('off')
    else:
            axes[-1].axis('off')


    cax = fig.add_axes([(len(tlms[0]) + 0.1) / (len(tlms[0]) + 1.0), 0.1, 0.02, 0.8])
    cbar = ColorbarBase(cax, cmap=cm.hot_r)
    cbar.set_label("Detection Probability")


    plt.tight_layout()

    if output_path:
        plt.savefig(output_path + ".det-stats.png", dpi=300) 
    
    if show_fig:
        plt.show()


def plot_network_performance(info_file, freq, W0, det_cnt_min, lat_min, lat_max, lon_min, lon_max, resol, output_path, show_fig):

    # Extract network locations and unique TLM/dim info
    print("Loading network info from " + info_file + "..." + '\n')
    sta_locs, dims, model_files = [], [], []
    with open(info_file, 'r') as of:
        for line in of:
            temp = line.strip('()[]\n').replace(' ', '').split(',')
            sta_locs = sta_locs + [[float(item) for item in temp[:2]]]
            dims = dims + [int(temp[2])]
            model_files = model_files + [temp[3]]

    sta_locs = np.array(sta_locs)
    dims = np.array(dims)

    # Loop over unique TLM/dim combinations, build PDF(r, az) interpolation, and evaluate on the grid
    _, ims_ns_cdf = ims_noise_model()

    # check if frequency value is in TLM file name(s):
    try:
        freq = float(model_files[0].split("_")[-1].split("Hz")[0])
    except:
        print("Warning!  Couldn't extract frequency from TLM file name(s).  Using specified value (" + freq + ")")
        freq = float(freq)

    print("Computing spectral amplitude for source:")
    print('\t' + "Yield: " + str(W0 * 1.0e-3) + " ton eq. TNT")
    print('\t' + "Frequency: " + str(freq) + " Hz")

    src0 = 10.0 * np.log10(blastwave_spectrum(freq, 2.0 * W0, 1.0))
    P_vals = np.linspace(-50.0, src0 + 10.0, 200)

    lat_vals = np.linspace(lat_min, lat_max, resol)
    lon_vals = np.linspace(lon_min, lon_max, resol)

    grid_lats, grid_lons = np.meshgrid(lat_vals, lon_vals)
    grid_lats = grid_lats.flatten()
    grid_lons = grid_lons.flatten()

    print('\n' + "Evaluating detection statistics for each station...")
    det_stats = []
    for j, model_info in enumerate(model_files):
        print('\t' + "Loading TLM from '" + model_info + "' for station at (" + str(sta_locs[j][0]) + ", " + str(sta_locs[j][1]) + ")")

        tlm = TLossModel()
        tlm.load(model_info, verbose=False)

        rng_vals = np.linspace(1.0, tlm.rng_vals[-1], 100)
        P_grid, rng_grid = np.meshgrid(P_vals, rng_vals)

        duration = rng_grid * (0.34 - 0.27) / (0.34 * 0.27)
        duration[duration < 2.5] = 2.5
        
        det_prob = np.empty((len(rng_vals), tlm._az_bin_cnt))
        for az_index in range(tlm._az_bin_cnt):
            center = -180 + 360.0 * (az_index / tlm._az_bin_cnt)
            arrival = tlm.eval(rng_grid, P_grid - src0, np.ones_like(rng_grid) * center)
            det_prob[:, az_index] = simps(arrival * ims_ns_cdf(P_grid, freq, duration=duration, array_dim=dims[j]), P_vals, axis=1)

        az_vals = np.array([-180 + 360.0 * (j / tlm._az_bin_cnt) for j in range(tlm._az_bin_cnt)])
        az_vals = np.append(az_vals, 180.0)
        det_prob = np.column_stack((det_prob, det_prob[:, 0]))

        interp = RectBivariateSpline(rng_vals, az_vals, det_prob)

        temp = np.array(sph_proj.inv(sta_locs[j][1] * np.ones_like(grid_lons), sta_locs[j][0] * np.ones_like(grid_lats), grid_lons, grid_lats, radians=False))
        output_azs, output_rngs = np.array(temp[0]) - 180.0, np.array(temp[2]) / 1000.0

        output_azs[output_azs > 180.0] -= 360.0
        output_azs[output_azs < -180.0] += 360.0

        output_rngs[output_rngs < 5.0] = 5.0

        det_stats = det_stats + [interp.ev(output_rngs, output_azs)]
   
    det_stats = np.array(det_stats)

    # Combine statistics and generate map
    print('\n' + "Combining station detection stats to determine network performance...")
    print('\t' + "Requiring " + str(det_cnt_min) + " detecting stations.")
    result = 1.0 - np.prod(1.0 - det_stats, axis=0)
    for n in range(1, det_cnt_min):
        for indices in itertools.combinations(list(range(len(sta_locs))), r=n):
            P1 = np.array([det_stats[j] for j in range(len(det_stats)) if j in indices])
            P2 = np.array([(1.0 - det_stats[j]) for j in range(len(det_stats)) if j not in indices])

            result = result - np.prod(P1, axis=0) * np.prod(P2, axis=0)


    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=map_proj)

    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)

    gl = ax.gridlines(crs=map_proj, draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    lat_tick, lon_tick = int((lat_max - lat_min) / 5), int((lon_max - lon_min) / 5)
    gl.xlocator = mticker.FixedLocator(np.arange(lon_min - np.ceil(lon_tick / 2), lon_max + lon_tick, lon_tick))
    gl.ylocator = mticker.FixedLocator(np.arange(lat_min - np.ceil(lat_tick / 2), lat_max + lat_tick, lat_tick))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Add features (coast lines, borders)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    if (lon_max - lon_min) < 20.0:
        ax.add_feature(cfeature.STATES, linewidth=0.5)
        ax.add_feature(cfeature.RIVERS, edgecolor='dodgerblue', alpha=0.3)
        ax.add_feature(cfeature.LAKES, facecolor='dodgerblue', edgecolor='dodgerblue', alpha=0.3)

    sc = ax.scatter(grid_lons, grid_lats, c=result, cmap=cm.YlOrRd, vmin=0.0, vmax=1.0)
    ax.plot(sta_locs[:, 1], sta_locs[:, 0], "^k")

    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar = plt.colorbar(sc, cax=ax_cb)
    cbar.set_label('Detection Probability')
    
    if output_path:
        plt.savefig(output_path + ".network.png", dpi=300) 
    
    if show_fig:
        plt.show()


def yield_hob_stats(yld_vals, hob_vals, infraga_config, output_path, channel_Cnt, local_temp_dir=None):

    _, ims_ns_cdf = ims_noise_model()

    freq_bands = [[0.02, 0.1],
                [0.05, 0.3],
                [0.2, 0.8],
                [0.5, 2.0],
                [1.0, 8.0]]


    with tempfile.TemporaryDirectory(prefix='infraga_') as tmpdirname:
        if local_temp_dir is not None:
            print("\n  Writing individual infraGA results into local directory: " + local_temp_dir + '\n')
            if not os.path.isdir(local_temp_dir):
                os.mkdir(local_temp_dir)
            tmpdirname = local_temp_dir
        else:
            print("\n  Created temp directory:", tmpdirname + '\n')
        
        for alt in hob_vals:
            if os.path.isfile(tmpdirname + "/" + str(np.round(alt, 0)) + "km.arrivals.dat"):
                print(tmpdirname + "/" + str(np.round(alt, 0)) + "km.arrivals.dat already exists.  Skipping infraGA run...")
            else:
                for yld in yld_vals:
                    command = "infraga sph eig_wvfrm --config-file " + infraga_config + " --src-alt " + str(alt) + " --wvfrm-yield " + str(np.round(yld, 2))
                    command = command + " --keep-eig-results true --output-id " + tmpdirname + "/" + str(np.round(alt, 0)) + "km"
                    print('\n' + "-" * 40 + '\n' + command)
                    subprocess.call(command, shell=True)

                    command = "mv " + tmpdirname + "/" + str(np.round(alt, 0)) + "km.wvfrms.dat "
                    command = command + tmpdirname + "/" + str(np.round(alt, 0)) + "km_" + str(np.round(yld / 1.0e3, 2)) + "t.wvfrms.dat"
                    print(command)
                    subprocess.call(command, shell=True)

        print('\nBuilding statistics...')
        file_out = open(output_path + ".det_stats.dat", 'w')


        infraga_cnfg = cnfg.ConfigParser()
        infraga_cnfg.read(infraga_config)

        print("# stochprop prop yld-hob results (stochprop v" + pkg_resources.get_distribution("stochprop").version + ")", file=file_out)
        print("#   infraga config: " + infraga_config, file=file_out)
        print("#   atmo file: " + infraga_cnfg['GENERAL']['atmo_file'], file=file_out)
        print("#   source lat/lon: " + infraga_cnfg['EIGENRAY']['src_lat'] + ", " + infraga_cnfg['EIGENRAY']['src_lon'], file=file_out) 
        print("#   receiver lat/lon: " + infraga_cnfg['EIGENRAY']['rcvr_lat'] + ", " + infraga_cnfg['EIGENRAY']['rcvr_lon'], file=file_out) 
        print('#\n# alt [km]\tyld [t eq. TNT]\tdet_prob\tdom_freq [Hz]\tpk_spec_amp [Pa/Hz]', file=file_out)

        for alt in hob_vals:
            print('\tAnalyzing waveforms for HOB: ' + str(alt) + " km...")
            for yld in yld_vals:
                wvfrm_file = tmpdirname + "/" + str(np.round(alt, 0)) + "km_" + str(np.round(yld / 1.0e3, 2)) + "t.wvfrms.dat"
                if os.path.isfile(wvfrm_file):
                    wvfrm_data = np.loadtxt(wvfrm_file)
                    t_vals, p_vals = wvfrm_data[:, 0], np.sum(wvfrm_data[:, 1:], axis=1)            

                    dt, N = abs(t_vals[1] - t_vals[0]), len(t_vals)
                    N_f = int(N / 2 + 1)

                    f = (1.0 / dt) * (np.arange(float(N_f)) / N)
                    p_fft = abs(np.fft.rfft(p_vals)) * dt

                    # extract dominant preiod and peak spectral amplitude
                    f_dom = f[np.argmax(p_fft)]
                    p_fft_pk = np.max(p_fft)

                    # Compute detection probability
                    duration = (t_vals[-1] - t_vals[0]) / 2.0
                    det_prob = []
                    for band in freq_bands:
                        f_mask = np.logical_and(band[0] < f, f < band[1])
                        det_prob = det_prob + [np.nanmean(ims_ns_cdf(f[f_mask], p_fft[f_mask], duration=duration, array_dim=channel_Cnt))]

                    det_prob = np.array(det_prob)
                    print(str(int(alt)) + '\t' + str(int(np.round(yld / 1.0e3))) + '\t' + str(np.nanmax(det_prob)) + '\t' + str(f_dom) + '\t' + str(p_fft_pk), file=file_out)
                else:
                    print('\t\tError loading waveform result: ' + wvfrm_file + ' (eigenrays may not exist or simulations might have failed).')



