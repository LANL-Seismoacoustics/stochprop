# propagation.py
#
# Run infraga-accel-3d or NCPAprop's modess to generate predictions to
# use for constructing stochastic propagation models
#
# Philip Blom (pblom@lanl.gov)

import sys
import fnmatch
import warnings
import os
import pickle
import imp
import itertools

import numpy as np

from pyproj import Geod

from scipy.interpolate import interp1d, interp2d, RectBivariateSpline
from scipy.optimize import curve_fit
from scipy.stats import norm, gaussian_kde
from scipy.signal import savgol_filter

import matplotlib.cm as cm
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

sph_proj = Geod(ellps='sphere')

##################################
#  Running infraga and NCPAprop  #
#     with multiple profiles     #
##################################
def run_infraga(profs_path, results_file, pattern="*.met", cpu_cnt=None, geom="3d", bounces=25, inclinations=[0.5, 45.0, 1.5], azimuths=[-180.0, 180.0, 6.0], freq=0.1, z_grnd=0.0, rng_max=1000.0, src_loc=[0.0, 0.0, 0.0], infraga_path=""):
    """
        Run the infraga -prop algorithm to compute path geometry
            statistics for BISL using a suite of specifications 
            and combining results into single file

        Parameters
        ----------
        profs_path : string
            Path to atmospheric specification files
        results_file : string
            Path and name of file where results will be written
        pattern : string
            Pattern identifying atmospheric specification within profs_path location
        cpu_cnt : int
            Number of threads to use in OpenMPI implementation.  None runs non-OpenMPI version of infraga
        geom : string
            Defines geometry of the infraga simulations ("2d", "3d", or "sph")
        bounces : int
            Maximum number of ground reflections to consider in ray tracing
        inclinations : iterable object
            Iterable of starting, ending, and step for ray launch inclination
        azimuths : iterable object
            Iterable of starting, ending, and step for ray launch azimuths
        freq : float
            Frequency to use for Sutherland Bass absorption calculation
        z_grnd : float
            Elevation of the ground surface relative to sea level
        rng_max : float
            Maximum propagation range for propagation paths
        src_loc : iterable object
            The horizontal (latitude and longitude) and altitude of the source
        """

    if os.path.isfile(results_file):
        print(results_file + " already exists  --->  Skipping infraGA/GeoAc ray tracing runs...")
    else:
        if geom is not ("3d" or "sph"):
            msg = "Incompatible geometry option for infraga: {}.  Options are 3d' and 'sph'".format(geom)
            warnings.warn(msg)

        open(results_file, 'w').close()

        dir_files = os.listdir(profs_path)
        for file_name in dir_files:
            if fnmatch.fnmatch(file_name, pattern):
                print("Generating ray paths for " + file_name)
                file_id = os.path.splitext(file_name)[0]

                if cpu_cnt:
                    command = "mpirun -np " + str(cpu_cnt) + " " + infraga_path + " infraga-accel-" + geom + " -prop "
                else:
                    command = infraga_path + " infraga-" + geom + " -prop "
            
                command = command + profs_path + "/" + file_name
                command = command + " incl_min=" + str(inclinations[0]) + " incl_max=" + str(inclinations[1]) + " incl_step=" + str(inclinations[2])
                command = command + " az_min=" + str(azimuths[0]) + " az_max=" + str(azimuths[1]) + " az_step=" + str(azimuths[2])
                if geom == "sph":
                    command = command + " src_lat=" + str(src_loc[0]) + " src_lon=" + str(src_loc[1])
                command = command + " src_alt=" + str(src_loc[2])
                command = command + " freq=" + str(freq) + " z_grnd=" + str(z_grnd) + " max_rng=" + str(rng_max)
                command = command + " calc_amp=False" + " bounces=" + str(bounces) + " write_rays=false" 

                print(command)
                os.system(command)
    
        command = "cat " + profs_path + "/*.arrivals.dat > " + results_file
        print(command)
        os.system(command)
    
        command = "rm "  + profs_path + "/*.dat"
        print(command)
        os.system(command)


def run_modess(profs_path, results_path, pattern="*.met", cpu_cnt=None, azimuths=[-180.0, 180.0, 6.0], freqs=[0.1, 1.0, 10], z_grnd=0.0, rng_max=1000.0, ncpaprop_path=""):
    """
        Run the NCPAprop normal mode methods to compute transmission
            loss values for a suite of atmospheric specifications at
            a set of frequency values

        Parameters
        ----------
        profs_path : string
            Path to atmospheric specification files
        results_file : string
            Path and name of file where results will be written
        pattern : string
            Pattern identifying atmospheric specification within profs_path location
        azimuths : iterable object
            Iterable of starting, ending, and step for propagation azimuths
        freqs : float
            Iterable of starting, ending, and step for frequencies
        z_grnd : float
            Elevation of the ground surface relative to sea level
        rng_max : float
            Maximum propagation range for propagation paths
        """

    dir_files = os.listdir(profs_path)
    for fn in np.logspace(np.log10(freqs[0]), np.log10(freqs[1]), freqs[2]):
        if os.path.isfile(results_path + "_%.3f" %fn + ".lossless.nm"):
            print(results_path + "_%.3f" %fn + ".lossless.nm alraedy exists  --->  Skipping NCPAprop modess runs...")
        else:
            for file_name in dir_files:
                if fnmatch.fnmatch(file_name, pattern):
                    file_id = os.path.splitext(file_name)[0]

                    print('\t' + "Running NCPAprop modess for " + profs_path + "/" + file_name + " at " + "%.3f" % fn + " Hz...")
                    command = ncpaprop_path + "Modess --atmosfile " + profs_path + "/" + file_name + " --atmosfileorder ztuvdp --skiplines 0 --freq " + str(fn)
                    command = command + " --Nby2Dprop --azimuth_start " + str(azimuths[0]) + " --azimuth_end " + str(azimuths[1]) + " --azimuth_step " + str(azimuths[2])
                    command = command + "  --maxrange_km " + str(rng_max) + "--zground_km " + str(z_grnd) + " > /dev/null"
                    os.system(command)

                    os.system("mv Nby2D_tloss_1d.lossless.nm " + profs_path + "/" + file_id + "_%.3f" % fn + "Hz.lossless.nm")
                    os.system("mv Nby2D_tloss_1d.nm " + profs_path + "/" + file_id + "_%.3f" % fn + "Hz.nm")

            print('\t' + "Combining transmission loss predictions..." + '\n')
            os.system("cat " + profs_path + "/*_%.3f" % fn + "Hz.lossless.nm > " + results_path + "_%.3f" %fn + ".lossless.nm")
            os.system("rm " + profs_path + "/*_%.3f" % fn + "Hz.lossless.nm")

            os.system("cat " + profs_path + "/*_%.3f" % fn + "Hz.nm > " + results_path + "_%.3f" %fn + ".nm")
            os.system("rm " + profs_path + "/*_%.3f" % fn + "Hz.nm")


# ############################ #
#          Stochastic          #
#      Propagation Models      #
# ############################ #
az_dirs = ['S', 'SW', 'W', 'NW', 'N', 'NE', 'E', 'SE']

def find_azimuth_bin(az, bin_cnt=8):
    reduced = np.degrees(np.arctan2(np.sin(np.radians(az)), np.cos(np.radians(az))))
    bins = np.arange(-180.0, 180.0, 360.0 / (bin_cnt * 2.0))

    result = np.asarray(np.digitize(reduced, bins) / 2)
    result[result >= bin_cnt] = 0
    return result.astype(int)


class PathGeometryModel(object):
    az_bin_cnt, az_bin_wdth = 8, 60.0

    rng_max = 1000.0

    default_az_dev_std = 4.0
    min_az_dev_std = 0.9

    tropo_strat_bnd = 1.0 / 0.32
    strat_therm_bnd = 1.0 / 0.26
    bnd_overlap = 0.05

    wts0 = np.array([0.0539, 0.0899, 0.8562])
    std0 = np.array([0.066, 0.08, 0.33])
    mns0 = np.array([1.0 / 0.33, 1.0 / 0.29, 1.0 / 0.26])
    
    win_min_pts = 25
    rcel_std_min = 0.05
    rcel_wt_min = 1.0e-2

    def __init__(self):
        self.rngs = np.array([])

        self.rcel_wts = []
        self.rcel_mns = []
        self.rcel_std = []

        self.az_dev_mns = []
        self.az_dev_std = []


    def eval_rcel_gmm(self, rng, rcel, az):
        rng_eval = np.array(rng)
        rng_eval[rng_eval > self.rng_max] = self.rng_max

        if len(np.atleast_1d(rng)) == 1:
            n_az = find_azimuth_bin(az, self.az_bin_cnt)
            fit_rcel_wts = np.array([f(rng_eval) for f in self.rcel_wts[n_az]])
            fit_rcel_mns = np.array([f(rng_eval) for f in self.rcel_mns[n_az]])
            fit_rcel_std = np.array([f(rng_eval) for f in self.rcel_std[n_az]])
            result = np.sum(fit_rcel_wts / fit_rcel_std * norm.pdf((rcel - fit_rcel_mns) / fit_rcel_std))
        else:
            mn = np.empty((len(rng_eval), 3))
            vr = np.empty((len(rng_eval), 3))
            wt = np.empty((len(rng_eval), 3))

            az_indices = find_azimuth_bin(az, self.az_bin_cnt)
            for n_az in range(self.az_bin_cnt):
                mask = tuple([az_indices==n_az])
                if np.any(mask):
                    mn[mask] = np.array([f(rng_eval[mask]) for f in self.rcel_mns[n_az]]).T
                    vr[mask] = np.array([f(rng_eval[mask]) for f in self.rcel_std[n_az]]).T
                    wt[mask] = np.array([f(rng_eval[mask]) for f in self.rcel_wts[n_az]]).T

            result = np.sum(wt / vr * norm.pdf((np.array([rcel] * 3).T - mn) / vr), axis=1)
        return result


    def eval_az_dev_mn(self, rng, az):
        if len(np.atleast_1d(rng)) == 1:
            return self.az_dev_mns[find_azimuth_bin(az, self.az_bin_cnt)](min(rng, self.rng_max))
        else:
            rng_eval = np.array(rng)
            rng_eval[rng_eval > self.rng_max] = self.rng_max
            az_indices = find_azimuth_bin(az, self.az_bin_cnt)

            mn = np.empty_like(rng_eval)
            for n_az in range(self.az_bin_cnt):
                mask = az_indices==n_az
                if np.any(mask):
                    mn[mask] = self.az_dev_mns[n_az](rng_eval[mask])
            return mn

    def eval_az_dev_vr(self, rng, az):
        if len(np.atleast_1d(rng)) == 1:
            return self.az_dev_std[find_azimuth_bin(az, self.az_bin_cnt)](min(rng, self.rng_max))
        else:
            rng_eval = np.array(rng)
            rng_eval[rng_eval > self.rng_max] = self.rng_max
            az_indices = find_azimuth_bin(az, self.az_bin_cnt)

            vr = np.empty_like(rng_eval)
            for n_az in range(self.az_bin_cnt):
                mask = az_indices==n_az
                if np.any(mask):
                    vr[mask] = self.az_dev_std[n_az](rng_eval[mask])
            return vr

    def build(self, arrivals_file, output_file, show_fits=False, rng_width=50.0, rng_spacing=10.0, geom="3d", src_loc=[0.0, 0.0, 0.0]):
        if os.path.isfile(output_file):
            print(output_file + " already exists  --->  Skipping path geometry model construction...")
        else:
            if geom is not ("3d" or "sph"):
                msg = "Incompatible geometry option for infraga: {}.  Options are 3d' and 'sph'".format(geom)
                warnings.warn(msg)
            
            print('Builing celerity and azimuth priors from file:', arrivals_file)

            # define range bins and parameter arrays
            rng_bins = np.arange(0.0, self.rng_max, rng_spacing)
            rng_cnt = len(rng_bins)

            az_dev_mns = np.empty((self.az_bin_cnt, rng_cnt))
            az_dev_std = np.empty((self.az_bin_cnt, rng_cnt))

            rcel_wts = np.empty((self.az_bin_cnt, rng_cnt, 3))
            rcel_mns = np.empty((self.az_bin_cnt, rng_cnt, 3))
            rcel_std = np.empty((self.az_bin_cnt, rng_cnt, 3))

            # load infraGA/GeoAc predictions
            if geom == "3d":
                theta, phi, n, x, y, t, cel, z_max, incl, back_az, amp_geo, amp_atmo = np.loadtxt(arrivals_file, unpack=True)

                rngs = np.sqrt(x**2 + y**2)
                az = 90.0 - np.degrees(np.arctan2(y, x))
                az_dev = (90.0 - np.degrees(np.arctan2(-y, -x))) - back_az
            else:
                theta, phi, n, lats, lons, t, cel, z_max, incl, back_az, amp_geo, amp_atmo = np.loadtxt(arrivals_file, unpack=True)

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

            # Cycle through azimuth bins creating fit
            az_wts = np.empty(self.az_bin_cnt)
            for n_az in range(self.az_bin_cnt):
                if n_az == 0:
                    az_mask = np.logical_or(az > 180.0 - self.az_bin_wdth / 2.0, az < -180.0 + self.az_bin_wdth / 2.0)
                else:
                    center = -180 + (360.0 / self.az_bin_cnt) * n_az
                    az_mask = np.logical_and(center - self.az_bin_wdth / 2.0 <= az, az <= center + self.az_bin_wdth / 2.0)
                az_wts[n_az] = float(len(rngs[az_mask])) / float(len(rngs))

                # display work if show_fits is enabled
                if show_fits:
                    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 11))

                    ax1.set_xlim([0.0, self.rng_max])
                    ax1.set_ylim([0.2, 0.4])

                    ax3.set_xlim([0.0, self.rng_max])
                    ax3.set_ylim([-12.0, 12.0])

                    ax1.set_xlabel('Range [km]')
                    ax1.set_ylabel('Celerity [km/s]')

                    ax2.set_xlabel('Celerity [km/s]')
                    ax2.set_ylabel('Probability')

                    ax3.set_xlabel('Range [km]')
                    ax3.set_ylabel('Azimuth Deviation [degrees]')

                    ax4.set_xlabel('Azimuth Deviation [degrees]')
                    ax4.set_ylabel('Probability')

                    plt.suptitle('Path Geometry Statistics (' + az_dirs[n_az] + ')', fontsize=18)
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
                    if az_dev[rng_az_mask].shape[0] > self.win_min_pts:
                        az_dev_mns[n_az][nr] = np.mean(az_dev[rng_az_mask], dtype=np.float64)
                        az_dev_std[n_az][nr] = max(np.std(az_dev[rng_az_mask], dtype=np.float64), self.min_az_dev_std)
                    else:
                        if nr == 0:
                            az_dev_mns[n_az][nr] = 0.0
                            az_dev_std[n_az][nr] = self.default_az_dev_std
                        else:
                            az_dev_mns[n_az][nr] = az_dev_mns[n_az][nr - 1] * 0.75
                            az_dev_std[n_az][nr] = self.default_az_dev_std + (az_dev_std[n_az][nr - 1] - self.default_az_dev_std) * 0.5

                    # set tropospheric contribution to reciprocal celerity fit
                    combo_mask = np.logical_and(rng_az_mask, rcel < self.tropo_strat_bnd + self.bnd_overlap)

                    wt1 = float(len(rngs[combo_mask]))
                    if rcel[combo_mask].shape[0] > self.win_min_pts:
                        rcel_mns[n_az][nr][0] = np.mean(rcel[combo_mask], dtype=np.float64)
                        rcel_std[n_az][nr][0] = max(np.std(rcel[combo_mask], dtype=np.float64) * 1.25, self.rcel_std_min)
                    else:
                        if nr == 0:
                            rcel_mns[n_az][nr][0] = self.mns0[0]
                            rcel_std[n_az][nr][0] = self.std0[0]
                        else:
                            rcel_mns[n_az][nr][0] = self.mns0[0] + (rcel_mns[n_az][nr - 1][0] - self.mns0[0]) * 0.25
                            rcel_std[n_az][nr][0] = self.std0[0] + (rcel_std[n_az][nr - 1][0] - self.std0[0]) * 0.25

                    # set stratospheric contribution to reciprocal celerity fit
                    strat_mask = np.logical_and(self.tropo_strat_bnd - self.bnd_overlap <= rcel, rcel <= self.strat_therm_bnd + self.bnd_overlap)
                    combo_mask = np.logical_and(rng_az_mask, strat_mask)

                    wt2 = float(len(rngs[combo_mask]))
                    if rcel[combo_mask].shape[0] > self.win_min_pts:
                        rcel_mns[n_az][nr][1] = np.mean(rcel[combo_mask], dtype=np.float64)
                        rcel_std[n_az][nr][1] = max(np.std(rcel[combo_mask], dtype=np.float64) * 1.25, self.rcel_std_min)
                    else:
                        if nr == 0:
                            rcel_mns[n_az][nr][1] = self.mns0[1]
                            rcel_std[n_az][nr][1] = self.std0[1]
                        else:
                            rcel_mns[n_az][nr][1] = self.mns0[1] + (rcel_mns[n_az][nr - 1][1] - self.mns0[1]) * 0.25
                            rcel_std[n_az][nr][1] = self.std0[1] + (rcel_std[n_az][nr - 1][1] - self.std0[1]) * 0.25

                    # set thermospheric contribution to reciprocal celerity fit
                    combo_mask = np.logical_and(rng_az_mask, self.strat_therm_bnd - self.bnd_overlap < rcel)

                    wt3 = float(len(rngs[combo_mask]))
                    if rcel[combo_mask].shape[0] > self.win_min_pts:
                        rcel_mns[n_az][nr][2] = np.mean(rcel[combo_mask], dtype=np.float64)
                        rcel_std[n_az][nr][2] = max(np.std(rcel[combo_mask], dtype=np.float64) * 1.25, self.rcel_std_min)
                    else:
                        if nr == 0:
                            rcel_mns[n_az][nr][2] = self.mns0[2]
                            rcel_std[n_az][nr][2] = self.std0[2]
                        else:
                            rcel_mns[n_az][nr][2] = self.mns0[2] + (rcel_mns[n_az][nr - 1][2] - self.mns0[2]) * 0.25
                            rcel_std[n_az][nr][2] = self.std0[2] + (rcel_std[n_az][nr - 1][2] - self.std0[2]) * 0.25

                    # set weights of reciprocal celerity distribution
                    if len(rngs[rng_az_mask]) > self.win_min_pts:
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

                for nr in range(rng_cnt):
                    rcel_wts[n_az][nr] *= rng_wts[nr] / np.sum(rng_wts)

                plt.close('all')

            # Normalize weights by total arrivals at all azimuths
            for n_az in range(self.az_bin_cnt):
                rcel_wts[n_az] *= az_wts[n_az] / np.sum(az_wts)

            priors = [0] * 6
            priors[0] = rng_bins
            priors[1] = az_dev_mns
            priors[2] = az_dev_std
            priors[3] = rcel_mns
            priors[4] = rcel_std
            priors[5] = rcel_wts

            pickle.dump(priors, open(output_file, "wb"))

    def load(self, model_file, smooth=None):
        fit_params = pickle.load(open(model_file, "rb"), encoding='latin1')

        self.rng_max = max(fit_params[0])

        self.az_dev_mns = [0] * self.az_bin_cnt
        self.az_dev_std = [0] * self.az_bin_cnt

        self.rcel_wts = [0] * self.az_bin_cnt
        self.rcel_mns = [0] * self.az_bin_cnt
        self.rcel_std = [0] * self.az_bin_cnt

        for n_az in range(self.az_bin_cnt):
            self.rcel_mns[n_az] = [0] * 3
            self.rcel_std[n_az] = [0] * 3
            self.rcel_wts[n_az] = [0] * 3

        if smooth:
            print("Loading propagation model parameters from " + model_file + " with smoothing.")
            for n_az in range(self.az_bin_cnt):
                self.az_dev_mns[n_az] = interp1d(fit_params[0], savgol_filter(fit_params[1][n_az], 9, 3), kind='cubic')
                self.az_dev_std[n_az] = interp1d(fit_params[0], savgol_filter(fit_params[2][n_az], 9, 3), kind='cubic')

                self.rcel_mns[n_az] = [0] * 3
                self.rcel_std[n_az] = [0] * 3
                self.rcel_wts[n_az] = [0] * 3

                for j in range(3):
                    self.rcel_mns[n_az][j] = interp1d(fit_params[0], savgol_filter(fit_params[3][n_az][:, j], 9, 3), kind='cubic')
                    self.rcel_std[n_az][j] = interp1d(fit_params[0], savgol_filter(fit_params[4][n_az][:, j], 9, 3), kind='cubic')
                    self.rcel_wts[n_az][j] = interp1d(fit_params[0], savgol_filter(fit_params[5][n_az][:, j], 9, 3), kind='cubic')
        else:
            print("Loading propagation model parameters from " + model_file + " without smoothing.")
            for n_az in range(self.az_bin_cnt):
                self.az_dev_mns[n_az] = interp1d(fit_params[0], fit_params[1][n_az], kind='cubic')
                self.az_dev_std[n_az] = interp1d(fit_params[0], fit_params[2][n_az], kind='cubic')

                self.rcel_mns[n_az] = [0] * 3
                self.rcel_std[n_az] = [0] * 3
                self.rcel_wts[n_az] = [0] * 3

                for j in range(3):
                    self.rcel_mns[n_az][j] = interp1d(fit_params[0], fit_params[3][n_az][:, j], kind='cubic')
                    self.rcel_std[n_az][j] = interp1d(fit_params[0], fit_params[4][n_az][:, j], kind='cubic')
                    self.rcel_wts[n_az][j] = interp1d(fit_params[0], fit_params[5][n_az][:, j], kind='cubic')

    def display(self, file_id=None, subtitle=None):
        resol = 100
        rngs = np.linspace(0.0, 1000.0, resol)
        bias = np.empty([resol])
        width = np.empty([resol])

        bias_color = 'Blue'
        var_color = 'LightBlue'

        compass_file = imp.find_module('stochprop')[1] + '/compass.png'

        f1, ax = plt.subplots(3, 3, figsize=(12, 9))

        for n1, n2 in itertools.product(list(range(3)), repeat=2):
            if n1 != 1 or n2 != 1:
                ax[n1, n2].set_ylim([-10, 10])
                ax[n1, n2].set_xlim([0, 1000])
                ax[n1, n2].set_xticks([0, 250, 500, 750, 1000])
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
            plt.suptitle("Azimuth Deviation Statistics" , fontsize=22)
        plt.show(block=False)

        def plot_az_stats(axis_id, azimuth):
            bias = self.eval_az_dev_mn(rngs, [azimuth] * len(rngs))
            width = self.eval_az_dev_vr(rngs, [azimuth] * len(rngs))
            axis_id.fill_between(rngs, bias - 2.0 * width, bias + 2.0 * width, facecolor=var_color)
            axis_id.plot(rngs, bias, linewidth=2.0, color=bias_color)
            axis_id.plot(rngs, [0] * len(rngs), 'k:')
            plt.pause(0.1)

        plot_az_stats(ax[0, 0],  -45.0)
        plot_az_stats(ax[0, 1],    0.0)
        plot_az_stats(ax[0, 2],   45.0)
        plot_az_stats(ax[1, 0],  -90.0)
        plot_az_stats(ax[1, 2],   90.0)
        plot_az_stats(ax[2, 0], -135.0)
        plot_az_stats(ax[2, 1],  180.0)
        plot_az_stats(ax[2, 2],  135.0)

        if file_id:
            plt.savefig(file_id + "_az-dev.png", bbox_inches='tight', dpi=200)

        # Plot celerity-range statistics
        cels = np.linspace(0.2, 0.4, resol)
        rngs = np.linspace(0, 1000.0, resol)
        R, V = np.meshgrid(rngs, cels)
        R = R.flatten()
        V = V.flatten()

        pdf = np.empty([resol, resol])

        palette = cm.nipy_spectral_r

        az_dirs = [-180.0, -135.0, -90.0, -45.0, 0.0, 45.0, 90.0, 135.0]
        pdf_max = 0.0
        for az in az_dirs:
            pdf_max = max(pdf_max, max(self.eval_rcel_gmm(R, 1.0 / V, [az] * len(R))))

        f2, ax = plt.subplots(3, 3, figsize=(12, 9))

        for n1, n2 in itertools.product(list(range(3)), repeat=2):
            if n1 != 1 or n2 != 1:
                ax[n1, n2].set_ylim([0.2, 0.4])
                ax[n1, n2].set_xlim([0, 1000])
                ax[n1, n2].set_xticks([0, 250, 500, 750, 1000])
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
            plt.suptitle("Celerity-Range Statistics" , fontsize=22)
        plt.show(block=False)

        pdf = self.eval_rcel_gmm(R, 1.0 / V, [-45.0] * len(R))
        rngcel_plot = ax[0, 0].scatter(R, V, c=pdf, cmap=palette, marker=".", alpha=1.0, edgecolor='none', vmin=0.0, vmax=pdf_max)
        f1.colorbar(rngcel_plot, ax=[ax[0,2], ax[1,2], ax[2,2]], label="Probability")
        plt.pause(0.1)
        
        def plot_celerity_stats(axis_id, azimuth):
            pdf = self.eval_rcel_gmm(R, 1.0 / V, [azimuth] * len(R))
            axis_id.scatter(R, V, c=pdf, cmap=palette, marker=".", alpha=1.0, edgecolor='none', vmin=0.0, vmax=pdf_max)
            plt.pause(0.1)

        plot_celerity_stats(ax[0, 1],    0.0)
        plot_celerity_stats(ax[0, 2],   45.0)
        plot_celerity_stats(ax[1, 0],  -90.0)
        plot_celerity_stats(ax[1, 2],   90.0)
        plot_celerity_stats(ax[2, 0], -135.0)
        plot_celerity_stats(ax[2, 1],  180.0)
        plot_celerity_stats(ax[2, 2],  135.0)

        if file_id:
            plt.savefig(file_id + "_cel-rng.png", bbox_inches='tight', dpi=200)


class TLossModel(object):
    az_bin_cnt = 8
    az_bin_wdth = 60.0
    
    def __init__(self):
        self.rng_vals = [0]
        self.tloss_vals = [0]
        self.pdf_vals = [0] * self.az_bin_cnt
        self.pdf_fits = [0] * self.az_bin_cnt
    

    def build(self, tloss_file, output_file, show_fits=False):
        if os.path.isfile(output_file):
            print(output_file + " already exists  --->  Skipping transmission loss model construction...")
        else:
            print('Builing transmission loss models from file:', tloss_file)

            # read in data, convert tloss to dB relative to 1 km, and wrap azimuths to [-180.0:180.0]
            print('\t' + "Reading in data...")
            rngs, az, tloss_re, tloss_im, tloss_coh = np.loadtxt(tloss_file, unpack=True)
        
            output_rngs = np.sort(np.unique(rngs)[::2])
        
            az[az > 180.0] -= 360.0
            az[az < -180.0] += 360.0
        
            tloss = 10.0 * np.log10(np.sqrt(tloss_re**2 + tloss_im**2) * 1000.0)
            # tloss = 10.0 * np.log10(tloss_coh * 1000.0)

            tloss[np.isneginf(tloss)] = min(tloss[np.isfinite(tloss)])
            tloss[np.isposinf(tloss)] = max(tloss[np.isfinite(tloss)])
        
            tloss_vals = np.linspace(-75.0, 0.0, len(output_rngs))
            pdf_vals = np.empty((self.az_bin_cnt, len(output_rngs), len(tloss_vals)))
        
            for az_index in range(self.az_bin_cnt):
                center = -180 + 360.0 / self.az_bin_cnt * az_index
                if az_index == 0:
                    az_mask = np.logical_or(az >= 180.0 - self.az_bin_wdth / 2.0, az <= -180.0 + self.az_bin_wdth / 2.0)
                else:
                    az_mask = np.logical_and(center - self.az_bin_wdth / 2.0 <= az, az <= center + self.az_bin_wdth / 2.0)
        
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
                
                    plt.suptitle("Stochastic Transmission Loss Model \n Azimuth: " + az_dirs[az_index], fontsize=18)
                    plt.show(block=False)
                
                    ax1.plot(rngs[az_mask][::11], tloss[az_mask][::11], 'ko', markersize=1)
                    plt.pause(0.001)
            
                print('\t' + "Propagation direction (" + az_dirs[az_index] + ")...")
            
                # Define tloss pdf at each range point from KDE
                for nr, rng_val in enumerate(output_rngs):
                    masked_tloss = tloss[np.logical_and(az_mask, rngs == rng_val)]
                
                    if np.std(masked_tloss) < 0.025:
                        pdf_vals[az_index][nr] = norm.pdf(tloss_vals, loc=np.mean(masked_tloss), scale=0.025)
                    else:
                        kernel = gaussian_kde(masked_tloss)
                        pdf_vals[az_index][nr] = kernel.evaluate(tloss_vals)
            
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
    
    
    def load(self, model_file):
        priors = pickle.load(open(model_file, "rb"), encoding='latin1')
        
        self.rng_vals = priors[0]
        self.tloss_vals = priors[1]
        
        for az_index in range(self.az_bin_cnt):
            self.pdf_vals[az_index] = priors[2][az_index]
            self.pdf_fits[az_index] = RectBivariateSpline(self.rng_vals, self.tloss_vals, self.pdf_vals[az_index])


    def eval(self, rng, tloss, az):
        az_index = find_azimuth_bin(az, self.az_bin_cnt)
    
        if len(np.atleast_1d(rng)) == 1:
            result = self.pdf_fits[az_index].ev(rng, tloss)
        else:
            result = np.empty(len(rng))
            for n_az in range(self.az_bin_cnt):
                mask = az_index==n_az
                if np.any(mask):
                    result[mask] = self.pdf_fits[n_az].ev(np.array(rng)[mask], np.array(tloss)[mask])
    
        return result


    def display(self, title="Transmission Loss Statistics", file_id=None):
        scale_max = 0.1
        compass_file = imp.find_module('stochprop')[1] + '/compass.png'
        
        resol = 100
        rngs = np.linspace(0.0, 1000.0, resol)
        tloss_min, tloss_max = -60.0, 0.0
        tloss = np.linspace(tloss_min, tloss_max, resol)
        
        R, TL = np.meshgrid(rngs, tloss)
        R = R.flatten()
        TL = TL.flatten()
        
        palette = cm.nipy_spectral_r
        f1, ax = plt.subplots(3, 3, figsize=(12, 9))
        
        for n1, n2 in itertools.product(range(3), repeat=2):
            if n1 != 1 or n2 != 1:
                ax[n1, n2].set_ylim([tloss_min, tloss_max])
                ax[n1, n2].set_xlim([0, 1000])
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
        f1.colorbar(tloss_plot, ax=[ax[0,2], ax[1,2], ax[2,2]], label="Probability")
        plt.pause(0.01)

        def plot_tloss_stats(axis_id, azimuth):
            pdf = self.eval(R, TL, np.array([azimuth] * len(R)))
            axis_id.scatter(R, TL, c=pdf, cmap=palette, marker='o', s=[12.5] * len(R), alpha=0.5, edgecolor='none', vmin=0.0, vmax=scale_max)
            plt.pause(0.01)

        plot_tloss_stats(ax[0, 1],    0.0)
        plot_tloss_stats(ax[0, 2],   45.0)
        plot_tloss_stats(ax[1, 0],  -90.0)
        plot_tloss_stats(ax[1, 2],   90.0)
        plot_tloss_stats(ax[2, 0], -135.0)
        plot_tloss_stats(ax[2, 1],  180.0)
        plot_tloss_stats(ax[2, 2],  135.0)

        if file_id:
            plt.savefig(file_id + "_tloss.png", bbox_inches='tight', dpi=200)

        plt.close('all')


