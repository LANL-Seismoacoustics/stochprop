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
import imp
import itertools
import subprocess

import numpy as np

from pyproj import Geod

from scipy.interpolate import interp1d, RectBivariateSpline
from scipy.stats import norm, gaussian_kde
from scipy.signal import savgol_filter

import matplotlib.cm as cm
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

sph_proj = Geod(ellps='sphere')

plt.rcParams.update({'font.size': 18})


# ############################## #
#  Running infraga and NCPAprop  #
#     with multiple profiles     #
# ############################## #
def run_infraga(profs_path, results_file, pattern="*.met", cpu_cnt=None, geom="3d", bounces=25, inclinations=[1.0, 60.0, 1.0], azimuths=[-180.0, 180.0, 3.0], freq=0.1, z_grnd=0.0, rng_max=1000.0, src_loc=[0.0, 0.0, 0.0], infraga_path="", clean_up=False):
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
        clean_up: boolean
            Flag to remove individual [..].arrival.dat files after combining
        """

    if os.path.isfile(results_file):
        print(results_file + " already exists  --->  Skipping infraGA/GeoAc ray tracing runs...")
    else:
        if not ((geom is "3d") or (geom is "sph")):
            msg = "Incompatible geometry option for infraga: {}.  Options are 3d' and 'sph'".format(geom)
            warnings.warn(msg)
        else:
            dir_files = np.sort(os.listdir(profs_path))
            for file_name in dir_files:
                if fnmatch.fnmatch(file_name, pattern):
                    file_id = os.path.splitext(file_name)[0]

                    if os.path.isfile(profs_path + "/" + file_id + ".arrivals.dat"):
                        print("InfraGA/GeoAc arrivals for " + file_name + " already finished...")
                    else:
                        print("Generating ray paths for " + file_name)
                        if cpu_cnt:
                            command = "mpirun --oversubscribe -np " + str(cpu_cnt) + " " + infraga_path + " infraga-accel-" + geom + " -prop "
                        else:
                            command = infraga_path + " infraga-" + geom + " -prop "

                        command = command + profs_path + "/" + file_name
                        command = command + " incl_min=" + str(inclinations[0]) + " incl_max=" + str(inclinations[1]) + " incl_step=" + str(inclinations[2])
                        command = command + " az_min=" + str(azimuths[0]) + " az_max=" + str(azimuths[1]) + " az_step=" + str(azimuths[2])
                        if geom == "sph":
                            command = command + " src_lat=" + str(src_loc[0]) + " src_lon=" + str(src_loc[1])
                        command = command + " src_alt=" + str(src_loc[2])
                        command = command + " freq=" + str(freq) + " z_grnd=" + str(z_grnd) + " max_rng=" + str(rng_max)
                        command = command + " calc_amp=False" + " bounces=" + str(bounces) + " write_rays=false" # + " > /dev/null &"

                        print(command)
                        subprocess.call(command, shell=True)

            command = "cat " + profs_path + "/*.arrivals.dat > " + results_file
            print(command)
            subprocess.call(command, shell=True)
            
            if clean_up:
                command = "rm "  + profs_path + "/*.dat"
                print(command)
                subprocess.call(command, shell=True)


def run_modess(profs_path, results_path, pattern="*.met", azimuths=[-180.0, 180.0, 3.0], freq=0.1, z_grnd=0.0, rng_max=1000.0, ncpaprop_path="", clean_up=False):
    """
        Run the NCPAprop normal mode methods to compute transmission
        loss values for a suite of atmospheric specifications at
        a set of frequency values

        Parameters
        ----------
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
        clean_up: boolean
            Flag to remove individual .nm files after combining
        """
# 
    dir_files = np.sort(os.listdir(profs_path))

    if os.path.isfile(results_path + ".lossless.nm"):
        print(results_path + ".lossless.nm already exists  --->  Skipping NCPAprop modess runs...")
    else:
        for file_name in dir_files:
            if fnmatch.fnmatch(file_name, pattern):
                file_id = os.path.splitext(file_name)[0]

                print('\t' + "Running NCPAprop modess for " + profs_path + "/" + file_name + " at " + "%.3f" % freq + " Hz...")
                command = ncpaprop_path + "Modess --multiprop --atmosfile " + profs_path + "/" + file_name + " --freq " + str(freq)
                command = command + " --azimuth_start " + str(azimuths[0]) + " --azimuth_end " + str(azimuths[1]) + " --azimuth_step " + str(azimuths[2])
                command = command + "  --maxrange_km " + str(rng_max) + " --zground_km " + str(z_grnd) # + " > /dev/null &"
                print(command)
                subprocess.call(command, shell=True)

                subprocess.call("mv Nby2D_tloss_1d.lossless.nm " + profs_path + "/" + file_id + "_%.3f" % freq + "Hz.lossless.nm", shell=True)
                subprocess.call("mv Nby2D_tloss_1d.nm " + profs_path + "/" + file_id + "_%.3f" % freq + "Hz.nm", shell=True)
        
        print('\t' + "Combining transmission loss predictions..." + '\n')
        subprocess.call("cat " + profs_path + "/*_%.3f" % freq + "Hz.lossless.nm > " + results_path + ".lossless.nm", shell=True)
        subprocess.call("cat " + profs_path + "/*_%.3f" % freq + "Hz.nm > " + results_path + ".nm", shell=True)

        if clean_up:
            subprocess.call("rm " + profs_path + "/*_%.3f" % freq + "Hz.lossless.nm", shell=True)
            subprocess.call("rm " + profs_path + "/*_%.3f" % freq + "Hz.nm", shell=True)


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

    def build(self, arrivals_file, output_file, show_fits=False, rng_width=50.0, rng_spacing=10.0, geom="3d", src_loc=[0.0, 0.0, 0.0], min_turning_ht=0.0, az_bin_cnt=16, az_bin_wdth=30.0):
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

        """
        if os.path.isfile(output_file):
            print(output_file + " already exists  --->  Skipping path geometry model construction...")
        else:
            if not ((geom is "3d") or (geom is "sph")):
                msg = "Incompatible geometry option for infraga: {}.  Options are 3d' and 'sph'".format(geom)
                warnings.warn(msg)
            else:
                print('Building celerity and azimuth statistics from file:', arrivals_file)

                self._az_bin_cnt = az_bin_cnt

                # define range bins and parameter arrays
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

                    for nr in range(rng_cnt):
                        rcel_wts[n_az][nr] *= rng_wts[nr] / np.sum(rng_wts)

                    plt.close('all')

                # Normalize weights by total arrivals at all azimuths
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
        print('\t' + "Maximum range: " + str(self._rng_max) + '\n')

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
            print('\t' + "Loading data into GMM with smoothing...")
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
            print('\t' + "Loading data into GMM without smoothing...")
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

    def display(self, file_id=None, subtitle=None, show_colorbar=True):
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
        rngs = np.linspace(0.0, 1000.0, resol)
        bias = np.empty([resol])
        width = np.empty([resol])

        bias_color = 'Blue'
        var_color = 'LightBlue'

        compass_file = imp.find_module('stochprop')[1] + '/resources/compass.png'

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
        rngs = np.linspace(0, 1000.0, resol)
        R, V = np.meshgrid(rngs, cels)
        R = R.flatten()
        V = V.flatten()

        pdf = np.empty([resol, resol])

        palette = cm.nipy_spectral_r

        pdf_max = 0.0
        for az in np.arange(-180.0, 180.0, 45.0):
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


class TLossModel(object):
    _az_bin_cnt = 16

    def __init__(self):
        self.rng_vals = [0]
        self.tloss_vals = [0]

    def build(self, tloss_file, output_file, show_fits=False, use_coh=False, az_bin_cnt=16, az_bin_wdth=30.0):
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
            rngs, az, tloss_re, tloss_im, tloss_coh = np.loadtxt(tloss_file, unpack=True)

            output_rngs = np.sort(np.unique(rngs)[::2])

            az[az > 180.0] -= 360.0
            az[az < -180.0] += 360.0

            if use_coh:
                tloss = 10.0 * np.log10(tloss_coh)
            else:
                tloss = 10.0 * np.log10(np.sqrt(tloss_re**2 + tloss_im**2))
            
            tloss[np.isneginf(tloss)] = min(tloss[np.isfinite(tloss)])
            tloss[np.isposinf(tloss)] = max(tloss[np.isfinite(tloss)])

            tloss_vals = np.linspace(-75.0, 0.0, len(output_rngs))
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

                print('\t' + "Propagation direction (" + str(center) + ")...")

                norm_mask = np.logical_and(az_mask, rngs == min(output_rngs))
                if use_coh:
                    tloss_norm = 10.0 * np.log10(np.mean(tloss_coh[norm_mask]))
                else:
                    tloss_norm = 10.0 * np.log10(np.mean(np.sqrt(tloss_re**2 + tloss_im**2)[norm_mask]))

                # Define tloss pdf at each range point from KDE
                for nr, rng_val in enumerate(output_rngs):
                    masked_tloss = tloss[np.logical_and(az_mask, rngs == rng_val)] - tloss_norm

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

    def display(self, file_id=None, title="Transmission Loss Statistics", show_colorbar=True):
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
        compass_file = imp.find_module('stochprop')[1] + '/resources/compass.png'

        resol = 100
        rngs = np.linspace(0.0, 1000.0, resol)
        tloss = np.linspace(-60.0, 0.0, resol)

        R, TL = np.meshgrid(rngs, tloss)
        R = R.flatten()
        TL = TL.flatten()

        palette = cm.nipy_spectral_r
        f1, ax = plt.subplots(3, 3, figsize=(12, 9))

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

        if show_colorbar:
            cbar = f1.colorbar(tloss_plot, ax=[ax[0, 2], ax[1, 2], ax[2, 2]], label="Probability", aspect=40)
            cbar.ax.set_yticklabels([]) 

        plt.pause(0.01)

        def plot_tloss_stats(axis_id, azimuth):
            pdf = self.eval(R, TL, np.array([azimuth] * len(R)))
            axis_id.scatter(R, TL, c=pdf, cmap=palette, marker='o', s=[12.5] * len(R), alpha=0.5, edgecolor='none', vmin=0.0, vmax=scale_max)
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

        plt.close('all')
