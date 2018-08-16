# propagation.py
#
# Run infraga-accel-3d or NCPAprop's modess to generate predictions to
# use for constructing stochastic propagation models
#
# Philip Blom (pblom@lanl.gov)

import sys
import os

import numpy as np

from infrapy.utils import prog_bar

ncpa_prop_dir = "/Users/pblom/Research/Coding/Packages/ncpaprop-1.3.2/ncpaprop-1.3.2/"

def calc_rays(profs_path, prof_rng, results_file, freq=0.1, z_grnd=0.0, inclinations=[0.5, 45.0, 0.5], azimuths=[-180.0, 180.0, 3.0]):
    open(results_file, 'w').close()
    for n in prof_rng:
        filename = profs_path + "%02d" % n + ".met"
        print("Generating ray paths for " + filename)
        os.system("mpirun -np 7 infraga-accel-3d -prop " + filename + " bounces=0 incl_min=" + str(inclinations[0]) + " incl_max=" + str(inclinations[1]) + " incl_step=" + str(inclinations[2]) + " az_min=" + str(azimuths[0]) + " az_max=" + str(azimuths[1]) + " az_step=" + str(azimuths[2]) + " freq=" + str(freq) + " z_grnd=" + str(z_grnd) + " calc_amp=False")

        os.system("cat " + filename[:-4] + ".results.dat >> " +  results_file)
        os.system("rm " + filename[:-4] + "*.dat")


'''
run_nm requires making modifications to the
NCPAprop source code as follows:

In the file "SolveModNB.cpp", at line 1490, replace

  FILE *tloss_1d    = fopen("tloss_1d.nm", "w");
  FILE *tloss_ll_1d = fopen("tloss_1d.lossless.nm","w");

with

  char output_buffer [60];
  sprintf(output_buffer, "tloss_1d-%.3fHz.nm", freq);
  FILE *tloss_1d    = fopen(output_buffer,"w");

  sprintf(output_buffer, "tloss_1d-%.3fHz.lossless.nm", freq);
  FILE *tloss_ll_1d = fopen(output_buffer,"w");
'''
def run_nm(profile, azimuth, freq, id, prog_step=0, z_grnd = 0.0):
    os.system(ncpa_prop_dir + "bin/Modess --atmosfile " + profile + " --atmosfileorder ztuvdp --skiplines 0 --azimuth " + str(azimuth) + " --freq " + str(freq) + " --zground_km " + str(z_grnd) + " > /dev/null")
    os.system("mv tloss_1d-%.3f.nm result." % freq + str(id) + ".dat")
    os.system("mv tloss_1d-%.3f.lossless.nm result.lossless." % freq + str(id) + ".dat")
    prog_bar.increment(prog_step)


def run_nm_wrapper(args): return run_nm(*args)


def calc_tloss(profs_path, prof_rng, results_id, freq_vals, pool, azimuths=[-180.0, 180.0, 3.0], z_grnd=0.0, crit_angle_lim = 2.0, prog_scalar = 3):
    for n in prof_rng:
        filename = profs_path + "%02d" % n + ".met"
        print("Generating ray paths for " + filename)

        # Read in profile to check for tropospheric/stratospheric duct
        profile = np.loadtxt(filename)

        alt_vals = profile[:, 0]
        snd_spd_vals = np.sqrt(1.4 * 287.0 * profile[:, 1])
        u_vals = profile[:, 2]
        v_vals = profile[:, 3]


        # Define the tloss output files at each frequency
        tloss1_out = [''] * len(freq_vals)
        tloss2_out = [''] * len(freq_vals)
        for nf, freq in enumerate(freq_vals):
            tloss1_out[nf] = open(results_id + "-%02d" % n + "_%.3f" % freq + "Hz.dat", 'w')
            tloss2_out[nf] = open(results_id + "-%02d" % n + "_%.3f" % freq + "Hz.lossless.dat", 'w')

        # Cycle through azimuths computing the transmission loss and writing it to file
        for azimuth in np.arange(azimuths[0], azimuths[1], azimuths[2]):
            print('\t' + "Azimuth: " + str(azimuth) + " degrees..." + '\t', end=' ')

            c_eff = snd_spd_vals + u_vals * np.sin(np.radians(azimuth)) + v_vals * np.cos(np.radians(azimuth))

            c_eff_grnd = c_eff[np.argmin(abs(alt_vals - z_grnd))]
            c_eff_max = max(c_eff[np.logical_and(z_grnd + 1.0 < alt_vals, alt_vals < 70.0)])

            print("c_eff ratio (max vs. grnd): ", c_eff_max / c_eff_grnd, '\t', end= ' ')

            # Calculate transmission loss if the critical angle is at least the critical angle limit
            if c_eff_grnd / c_eff_max < np.cos(np.radians(crit_angle_lim)):
                print("Critical angle: ", np.degrees(np.arccos(c_eff_grnd / c_eff_max)), '\t', end = ' ')

                args_list = [0] * len(freq_vals)
                for nf, freq in enumerate(freq_vals):
                    args_list[nf] = (filename, azimuth, freq, nf, prog_scalar, z_grnd)

                prog_bar.prep(len(freq_vals) * prog_scalar)
                pool.map(run_nm_wrapper, args_list)
                prog_bar.close()

                # Open each of the results files, convert the format from [r, re(A), im(A)] to [r, az, |A|], and write to tloss_out
                # using scaling to reference point at 1 km (first entry of the file)
                for nf in range(len(freq_vals)):
                    results = np.loadtxt("result." + str(nf) + ".dat")
                    ref_amp = np.sqrt(results[0][1]**2 + results[0][2]**2)
                    for line in results:
                        print(line[0], azimuth, np.sqrt(line[1]**2 + line[2]**2) / ref_amp, file=tloss1_out[nf])
                    print('', file=tloss1_out[nf])

                    results = np.loadtxt("result.lossless." + str(nf) + ".dat")
                    ref_amp = np.sqrt(results[0][1]**2 + results[0][2]**2)
                    for line in results:
                        print(line[0], azimuth, np.sqrt(line[1]**2 + line[2]**2) / ref_amp, file=tloss2_out[nf])
                    print(' ', file=tloss2_out[nf])
                # Remove result files and .nm files
                os.system("rm result.*.dat *.nm")
            else:
                print("No low- or mid-altitude waveguide or critical angle is too shallow (< " + str(crit_angle_lim) + " degrees).")

        # Close completed tloss files
        for nf in range(len(freq_vals)):
            tloss1_out[nf].close()
            tloss2_out[nf].close()

