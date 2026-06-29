# cli.py
#
# Methods to compute the gravity wave spectra that perturb a G2S specification
# based on the publication by Drob et al. (2013).  Additional insights from the
# Lalande & Waxler (2016) analysis and Fritts & Alexander (2003) manuscripts on
# atmospheric gravity waves.  The source and saturation spectra are from the 
# Warner & McIntyre (1996) work.
#  
#
# Philip Blom (pblom@lanl.gov)

import os
import click

import gzip
import json

import numpy as np

import matplotlib.pyplot as plt 

from scipy.stats import gaussian_kde, norm
from scipy.optimize import curve_fit

from obspy import UTCDateTime

from infrapy.utils import data_io

from stochprop import utils as sp_utils


@click.command('eig-wvfrm2json', short_help="Create JSON detections from infraGA results")
@click.option("--wvfrm-file", help="InfraGA waveform output", prompt="infraGA waveform output: ")
@click.option("--output-label", help="Output label to write JSON file", prompt="Detection JSON file label: ")
@click.option("--origin-time", help="Origin time of synthetic event", prompt="Event origin time: ")
@click.option("--sta", help="Station code for synthetic detection(s)", prompt="Station code: ")
@click.option("--loc", help="Station code for synthetic detection(s)", default="")
@click.option("--cha", help="Station code for synthetic detection(s)", default="BDF")
@click.option("--c0", help="Sound speed to relate inclination to trace velocity", default=340.0)
@click.option("--ns-option", help="'min', 'mid', or 'max' IMS noise statistic", default="min")
@click.option("--array-dim", help="Array dimension (station size)", default=6)
def eig_wvfrm2json(wvfrm_file, output_label, origin_time, sta, loc, cha, c0, ns_option, array_dim):
    '''
    \b
    stochprop plot eof-fit
    -----------------------
    \b
    Example Usage:
    \t stochprop utils eig-wvfrm2json --wvfrm-file test.wvfrms.dat --output-label "test" --origin-time "2010-01-01T12:10:00" --sta "STA"
    
    '''

    click.echo("")
    click.echo("########################################")
    click.echo("##                                    ##")
    click.echo("##             stochprop              ##")
    click.echo("##           Utility Tools            ##")
    click.echo("##  infraGA eig-wvfrm 2 infraPy JSON  ##")
    click.echo("##                                    ##")
    click.echo("########################################")
    click.echo("")  

    click.echo('\n' + "Data I/O and Configuration Summary:")
    click.echo("  infraGA eigenray + waveform results file: " + wvfrm_file)
    click.echo("  Output JSON detection file: " + output_label + ".dets.json.gz")
    click.echo("  Event origin time: " + origin_time)
    click.echo("  Reference trace ID: SY." + sta + "." + loc + "." + cha)
    click.echo("  Ground sound speed [m/s]: " + str(c0))
    click.echo("  IMS noise option ('min', 'mid', or 'max'): " + ns_option)
    click.echo("  Array dimension: " + str(array_dim))
    click.echo("")


    dets_dict = sp_utils.eig_wvfrm2json(wvfrm_file, UTCDateTime(origin_time), c0, ns_option, array_dim, sta, loc, cha)
    click.echo("Writing " + str(len(dets_dict)) + " detections into JSON file." + '\n')

    det_output = {'det_info' : dets_dict}
    with gzip.open(output_label + ".dets.json.gz", 'wt', encoding='UTF-8') as zipfile:
        json.dump(det_output, zipfile, indent=4, cls=data_io.Infrapy_Encoder)



@click.command('celerity_gmm', short_help="Generate a GMM celerity model")
@click.option("--data-file", help="File containing celerity information", default=None, prompt="Arrivals data file: ")
@click.option("--output-label", help="File label for output RCG file", default=None, prompt="Output file label: ")
@click.option("--cel-index", help="Column index of celerity values", default=6)
@click.option("--atten-index", help="Column index of attenuation values", default=11)
@click.option("--atten-lim", help="Attenuation limit", default=None, type=float)
@click.option("--turn-ht-index", help="Column index of turning height values", default=7)
@click.option("--turn-ht-min", help="Turning height minimum", default=None, type=float)
@click.option("--broadening", help="Broadening factor", default=0.0, type=float)

def celerity_gmm(data_file, output_label, cel_index, atten_index, atten_lim, turn_ht_index, turn_ht_min, broadening):
    '''
    Compute a KDE of celerity values and generate parameters for a reciprocal celerity model

    \b
    Example usage (requires a data file with celerities):
    \tstochprop utils fit-celerity --data-file ToyAtmo.arrivals.dat

    '''

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##          stochprop          ##")
    click.echo("##        Utility Tools        ##")
    click.echo("##        celerity-gmm         ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")   

    celerity_gmm = {}
    celerity_gmm['arrivals file'] = data_file
    celerity_gmm['atten_lim'] = atten_lim
    celerity_gmm['turn_ht_min'] = turn_ht_min
    celerity_gmm['broadening'] = broadening

    click.echo("  Loading data from " + data_file)
    data = np.loadtxt(data_file)
    cel_data = data[:, cel_index]

    select_mask = np.ones_like(data[:, 0], dtype=bool)
    if atten_lim is not None:
        click.echo("  Applying mask to remove arrivals with greater than " + str(atten_lim) + " dB Sutherland & Bass attenuation prediction.")
        select_mask = np.logical_and(select_mask, data[:, atten_index] > atten_lim)

    if turn_ht_min is not None:
        click.echo("  Applying mask to remove arrivals with turning heights below " + str(turn_ht_min) + " km. ")
        select_mask = np.logical_and(select_mask, data[:, turn_ht_index] > turn_ht_min)

    click.echo("  Computing KDE of arrival celerities and fitting...")
    cel_kernel = gaussian_kde(1.0 / cel_data[select_mask])

    cel_vals = np.linspace(0.38, 0.18, 200)
    rcel_pdf = cel_kernel(1.0 / cel_vals)

    def rcel_func(rcel, wt1, wt2, wt3, mn1, mn2, mn3, std1, std2, std3):
        result = (wt1 / std1) * norm.pdf((rcel - mn1) / std1)
        result = result + (wt2 / std2) * norm.pdf((rcel - mn2) / std2)
        result = result + (wt3 / std3) * norm.pdf((rcel - mn3) / std3)

        return result
    
    gmm_p0 = [0.0539, 0.0899, 0.8562, 
                1.0 / 0.327, 1.0 / 0.293, 1.0 / 0.26,
                0.066, 0.08, 0.33]
    
    lower_bounds = [0.0, 0.0, 0.0,
                    1.0 / 0.36, 1.0 / 0.31, 1.0 / 0.29,
                    1.0e-5, 1.0e-5, 1.0e-5]
    
    upper_bounds = [1.0, 1.0, 1.0,
                    1.0 / 0.29, 1.0 / 0.26, 1.0 / 0.20,
                    1.0, 1.0, 1.0]

    gmm_bnds = (lower_bounds, upper_bounds)

    popt, _ = curve_fit(rcel_func, 1.0 / cel_vals, rcel_pdf,
                         p0=gmm_p0, bounds=gmm_bnds)
    
    popt[6:] = popt[6:] * (1.0 + broadening)
    popt = np.round(popt, 4)

    celerity_gmm['weights'] = list(popt[:3])
    celerity_gmm['means'] = list(popt[3:6])
    celerity_gmm['stdevs'] = list(popt[6:])

    print("  Writing Reciprocal Celerity GMM to " + output_label + ".rcg.json")
    with open(output_label + ".rcg.json", 'w') as json_file:
        json.dump(celerity_gmm, json_file, indent=4, cls=data_io.Infrapy_Encoder)

    # print info to screen and plot...
    click.echo('\n' + "celerity_gmm:")
    for key in celerity_gmm.keys():
        if celerity_gmm[key] is not None:
            click.echo("  " + key + ": " + str(celerity_gmm[key]))

    click.echo("    Note: mean reciprocal celerities: 1.0/" + str(np.round(1.0 / popt[3], 3)) + ", 1.0/" + str(np.round(1.0 / popt[4], 3)) + ", 1.0/" + str(np.round(1.0 / popt[5], 3)) + '\n')

    plt.figure(figsize=(7, 4))
    plt.plot(cel_vals, rcel_pdf, '-k', linewidth=4.0, label="Data KDE")
    plt.plot(cel_vals, rcel_func(1.0 / cel_vals, popt[0], popt[1], popt[2],popt[3], popt[4], popt[5], 
                                 popt[6], popt[7], popt[8]), '--r', linewidth=2.0, label="GMM Fit")
    plt.xlabel("Celerity [km/s]")
    plt.ylabel("Probability")
    plt.legend()
    plt.show()
