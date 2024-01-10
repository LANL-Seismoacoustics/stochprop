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
import re
import click
import fnmatch

import numpy as np

import matplotlib.pyplot as plt 
import matplotlib.cm as cm 
from matplotlib.colorbar import ColorbarBase

from datetime import datetime

from stochprop import propagation as sp_prop
from stochprop import eofs as sp_eofs

def parse_option_list(input):
    if input is None:
        return input
    else:

        if "," in input:
            # remove white space and brackets/parentheses, then split by commas
            temp_list = []
            for item in input.replace(" ", "").strip('([])').split(','):
                if ':' in item:
                    str_format = "%0" + str(len(item.strip('([])').split(":")[0])) + "d"
                    temp_list = temp_list + [str_format % val for val in range(int(item.strip('([])').split(":")[0]),
                                                                                int(item.strip('([])').split(":")[1]) + 1)]
                else:
                    temp_list = temp_list + [item]
            return temp_list
        elif ':' in input:
            str_format = "%0" + str(len(input.strip('([])').split(":")[0])) + "d"
            return [str_format % val for val in range(int(input.strip('([])').split(":")[0]),
                                                        int(input.strip('([])').split(":")[1]) + 1)]
        else:
            return input


@click.command('ess-ratio', short_help="Compute seasonality via the effective sound speed ratio")
@click.option("--atmo-dir", help="Directory of atmospheric specifications (required)", prompt="Atmospheric specifications: ")
@click.option("--results-path", help="Output path and prefix", default=None)
@click.option("--atmo-pattern", help="Specification file pattern (default: '*.dat')", default='*.dat')
@click.option("--atmo-format", help="Specification format (default: 'zTuvdp')", default='zTuvdp')
@click.option("--year-selection", help="Limit analysis to specific year(s) (default: None)", default=None)
@click.option("--include-NS", help="Option to include north/south analysis", default=False)
@click.option("--show-title", help="Option to display title text", default=True)
@click.option("--xkcd-mode", help="Use XKCD plotting style", default=False, hidden=True)
def ess_ratio(atmo_dir, results_path, atmo_pattern, atmo_format, year_selection, include_ns, show_title, xkcd_mode):
    '''
    \b
    stochprop plot ess-ratio
    -----------------------
    \b
    Example Usage:
    \t stochprop plot ess-ratio --atmo-dir profs/ --results-path example
    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##            stochprop            ##")
    click.echo("##      Visualization Methods      ##")
    click.echo("##   ESS Ratio Seasonal Analysis   ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("")  

    years_list = parse_option_list(year_selection)
    
    click.echo('\n' + "Run summary:")
    click.echo("  Source directory: " + str(atmo_dir))
    click.echo("  Specification pattern: " + str(atmo_pattern))
    click.echo("  Specification format: " + str(atmo_format))
    if results_path is not None:
        click.echo("  Output path: " + str(results_path))
    if years_list is not None:
        click.echo("  Limited years: " + str(years_list))
    if include_ns:
        click.echo("  Include NS: True")
    click.echo("")
    
    A, z0, datetimes = sp_eofs.build_atmo_matrix(atmo_dir, pattern=atmo_pattern, prof_format=atmo_format, years=years_list, return_datetime=True)
    for line in open(atmo_dir + [file for file in os.listdir(atmo_dir) if fnmatch.fnmatch(file, atmo_pattern)][0], 'r'):
        if "Ground Height" in line:
            grnd_ht = float(line[18:])
            break
    grnd_index = np.argmin(abs(z0 - grnd_ht))
    click.echo('\t\t' + "Extracted ground elevation: " + str(grnd_ht))

    if results_path is not None:
        out_header = "Summary of 'stochprop prop season-trends' analysis" + '\n'
        out_header = out_header + "  Source directory: " + str(atmo_dir) + '\n'
        out_header = out_header + "  Specification pattern: " + str(atmo_pattern) + '\n'
        out_header = out_header + "  Specification format: " + str(atmo_format) + '\n'
        if years_list is not None:
            out_header = out_header + "  Limited years: " + str(years_list)

    if xkcd_mode:
        with plt.xkcd():
            sp_prop._plot_ess_ratio(A, z0, datetimes, grnd_index, results_path, include_ns, show_title, out_header)
    else:
        sp_prop._plot_ess_ratio(A, z0, datetimes, grnd_index, results_path, include_ns, show_title, out_header)


@click.command('eofs', short_help="Plot EOF results")
@click.option("--eofs-path", help="EOF output path and prefix (required)", prompt="EOF path: ")
@click.option("--eof-cnt", help="Number of EOFs to visualize (default: 5)", default=5)
@click.option("--xkcd-mode", help="Use XKCD plotting style", default=False, hidden=True)
def eofs(eofs_path, eof_cnt, xkcd_mode):
    '''
    \b
    stochprop plot eofs
    -----------------------
    \b
    Example Usage:
    \t stochprop plot eofs --eofs-path eofs/example
    
    '''
    click.echo("")
    click.echo("#######################################")
    click.echo("##                                   ##")
    click.echo("##             stochprop             ##")
    click.echo("##       Visualization Methods       ##")
    click.echo("##        EOF Analysis Results       ##")
    click.echo("##                                   ##")
    click.echo("#######################################")
    click.echo("")  

    if xkcd_mode:
        with plt.xkcd():
            sp_eofs._plot_eofs(eofs_path, eof_cnt)
    else:
        sp_eofs._plot_eofs(eofs_path, eof_cnt)


@click.command('eof-fit', short_help="Plot EOF fit to a reference atmosphere")
@click.option("--atmo-file", help="Reference atmospheric specification (required)", prompt="Reference atmospheric specification: ")
@click.option("--eofs-path", help="EOF output path and prefix (required)", prompt="EOF path: ")
@click.option("--eof-cnt", help="Number of EOFs to visualize (default: 5)", default=5)
@click.option("--output-file", help="Output file to save fit (optional)", default=None)
def eof_fit(atmo_file, eofs_path, eof_cnt, output_file):
    '''
    \b
    stochprop plot eof-fit
    -----------------------
    \b
    Example Usage:
    \t stochprop plot eof-fit --atmo-file profs/g2stxt_2011010118_39.1026_-84.5123.dat --eofs-path eofs/example_winter --eof-cnt 25
    
    '''

    click.echo("")
    click.echo("#######################################")
    click.echo("##                                   ##")
    click.echo("##             stochprop             ##")
    click.echo("##       Visualization Methods       ##")
    click.echo("##       EOF Fit to Atmosphere       ##")
    click.echo("##                                   ##")
    click.echo("#######################################")
    click.echo("")  

    sp_eofs.fit_atmo(atmo_file, eofs_path, output_path=output_file, eof_cnt=eof_cnt)



@click.command('atmo-ensemble', short_help="Plot an ensemble of atmospheric states")
@click.option("--atmo-dir", help="Directory of atmospheric specifications (required)", prompt="Atmospheric specifications: ")
@click.option("--atmo-pattern", help="Specification file pattern (default: '*.dat')", default='*.dat')
@click.option("--atmo-format", help="Specification format (default: 'zTuvdp')", default='zTuvdp')
@click.option("--month-selection", help="Limit analysis to specific month(s) (default: None)", default=None)
@click.option("--week-selection", help="Limit analysis to specific week(s) (default: None)", default=None)
@click.option("--year-selection", help="Limit analysis to specific year(s) (default: None)", default=None)
@click.option("--max-alt", help="Maximum altitude for trimming data (default: None)", default=None)
@click.option("--plot-cnt", help="Number of ensemble samples to plot (default: 25)", default=25)
@click.option("--ref-atmo", help="Reference atmosphere to plot", default=None)
def atmo_ensemble(atmo_dir, atmo_pattern, atmo_format, month_selection, week_selection, year_selection, max_alt, plot_cnt, ref_atmo):
    '''
    \b
    stochprop plot atmo-ensemble
    -----------------------
    \b
    Example Usage:
    \t stochprop plot atmo-ensemble --atmo-dir samples/winter/
    
    '''

    click.echo("")
    click.echo("#######################################")
    click.echo("##                                   ##")
    click.echo("##             stochprop             ##")
    click.echo("##       Visualization Methods       ##")
    click.echo("##        Atmosphere Ensemble        ##")
    click.echo("##                                   ##")
    click.echo("#######################################")
    click.echo("") 


    months_list = parse_option_list(month_selection)
    weeks_list = parse_option_list(week_selection)
    years_list = parse_option_list(year_selection)

    click.echo('\n' + "Run summary:")
    click.echo("  Source directory: " + str(atmo_dir))
    click.echo("  Specification pattern: " + str(atmo_pattern))
    click.echo("  Specification format: " + str(atmo_format))
    if months_list is not None:
        click.echo("  Limited months: " + str(months_list))
    if weeks_list is not None:
        click.echo("  Limited weeks: " + str(weeks_list))
    if years_list is not None:
        click.echo("  Limited years: " + str(years_list))
    if max_alt is not None:
        click.echo("  max_alt: " + str(max_alt))
        max_alt = float(max_alt)
    if ref_atmo is not None:
        click.echo("  Reference atmosphere:" + ref_atmo)

    click.echo("")
    
    A, z0 = sp_eofs.build_atmo_matrix(atmo_dir, atmo_pattern, prof_format=atmo_format, months=months_list, weeks=weeks_list, years=years_list, max_alt=max_alt)

    # compute the differences to fit with EOFs
    file_len = int(A.shape[1] / 5)
    u_vals = A[:, file_len * 1:file_len * 2]
    v_vals = A[:, file_len * 2:file_len * 3]

    d_vals = A[:, file_len * 3:file_len * 4]
    p_vals = A[:, file_len * 4:file_len * 5]
    c_vals = np.sqrt((sp_eofs.gam / 10.0) * p_vals / d_vals)

    c_mean = np.mean(c_vals, axis=0)
    u_mean = np.mean(u_vals, axis=0)
    v_mean = np.mean(v_vals, axis=0)

    # Plot the atmospheric data
    plt.rcParams.update({'font.size': 12})
    fig, ax = plt.subplots(1, 3, figsize=(8, 6), sharey=True)
    ax[0].set_ylabel("Altitude [km]")
    ax[0].set_ylim(min(z0), max(z0))
    
    ax[0].set_xlabel("c [m/s]")
    for j in range(plot_cnt):
        ax[0].plot(c_vals[j], z0, color='0.5', linewidth=0.5)

    ax[1].set_xlabel("u [m/s]")
    for j in range(plot_cnt):
        ax[1].plot(u_vals[j], z0, color="xkcd:lightblue", linewidth=0.5)

    ax[2].set_xlabel("v [m/s]")
    for j in range(plot_cnt):
        ax[2].plot(v_vals[j], z0, color="xkcd:red", linewidth=0.5)

    if ref_atmo is not None:
        ref_atmo_data = np.loadtxt(ref_atmo)
        ax[0].plot(np.sqrt((sp_eofs.gam / 10.0) * ref_atmo_data[:, 5] / ref_atmo_data[:, 4]), ref_atmo_data[:, 0], '-k', linewidth=1.5)
        ax[1].plot(ref_atmo_data[:, 2], ref_atmo_data[:, 0], color="xkcd:blue", linewidth=1.5)
        ax[2].plot(ref_atmo_data[:, 3], ref_atmo_data[:, 0], color="darkred", linewidth=1.5)
    else:
        ax[0].plot(c_mean, z0, '-k', linewidth=1.5)
        ax[1].plot(u_mean, z0, color="xkcd:blue", linewidth=1.5)
        ax[2].plot(v_mean, z0, color="darkred", linewidth=1.5)

    plt.tight_layout()
    plt.show()

    return 0


@click.command('coeff-overlap', short_help="Plot EOF coefficient overlap for seasonality")
@click.option("--overlap", help="Coefficients output path and prefix (required)")
def coeff_overlap(overlap):
    '''
    \b
    stochprop plot coeff-overlap
    -----------------------
    \b
    Example Usage:
    \t stochprop plot coeff-overlap --overlap coeff/example-overlap.npy
    
    '''

    click.echo("")
    click.echo("#######################################")
    click.echo("##                                   ##")
    click.echo("##             stochprop             ##")
    click.echo("##       Visualization Methods       ##")
    click.echo("##  Coefficient Overlap Seasonality  ##")
    click.echo("##                                   ##")
    click.echo("#######################################")
    click.echo("") 


    click.echo("")
    sp_eofs.compute_seasonality(overlap, overlap[:-12])


@click.command('prop-model', short_help="Visualize a PGM or TLM")
@click.option("--model-file", help="stochprop PGM or TLM file (required)", prompt="PGM or TLM file")
@click.option("--output-id", help="File ID for output image file", default=None)
@click.option("--cmap-max", help="Value maximum for colormap (default='auto')", default=None)
def prop_model(model_file, output_id, cmap_max):
    '''
    \b
    stochprop plot prop-model
    ---------------------
    \b
    Example Usage:
    \t stochprop plot prop-model --model-file prop/winter/winter.pgm --output-id winter
    \t stochprop plot prop-model --model-file prop/winter/winter_0.200Hz.tlm

    '''
    click.echo("")
    click.echo("#######################################")
    click.echo("##                                   ##")
    click.echo("##             stochprop             ##")
    click.echo("##       Visualization Methods       ##")
    click.echo("##         Propagation Model         ##")
    click.echo("##                                   ##")
    click.echo("#######################################")
    click.echo("")  

    if model_file[-3:] == 'pgm':
        model = sp_prop.PathGeometryModel()
        model.load(model_file)
        model.display(file_id=output_id, hold_fig=True, cmap_max=cmap_max)
    elif model_file[-3:] == 'tlm':
        model = sp_prop.TLossModel()
        model.load(model_file)
        model.display(file_id=output_id, hold_fig=True)
    else:
        click.echo("Error: invalid model file.")


@click.command('detection-stats', short_help="Visualize single station detection statistics")
@click.option("--tlm-files", help="Path and name(s) for TLM file(s)", prompt="Enter tlm file info: ")
@click.option("--yield-vals", help="List of yield values (tons eq. TNT)", prompt="Enter yield value(s): ")
@click.option("--array-dim", help="Array dimension (number of sensors)", default=1, type=int)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--show-figure", help="Print figure to screen", default=True)
def detection_stats(tlm_files, yield_vals, array_dim, figure_out, show_figure):
    '''
    \b
    stochprop prop detection_stats 
    ---------------------
    \b
    Example Usage:
    \t stochprop plot detection-stats --tlm-files 'prop/US_RM/US_RM-winter_*.tlm' --yield-vals '1, 10, 100' --array-dim 5

    '''
    click.echo("")
    click.echo("#######################################")
    click.echo("##                                   ##")
    click.echo("##             stochprop             ##")
    click.echo("##        Propagation Methods        ##")
    click.echo("##       Detection Statistics        ##")
    click.echo("##                                   ##")
    click.echo("#######################################")
    click.echo("")  

    # Load TLMs
    if "*" in tlm_files:
        tlm_dir = os.path.dirname(tlm_files)   
        tlm_pattern = tlm_pattern = tlm_files.split("/")[-1]
        tlm_file_list = [tlm_dir + "/" + file_name for file_name in np.sort(os.listdir(tlm_dir)) if fnmatch.fnmatch(file_name, tlm_pattern + "*")]
    elif "," in tlm_files:
        tlm_file_list = tlm_files.replace(" ","").split(",")
    else:
        tlm_file_list = [tlm_files]

    models = [0] * 2
    models[0] = [float(re.compile(r'\d+\.\d+').findall(file_name)[-1]) for file_name in tlm_file_list]
    models[1] = [0] * len(tlm_file_list)
    for n in range(len(tlm_file_list)):
        models[1][n] = sp_prop.TLossModel()
        models[1][n].load(tlm_file_list[n])

    # Generate plot
    yield_vals = [float(val)*1.0e3 for val in yield_vals.strip(' ()[]').split(',')]
    sp_prop.plot_detection_stats(models, yield_vals, array_dim=array_dim, show_fig=show_figure, output_path=figure_out)



@click.command('network-performance', short_help="Visualize detection statistics for a network")
@click.option("--network-info", help="Text file of network info")
@click.option("--freq", help="Frequency for analysis", default=None)
@click.option("--src-yld", help="Explosive yield (tons eq. TNT) (default 10)", default=10, type=float)
@click.option("--min-det-cnt", help="Minimum detecting stations", default=3, type=int)
@click.option("--resol", help="Grid resolution", default=100, type=int)
@click.option("--lat-min", help="Minimum latitude of grid", default=30.0)
@click.option("--lat-max", help="Maximum latitude of grid", default=40.0)
@click.option("--lon-min", help="Minimum longitude of grid", default=-110.0)
@click.option("--lon-max", help="Maximum longitude of grid", default=-105.0)
@click.option("--figure-out", help="Destination for figure", default=None)
@click.option("--show-figure", help="Print figure to screen", default=True)
def network_performance(network_info, freq, src_yld, min_det_cnt, resol, lat_min, lat_max, lon_min, lon_max, figure_out, show_figure):
    '''
    \b
    stochprop prop detection_stats 
    ---------------------
    \b
    Example Usage:
    \t stochprop plot network-performance --network-info network_test.dat --lat-min 36 --lat-max 42 --lon-min -117.5 --lon-max -107.5

    '''
    click.echo("")
    click.echo("#######################################")
    click.echo("##                                   ##")
    click.echo("##             stochprop             ##")
    click.echo("##        Propagation Methods        ##")
    click.echo("##        Network Performance        ##")
    click.echo("##                                   ##")
    click.echo("########################################")
    click.echo("")  

    sp_prop.plot_network_performance(network_info, freq, src_yld * 1.0e3, min_det_cnt, lat_min, lat_max, lon_min, lon_max, resol, figure_out, show_figure)

