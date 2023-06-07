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
def ess_ratio(atmo_dir, results_path, atmo_pattern, atmo_format, year_selection, include_ns, show_title):
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

    plt.rcParams.update({'font.size': 18})

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
        print("Summary of 'stochprop prop season-trends' analysis", file=output_file)
        print("  Source directory: " + str(atmo_dir), file=output_file)
        print("  Specification pattern: " + str(atmo_pattern), file=output_file)
        print("  Specification format: " + str(atmo_format), file=output_file)
        if years_list is not None:
            print("  Limited years: " + str(years_list), file=output_file)

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

@click.command('eofs', short_help="Plot EOF results")
@click.option("--eofs-path", help="EOF output path and prefix (required)", prompt="EOF path: ")
@click.option("--eof-cnt", help="Number of EOFs to visualize (default: 5)", default=5)
def eofs(eofs_path, eof_cnt):
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

    sp_eofs._plot_eofs(eofs_path, eof_cnt)


@click.command('eof-fit', short_help="Plot EOF fit to a reference atmosphere")
@click.option("--atmo-file", help="Reference atmospheric specification (required)", prompt="Reference atmospheric specification: ")
@click.option("--eofs-path", help="EOF output path and prefix (required)", prompt="EOF path: ")
@click.option("--eof-cnt", help="Number of EOFs to visualize (default: 5)", default=5)
@click.option("--output-file", help="Output file to save fit (optional)", default=None)
def eof_fit(atmo_file, eofs_path, eof_cnt, output_file):
    '''
    \b
    stochprop plot eofs
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

    # NEED TO WRITE THIS VISUALIZATION FUNCTION

    sp_eofs.fit_atmo(atmo_file, eofs_path, output_path=output_file, eof_cnt=eof_cnt)

    return 0


@click.command('atmo-ensemble', short_help="Plot an ensemble of atmospheric states")
@click.option("--atmo-dir", help="Directory of atmospheric specifications (required)", prompt="Atmospheric specifications: ")
@click.option("--atmo-pattern", help="Specification file pattern (default: '*.dat')", default='*.dat')
@click.option("--atmo-format", help="Specification format (default: 'zTuvdp')", default='zTuvdp')
@click.option("--month-selection", help="Limit analysis to specific month(s) (default: None)", default=None)
@click.option("--week-selection", help="Limit analysis to specific week(s) (default: None)", default=None)
@click.option("--year-selection", help="Limit analysis to specific year(s) (default: None)", default=None)
@click.option("--max-alt", help="Maximum altitude for trimming data (default: None)", default=None)
@click.option("--plot-cnt", help="Number of ensemble samples to plot", default=25)
def atmo_ensemble(atmo_dir, atmo_pattern, atmo_format, month_selection, week_selection, year_selection, max_alt, plot_cnt):
    '''
    \b
    stochprop plot eofs
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
    ax[0].plot(c_mean, z0, '-k', linewidth=1.5)

    ax[1].set_xlabel("u [m/s]")
    for j in range(plot_cnt):
        ax[1].plot(u_vals[j], z0, color="xkcd:lightblue", linewidth=0.5)
    ax[1].plot(u_mean, z0, color="xkcd:blue", linewidth=1.5)

    ax[2].set_xlabel("v [m/s]")
    for j in range(plot_cnt):
        ax[2].plot(v_vals[j], z0, color="xkcd:red", linewidth=0.5)
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
def prop_model(model_file):
    '''
    \b
    stochprop plot prop-model
    ---------------------
    \b
    Example Usage:
    \t stochprop plot prop-model --model-file prop/winter/winter.pgm
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
        model.display(hold_fig=True)
    elif model_file[-3:] == 'tlm':
        model = sp_prop.TLossModel()
        model.load(model_file)
        model.display(hold_fig=True)
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
    \t stochprop plot detection-stats --tlm-files prop/US_RM/US_RM-winter_0.500Hz.tlm --yield-vals '1, 10, 100' --array-dim 5

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
@click.option("--src-yld", help="Explosive yield (kg eq. TNT)", default=10e3, type=float)
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
    \t stochprop prop network-performance --network-info network_test.dat --lat-min 36 --lat-max 42 --lon-min -117.5 --lon-max -107.5

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

    sp_prop.plot_network_performance(network_info, freq, src_yld, min_det_cnt, lat_min, lat_max, lon_min, lon_max, resol, figure_out, show_figure)

