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
import fnmatch

import numpy as np

from . import eofs
from . import gravity_waves as grav

def parse_option_list(input):
    if input is None:
        return input
    else:
        if "," in input:
            # remove white space and split by commas
            return input[1:-1].replace(" ", "").split(",")
        else:
            return input


@click.command('build', short_help="Build EOF info through SVD")
@click.option("--atmo-dir", help="Directory of atmospheric specifications (required)", prompt="Atmospheric specifications: ")
@click.option("--eofs-path", help="EOF output path and prefix (required)", prompt="Output path: ")
@click.option("--atmo-pattern", help="Specification file pattern (default: '*.dat')", default='*.dat')
@click.option("--atmo-format", help="Specification format (default: 'zTuvdp')", default='zTuvdp')
@click.option("--month-selection", help="Limit analysis to specific month(s) (default: None)", default=None)
@click.option("--week-selection", help="Limit analysis to specific week(s) (default: None)", default=None)
@click.option("--year-selection", help="Limit analysis to specific year(s) (default: None)", default=None)
@click.option("--save-datetime", help="Save date time info (default: False)", default=False)
@click.option("--max-alt", help="Maximum altitude for trimming data (default: None)", default=None)
@click.option("--eof-cnt", help="Number of EOFs to store (default: 100)", default=100)
def eof_build(atmo_dir, eofs_path, atmo_pattern, atmo_format, month_selection, week_selection, year_selection, save_datetime, max_alt, eof_cnt):
    '''
    \b
    stochprop eof build
    -----------------------
    \b
    Example Usage:
    \t stochprop eof build --atmo-dir profs/ --eofs-path eofs/example
    \t stochprop eof build --atmo-dir profs/ --eofs-path eofs/example_low_alt --max-alt 80.0 --eof-cnt 50
    \t stochprop eof build --atmo-dir profs/ --eofs-path eofs/example_winter --month-selection '[10, 11, 12, 01, 02, 03]'
    
    '''

    click.echo("")
    click.echo("##################################")
    click.echo("##                              ##")
    click.echo("##           stochprop          ##")
    click.echo("##          EOF Methods         ##")
    click.echo("##   Build SVD to Define EOFs   ##")
    click.echo("##                              ##")
    click.echo("##################################")
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

    click.echo("  EOF count: " + str(eof_cnt))
    click.echo("  Output path: " + str(eofs_path))
    if save_datetime:
        click.echo("  Saving date time info")
    click.echo("")
    
    A, z0, datetimes = eofs.build_atmo_matrix(atmo_dir, atmo_pattern, prof_format=atmo_format, months=months_list, weeks=weeks_list, years=years_list, return_datetime=True, max_alt=max_alt)
    if save_datetime:
        np.save(eofs_path + ".datetimes", datetimes)
    eofs.compute_eofs(A, z0, eofs_path, eof_cnt=eof_cnt)
    

@click.command('coeffs', short_help="Compute EOF coefficients")
@click.option("--atmo-dir", help="Directory of atmospheric specifications (required)", prompt="Atmospheric specifications: ")
@click.option("--eofs-path", help="EOF info path and prefix (required)", prompt="EOFs path: ")
@click.option("--coeff-path", help="Coefficients output path and prefix (required)", prompt="Output path: ")
@click.option("--atmo-pattern", help="Specification file pattern (default: '*.dat')", default='*.dat')
@click.option("--atmo-format", help="Specification format (default: 'zTuvdp'", default='zTuvdp')
@click.option("--month-selection", help="Limit analysis to specific month(s) (default=None)", default=None)
@click.option("--week-selection", help="Limit analysis to specific week(s) (default=None)", default=None)
@click.option("--year-selection", help="Limit analysis to specific year(s) (default=None)", default=None)
@click.option("--save-datetime", help="Save date time info (default: False)", default=False)
@click.option("--eof-cnt", help="Number of EOFs to use (default: 100)", default=100)
def eof_coeffs(atmo_dir, eofs_path, coeff_path, atmo_pattern, atmo_format, month_selection, week_selection, year_selection, save_datetime, eof_cnt):
    '''
    \b
    stochprop eof coeffs
    --------------------
    \b
    Example Usage (run in stochprop/examples/ after running 'stochprop eof build' example):
    \t stochprop eof coeffs --atmo-dir profs/ --eofs-path eofs/example_low_alt --coeff-path coeffs/example_low_alt_01 --month-selection '01' --eof-cnt 50
    \t stochprop eof coeffs --atmo-dir profs/ --eofs-path eofs/example_low_alt --coeff-path coeffs/example_low_alt_02 --month-selection '02' --eof-cnt 50

    '''

    click.echo("")
    click.echo("##################################")
    click.echo("##                              ##")
    click.echo("##           stochprop          ##")
    click.echo("##          EOF Methods         ##")
    click.echo("##    Coefficient Calculation   ##")
    click.echo("##                              ##")
    click.echo("##################################")
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
        click.echo("  Limited years:: " + str(years_list))

    click.echo("  EOFs path: " + str(eofs_path))
    click.echo("  EOF count: " + str(eof_cnt))
    click.echo("  Coefficient output path: " + str(coeff_path))

    if save_datetime:
        click.echo("  Saving date time info")

    click.echo("")

    A, z0, datetimes = eofs.build_atmo_matrix(atmo_dir, atmo_pattern, prof_format=atmo_format, months=months_list, weeks=weeks_list, years=years_list, return_datetime=True)
    if save_datetime:
        np.save(coeff_path + ".datetimes", datetimes)

    eofs.compute_coeffs(A, z0, eofs_path, coeff_path, eof_cnt=eof_cnt)


@click.command('seasonality', short_help="Compute EOF coefficient overlap and seasonal relations")
@click.option("--coeff-path", help="Coefficients output path and prefix (required)")
@click.option("--eofs-path", help="EOF info path and prefix (required)", prompt="EOFs path: ")
@click.option("--eof-cnt", help="Number of EOFs to use (default: 100)", default=100)
def eof_seasonality(coeff_path, eofs_path, eof_cnt):
    '''
    \b
    stochprop eof seasonality
    -------------------------
    \b
    Example Usage:
    \t stochprop eof seasonality --eofs-path eofs/example --coeff-path coeffs/example --eof-cnt 50

    '''

    click.echo("")
    click.echo("##################################")
    click.echo("##                              ##")
    click.echo("##           stochprop          ##")
    click.echo("##          EOF Methods         ##")
    click.echo("##     Seasonality Analysis     ##")
    click.echo("##                              ##")
    click.echo("##################################")
    click.echo("")  

    click.echo('\n' + "Run summary:")

    click.echo("  Coefficient path: " + str(coeff_path))
    click.echo("  EOFs path: " + str(eofs_path))
    click.echo("  EOF count: " + str(eof_cnt) + '\n')

    if os.path.isfile(coeff_path + "-overlap.npy"):
        click.echo("Existing overlap file identified...")
    else:
        coeff_vals = []
        if "/" in coeff_path:
            dir_files = os.listdir(os.path.dirname(coeff_path))
            for file in np.sort(dir_files):
                if fnmatch.fnmatch(file, os.path.basename(coeff_path) + '*-coeffs.npy'):
                    print("Loading " + os.path.dirname(coeff_path) + '/' + file)
                    coeff_vals = coeff_vals + [np.load(os.path.dirname(coeff_path) + '/' + file)]
        else:
            dir_files = os.listdir(".")
            for file in np.sort(dir_files):
                if fnmatch.fnmatch(file, coeff_path + '*-coeffs.npy'):
                    print("Loading " + file)
                    coeff_vals = coeff_vals + [np.load(file)]

        click.echo("")     
        overlap = eofs.compute_overlap(coeff_vals, eofs_path, eof_cnt=eof_cnt, method="kde") 
        np.save(coeff_path + "-overlap", overlap)

    click.echo("")
    eofs.compute_seasonality(coeff_path + "-overlap.npy", coeff_path)


@click.command('sample', short_help="Sample EOF coefficient KDEs to generate atmosphere realizations")
@click.option("--coeff-path", help="Coefficients output path and prefix (required)", prompt="coefficient paths:")
@click.option("--eofs-path", help="EOF info path and prefix (required)", prompt="EOFs path: ")
@click.option("--sample-path", help="Path for samples to be saved (required)", prompt="Sample (output) path:")
@click.option("--sample-cnt", help="Number of samples to generate (default: 100)", default=100)
@click.option("--eof-cnt", help="Number of EOFs to use (default: 100)", default=100)
def eof_sample(coeff_path, eofs_path, sample_path, sample_cnt, eof_cnt):
    '''
    \b
    stochprop eof sample 
    ---------------------
    \b
    Example Usage:
    \t stochprop eof sample --eofs-path eofs/example_winter --coeff-path coeffs/example_winter --sample-path samples/winter/example_winter --sample-cnt 50 --eof-cnt 50

    '''

    click.echo("")
    click.echo("##################################")
    click.echo("##                              ##")
    click.echo("##           stochprop          ##")
    click.echo("##          EOF Methods         ##")
    click.echo("##       Generate Samples       ##")
    click.echo("##                              ##")
    click.echo("##################################")
    click.echo("")  

    click.echo('\n' + "Run summary:")

    click.echo("  Coefficient path: " + str(coeff_path))
    click.echo("  EOFs path: " + str(eofs_path))
    click.echo("  Sample output path: " + str(sample_path))
    click.echo("  Sample count: " + str(sample_cnt))
    click.echo("  EOF count: " + str(eof_cnt) + '\n')

    if "-coeffs.npy" in coeff_path:
        coeff_vals = np.load(coeff_path)
    else:
        coeff_vals = np.load(coeff_path + "-coeffs.npy")

    eofs.sample_atmo(coeff_vals, eofs_path, sample_path, eof_cnt=eof_cnt, prof_cnt=sample_cnt)


@click.command('eof', short_help="Use EOFs to perturb a specification")
@click.option("--atmo-file", help="Reference atmospheric specification (required)", prompt="Reference atmospheric specification: ")
@click.option("--eofs-path", help="Path to EOF info (required)", prompt="EOF results: ")
@click.option("--out", help="Output prefix (required)", prompt="Output prefix: ")
@click.option("--std-dev", help="Standard deviation (default: 10 m/s)", default=10.0)
@click.option("--eof-max", help="Maximum EOF coefficient to use (default: 100)", default=100)
@click.option("--eof-cnt", help="Number of EOFs to use (default: 50)", default=50)
@click.option("--sample-cnt", help="Number of perturbed samples (default: 25)", default=25)
@click.option("--alt-weight", help="Altitude weighting power (default: 2.0)", default=2.0)
@click.option("--singular-value-weight", help="Sing. value weighting power (default: 0.25)", default=0.25)
def eof_perturb(atmo_file, eofs_path, out, std_dev, eof_max, eof_cnt, sample_cnt, alt_weight, singular_value_weight):
    '''
    \b
    stochprop perturb eof
    ---------------------
    \b
    Example Usage:
    \t stochprop perturb eof --atmo-file profs/g2stxt_2010010118_39.7393_-104.9900.dat --eofs-path eofs/example --out test

    '''

    click.echo("")
    click.echo("###################################")
    click.echo("##                               ##")
    click.echo("##           stochprop           ##")
    click.echo("##      Perturbation Methods     ##")
    click.echo("##       EOF Perturbations       ##")
    click.echo("##                               ##")
    click.echo("###################################")
    click.echo("") 

    click.echo('\n' + "Run summary:")
    click.echo("  Reference specification: " + str(atmo_file))
    click.echo("  EOFs path: " + str(eofs_path))
    click.echo("  EOF Max Index: " + str(eof_max))
    click.echo("  EOF Count: " + str(eof_cnt))
    click.echo("  Standard Deviation [m/s]: " + str(std_dev))
    click.echo("  Altitude Weighting: " + str(alt_weight))
    click.echo("  Singular Value weighting: " + str(singular_value_weight))
    click.echo("")

    click.echo("  Output Path: " + out)
    click.echo("  Sample Count: " + str(sample_cnt))
    click.echo("")

    eofs.perturb_atmo(atmo_file, eofs_path, out, stdev=std_dev, eof_max=eof_max, eof_cnt=eof_cnt, sample_cnt=sample_cnt, alt_wt_pow=alt_weight, sing_val_wt_pow=singular_value_weight)


@click.command('gravity-waves', short_help="Use gravity waves to perturb a specification")
@click.option("--atmo-file", help="Reference atmospheric specification (required)", prompt="Reference atmospheric specification: ")
@click.option("--out", help="Output prefix (required)", prompt="Output prefix: ")
@click.option("--sample-cnt", help="Number of perturbed samples (default: 25)", default=25)
@click.option("--t0", help="Propagation time from source [hr] (default: 8)", default=8.0)
@click.option("--dx", help="Horizontal wavenumber scale [km] (default: 4.0)", default=4.0)
@click.option("--dz", help="Altitude resolution [km] (default: 0.2)", default=0.2)
@click.option("--nk", help="Horizontal wavenumber resolution (default: 128)", default=128)
@click.option("--nom", help="Frequency resolution (default: 5)", default=5)
@click.option("--random-phase", help="Randomize phase at source [bool] (default: False)", default=False)
@click.option("--z-src", help="Source altitude [km] (default: 20.0)", default=20.0)
@click.option("--m-star", help="Source spectrum peak [1/km] (default: (2 pi) / 2.5)", default=(2.0 * np.pi) / 2.5)
@click.option("--cpu-cnt", help="Number of CPUs in parallel analysis (default: None)", default=None, type=int)
@click.option("--debug-fig", help="Output for figures to aid in debugging (default: None)", default=None, type=str)
def gravity_waves(atmo_file, out, sample_cnt, t0, dx, dz, nk, nom, random_phase, z_src, m_star, cpu_cnt, debug_fig):
    '''
    \b
    Gravity wave perturbation calculation based on Drob et al. (2013) method.

    \b
    Example Usage:
    \t stochprop perturb gravity-waves --atmo-file profs/g2stxt_2010010118_39.7393_-104.9900.dat --out test_gw

    '''

    click.echo("")
    click.echo("####################################")
    click.echo("##                                ##")
    click.echo("##           stochprop            ##")
    click.echo("##      Perturbation Methods      ##")
    click.echo("##   Gravity Wave Perturbations   ##")
    click.echo("##                                ##")
    click.echo("####################################")
    click.echo("") 

    click.echo('\n' + "Run summary:")
    click.echo("  Reference specification: " + str(atmo_file))
    click.echo("  t0: " + str(t0))
    click.echo("  dx: " + str(dx))
    click.echo("  dz: " + str(dz))
    click.echo("  Nk: " + str(nk))
    click.echo("  Nom: " + str(nom))
    click.echo("  random_phase: " + str(random_phase))
    click.echo("  z_src: " + str(z_src))
    click.echo("  m_star: " + str(m_star))
    if cpu_cnt is not None:
        click.echo("  cpu_cnt: " + str(cpu_cnt))
    click.echo("")

    click.echo("  Output Path: " + out)
    click.echo("  Sample Count: " + str(sample_cnt))
    click.echo("")

    grav.perturb_atmo(atmo_file, out, sample_cnt=sample_cnt, t0=t0 * 3600.0, dx=dx, dz=dz, Nk=nk, N_om=nom, random_phase=random_phase, z_src=z_src, m_star=m_star, env_below=False, cpu_cnt=cpu_cnt, fig_out=debug_fig)


if __name__ == '__main__':
    main()

