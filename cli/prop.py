# stochprop/cli/prop.py
#
# Command line interface for building propagation statistics
#
# Philip Blom (pblom@lanl.gov)

import click

import numpy as np

import matplotlib.pyplot as plt 

from scipy.stats import gaussian_kde, norm
from scipy.optimize import curve_fit

from stochprop import propagation


@click.command('fit-celerity', short_help="Parameterize a GMM celerity model", hidden=True)
@click.option("--data-file", help="File containing celerity information", default=None)
@click.option("--cel-index", help="Column index of celerity values", default=6)
@click.option("--atten-index", help="Column index of attenuation values", default=11)
@click.option("--atten-lim", help="Attenuation limit", default=None, type=float)
def fit_celerity(data_file, cel_index, atten_index, atten_lim):
    '''
    Compute a KDE of celerity values and generate parameters for a reciprocal celerity model

    \b
    Example usage (requires a data file with celerities):
        stochprop utils fit-celerity --data-file prop/winter/winter.arrivals.dat

    '''

    click.echo("")
    click.echo("#################################")
    click.echo("##                             ##")
    click.echo("##          stochprop          ##")
    click.echo("##     Propagation Methods     ##")
    click.echo("##  Parameterize Celerity GMM  ##")
    click.echo("##                             ##")
    click.echo("#################################")
    click.echo("")   


    click.echo("  Loading data from " + data_file)
    data = np.loadtxt(data_file)
    cel_data = data[:, cel_index]

    if atten_lim is not None:
        click.echo("  Building KDE with limited arrivals (" + str(atten_lim) + " dB Sutherland & Bass attenuation limit)")
        atten_data = data[:, atten_index]
        cel_kernel = gaussian_kde(1.0 / cel_data[atten_data > atten_lim])
    else:
        click.echo("  Building KDE for all arrival celerities")
        cel_kernel = gaussian_kde(1.0 / cel_data)

    cel_vals = np.linspace(0.38, 0.18, 200)
    rcel_pdf = cel_kernel(1.0 / cel_vals)

    click.echo("  Generating fit to KDE...")
    def rcel_func(rcel, wt1, wt2, wt3, mn1, mn2, mn3, std1, std2, std3):
        result = (wt1 / std1) * norm.pdf((rcel - mn1) / std1)
        result = result + (wt2 / std2) * norm.pdf((rcel - mn2) / std2)
        result = result + (wt3 / std3) * norm.pdf((rcel - mn3) / std3)

        return result
    
    popt, _ = curve_fit(rcel_func, 1.0 / cel_vals, rcel_pdf,
                         p0=[0.0539, 0.0899, 0.8562, 
                             1.0 / 0.327, 1.0 / 0.293, 1.0 / 0.26,
                             0.066, 0.08, 0.33])
    popt = np.round(popt, 3)

    click.echo('\n' + "  Reciprocal celerity model parameters (CLI and config file formats):")
    click.echo("    --rcel-wts '" + str(popt[0]) + ", " + str(popt[1]) + ", " + str(popt[2]) + "' --rcel-mns '" + str(popt[3]) + ", " + str(popt[4]) + ", " + str(popt[5]) + "' --rcel-sds '" + str(popt[6]) + ", " + str(popt[7]) + ", " + str(popt[8]) + "'" + '\n')

    click.echo("    rcel_wts = '" + str(popt[0]) + ", " + str(popt[1]) + ", " + str(popt[2]) + "'")
    click.echo("    rcel_mns = '" + str(popt[3]) + ", " + str(popt[4]) + ", " + str(popt[5]) + "'")
    click.echo("    rcel_sds = '" + str(popt[6]) + ", " + str(popt[7]) + ", " + str(popt[8]) + "'" + '\n')

    click.echo("    Note: mean reciprocal celerities: 1.0/" + str(np.round(1.0 / popt[3], 3)) + ", 1.0/" + str(np.round(1.0 / popt[4], 3)) + ", 1.0/" + str(np.round(1.0 / popt[5], 3)) + '\n')

    plt.figure(figsize=(7, 4))
    plt.plot(cel_vals, rcel_pdf, '-k', linewidth=4.0, label="Data KDE")
    plt.plot(cel_vals, rcel_func(1.0 / cel_vals, popt[0], popt[1], popt[2],popt[3], popt[4], popt[5], 
                                 popt[6], popt[7], popt[8]), '--r', linewidth=2.0, label="GMM Fit")
    plt.xlabel("Celerity [km/s]")
    plt.ylabel("Probability")
    plt.legend()
    plt.show()




@click.command('build-pgm', short_help="Build a path geometry model (PGM)")
@click.option("--atmos-dir", help="Directory containing atmospheric specifications", prompt="Path to directory with atmospheric specifications")
@click.option("--atmos-pattern", help="Atmosphere file pattern (default: '*.met')", default="*.met")
@click.option("--output-path", help="Path and prefix for PGM output", prompt="Path and prefix for PGM output")
@click.option("--src-loc", help="Source location (lat, lon, alt)", prompt="Enter source location (lat, lon, alt)")
@click.option("--inclinations", help="Inclination min, max, and step (default: [2, 50, 2]", default = "[2.0, 50.0, 2.0]")
@click.option("--azimuths", help="Azimuth min, max, and step (default: [0, 360, 6]", default = "[0.0, 360.0, 6.0]")
@click.option("--bounces", help="Number of ground bounces for ray paths (default: 25)", default=25)
@click.option("--z-grnd", help="Ground elevation for simulations (default: 0.0 km)", default=0.0)
@click.option("--rng-max", help="Maximum range for simulations (default: 1000.0 km)", default=1000.0)
@click.option("--freq", help="Frequency for Sutherland-Bass atten. (default: 0.5 Hz)", default=0.5)
@click.option("--prof-format", help="Profile format (default: 'zTuvdp')", default='zTuvdp')
@click.option("--infraga-path", help="Path to infraGA/Geoac binaries", default='')
@click.option("--clean-up", help="Remove individual results after merge (default: True)", default=True)
@click.option("--cpu-cnt", help="Number of CPUs for propagation simulations", default=None)
@click.option("--rng-window", help="Range window in PGM (default: 50 km)", default=50.0)
@click.option("--rng-step", help="Range resolution in PGM (default: 10 km)", default=10.0)
@click.option("--az-bin-cnt", help="Number of azimuth bins in PGM (default: 16)", default=16)
@click.option("--az-bin-width", help="Azimuth bin width in PGM (default: 30 deg)", default=30.0)
@click.option("--verbose", help="Output analysis stages as they're done.",  default=False)
def build_pgm(atmos_dir, atmos_pattern, output_path, src_loc, inclinations, azimuths, bounces, z_grnd, rng_max,
                freq, prof_format, infraga_path, clean_up, cpu_cnt, rng_window, rng_step, az_bin_cnt, az_bin_width, verbose):
    '''
    \b
    stochprop prop build-pgm 
    ---------------------
    \b
    Example Usage:
    \t stochprop prop build-pgm --atmos-dir samples/winter/ --output-path prop/winter/winter --src-loc '[30.0, -120.0, 0.0]'  --cpu-cnt 8

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##            stochprop            ##")
    click.echo("##       Propagation Methods       ##")
    click.echo("##       Path Geometry Model       ##")
    click.echo("##                                 ##")
    click.echo("######################################")
    click.echo("")  

    click.echo('\n' + "Data IO summary:")
    click.echo("  Atmospheric specifications directory: " + atmos_dir)
    click.echo("  Specification pattern: " + atmos_pattern)
    click.echo("  Model output path: " + output_path)

    click.echo('\n' + "infraGA/GeoAc parameters:")
    click.echo("  Source location: " + src_loc)
    click.echo("  Inclination angles (min, max, step): " + inclinations)
    click.echo("  Azimuth angles (min, max, step): " + azimuths)
    click.echo("  Bounces: " + str(bounces))
    click.echo("  Ground elevation: " + str(z_grnd))
    click.echo("  Range max: " + str(rng_max))
    click.echo("  Frequency: " + str(freq))
    click.echo("  Clean up: " + str(clean_up))
    if cpu_cnt is not None:
        click.echo("  CPU count: " + str(cpu_cnt))

    click.echo('\n' + "Path Geometry Model (PGM) parameters:")
    click.echo("  Range window: " + str(rng_window))
    click.echo("  Range step: " + str(rng_step))
    click.echo("  Azimuth bin count: " + str(az_bin_cnt))
    click.echo("  Azimuth bin width: " + str(az_bin_width))
    click.echo("")

    src_loc = [float(val) for val in src_loc.strip(' ()[]').split(',')]
    inclinations = [float(val) for val in inclinations.strip(' ()[]').split(',')]
    azimuths = [float(val) for val in azimuths.strip(' ()[]').split(',')]
    
    propagation.run_infraga(atmos_dir, output_path + ".arrivals.dat", pattern=atmos_pattern, cpu_cnt=cpu_cnt, geom="sph", bounces=bounces, 
                    inclinations=inclinations, azimuths=azimuths, freq=freq, z_grnd=z_grnd, rng_max=rng_max, src_loc=src_loc, infraga_path=infraga_path, clean_up=clean_up, prof_format=prof_format, verbose=verbose)

    pgm = propagation.PathGeometryModel()
    pgm.build(output_path + ".arrivals.dat", output_path + ".pgm", geom="sph", src_loc=src_loc, rng_width=rng_window, 
                rng_spacing=rng_step, az_bin_cnt=az_bin_cnt, az_bin_wdth=az_bin_width)



@click.command('build-tlm', short_help="Build a transmission loss model (TLM)")
@click.option("--atmos-dir", help="Directory containing atmospheric specifications", prompt="Path to directory with atmospheric specifications")
@click.option("--output-path", help="Path and prefix for TLM output", prompt="Path and prefix for TLM output")
@click.option("--atmos-pattern", help="Atmosphere file pattern (default: '*.met')", default="*.met")
@click.option("--topo-label", help="Path and label for terrain files (optional)", default=None)
@click.option("--ncpaprop-method", help="NCPAprop method ('modess' or 'epape')", default='modess')
@click.option("--ncpaprop-path", help="Path to NCPAprop binaries (if not on path)", default="")
@click.option("--freq", help="Frequency for simulation (default: 0.5 Hz)", default=0.5)
@click.option("--azimuths", help="Azimuth min, max, and step (default: [0, 360, 3]", default = "[0.0, 360.0, 3.0]")
@click.option("--z-grnd", help="Ground elevation for simulations (default: 0.0)", default=0.0)
@click.option("--rng-max", help="Maximum range for simulations (default: 1000.0)", default=1000.0)
@click.option("--rng-resol", help="Range resolution for output (default: 1.0)", default=1.0)
@click.option("--clean-up", help="Clean up files after merge (default: True)", default=True)
@click.option("--verbose", help="Show NCPAprop output (default: False)", default=False)
@click.option("--cpu-cnt", help="Number of CPUs for propagation simulations", default=None)
@click.option("--az-bin-cnt", help="Number of azimuth bins in TLM (default: 16)", default=16)
@click.option("--az-bin-width", help="Azimuth bin width in TLM (default: 30 deg)", default=30.0)
@click.option("--rng-lims", help="Range limits in TLM (default: [1, 1000])", default='[1, 1000.0]')
@click.option("--rng-cnt", help="Range intervals in TLM (default: 100)", default=100)
@click.option("--rng-spacing", help="Option for range sampling ('linear' or 'log')", default='linear')
@click.option("--use-coherent-tl", help="Use coherent transmission loss (default: False", default=False)
def build_tlm(atmos_dir, output_path, atmos_pattern, topo_label, ncpaprop_method, ncpaprop_path, freq, azimuths, z_grnd,
              rng_max, rng_resol, clean_up, verbose, cpu_cnt, az_bin_cnt, az_bin_width, rng_lims, rng_cnt, rng_spacing, 
              use_coherent_tl):
    '''
    \b
    stochprop prop build-tlm
    ------------------------
    \b
    Example Usage:
    \t stochprop prop build-tlm --atmos-dir samples/winter/ --output-path prop/winter/winter --freq 0.2  --cpu-cnt 8

    '''

    click.echo("")
    click.echo("#######################################")
    click.echo("##                                   ##")
    click.echo("##             stochprop             ##")
    click.echo("##        Propagation Methods        ##")
    click.echo("##      Transmission Loss Model      ##")
    click.echo("##                                   ##")
    click.echo("########################################")
    click.echo("")  

    click.echo('\n' + "Data IO summary:")
    click.echo("  Atmospheric specifications directory: " + atmos_dir)
    click.echo("  Specification pattern: " + atmos_pattern)
    if topo_label is not None:
        click.echo("Topography files label: " + topo_label)
    click.echo("  Model output path: " + output_path)

    click.echo('\n' + "NCPAprop parameters:")
    if len(ncpaprop_path) > 0:
        click.echo("  NCPAprop path:", ncpaprop_path)
    click.echo("NCPAprop method: " + ncpaprop_method)
    click.echo("  Frequency: " + str(freq))
    click.echo("  Azimuth angles (min, max, step): " + azimuths)
    click.echo("  Ground elevation: " + str(z_grnd))
    click.echo("  Range max, resolution: " + str(rng_max) + ", " + str(rng_resol))
    click.echo("  Clean up: " + str(clean_up))
    if cpu_cnt is not None:
        click.echo("  CPU count: " + str(cpu_cnt))

    click.echo('\n' + "Transmission Loss Model (TLM) parameters:")
    click.echo("  Azimuth bin count: " + str(az_bin_cnt))
    click.echo("  Azimuth bin width: " + str(az_bin_width))
    click.echo("  Range limits: " + str(rng_lims))
    click.echo("  Range count: " + str(rng_cnt))
    click.echo("  Range spacing: " + str(rng_spacing))
    click.echo("  Use coherent transmission loss: " + str(use_coherent_tl))
    click.echo("")

    azimuths = [float(val) for val in azimuths.strip(' ()[]').split(',')]
    rng_lims = [float(val) for val in rng_lims.strip(' ()[]').split(',')]

    if cpu_cnt is None:
        cpu_cnt = 1
    else:
        cpu_cnt = int(cpu_cnt)

    propagation.run_ncpaprop(ncpaprop_method, atmos_dir, output_path + "_%.3f" %freq + "Hz", pattern=atmos_pattern, 
                             azimuths=azimuths, freq=freq, z_grnd=z_grnd, rng_max=rng_max, rng_resol=rng_resol, 
                             ncpaprop_path=ncpaprop_path, topo_path_label=topo_label, clean_up=clean_up, cpu_cnt=cpu_cnt,
                             verbose=verbose)

    tlm = propagation.TLossModel()
    tlm.build(output_path + "_%.3f" %freq + "Hz.nm", output_path + "_%.3f" %freq + "Hz.tlm", use_coh=use_coherent_tl, az_bin_cnt=az_bin_cnt,
                az_bin_wdth=az_bin_width, rng_lims=rng_lims, rng_cnt=rng_cnt, rng_smpls=rng_spacing)
