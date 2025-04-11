# stochprop/cli/prop.py
#
# Command line interface for building propagation statistics
#
# Philip Blom (pblom@lanl.gov)

import click

import numpy as np
import configparser as cnfg

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
@click.option("--atmo-dir", help="Directory containing atmospheric specifications", prompt="Path to directory with atmospheric specifications")
@click.option("--atmo-pattern", help="Atmosphere file pattern (default: '*.met')", default="*.met")
@click.option("--output-path", help="Path and prefix for PGM output", prompt="Path and prefix for PGM output")
@click.option("--src-loc", help="Source location (lat, lon, alt)", prompt="Enter source location (lat, lon, alt)")
@click.option("--inclinations", help="Inclination min, max, step (default: [2, 50, 2]", default = "[2.0, 50.0, 2.0]")
@click.option("--azimuths", help="Azimuth min, max, and step (default: [0, 360, 6]", default = "[0.0, 360.0, 6.0]")
@click.option("--bounces", help="Number of ground bounces for paths (default: 25)", default=25)
@click.option("--z-grnd", help="Ground elevation (default: 0.0 km)", default=0.0)
@click.option("--rng-max", help="Maximum range (default: 1000.0 km)", default=1000.0)
@click.option("--freq", help="Freq. for non-geometric atten. (default: 0.5 Hz)", default=0.5)
@click.option("--prof-format", help="Profile format (default: 'zTuvdp')", default='zTuvdp')
@click.option("--infraga-path", help="Path to infraGA/Geoac binaries", default='')
@click.option("--local-temp-dir", help="Local storage for individual infraGA results", default=None)
@click.option("--cpu-cnt", help="Number of CPUs for propagation simulations", default=None)
@click.option("--rng-window", help="Range window in PGM (default: 50 km)", default=50.0)
@click.option("--rng-step", help="Range resolution in PGM (default: 10 km)", default=10.0)
@click.option("--az-bin-cnt", help="Number of azimuth bins in PGM (default: 16)", default=16)
@click.option("--az-bin-width", help="Azimuth bin width in PGM (default: 30 deg)", default=30.0)
@click.option("--min-turning-ht", help="Minimum turning height altitude [km]", default=0.0)
@click.option("--station-centered", help="Build a station-centered model", default=False)
@click.option("--topo-file", help="Terrain file for propagation simulation", default=None)
@click.option("--verbose", help="Output analysis stages as they're done.",  default=False)
def build_pgm(atmo_dir, atmo_pattern, output_path, src_loc, inclinations, azimuths, bounces, z_grnd, rng_max,
                freq, prof_format, infraga_path, local_temp_dir, cpu_cnt, rng_window, rng_step, az_bin_cnt, az_bin_width,
                min_turning_ht, station_centered, topo_file, verbose):
    '''
    \b
    stochprop prop build-pgm 
    ---------------------
    \b
    Example Usage:
    \t stochprop prop build-pgm --atmo-dir samples/winter/ --output-path prop/winter/winter \\
             --src-loc '[30.0, -120.0, 0.0]'  --cpu-cnt 8 --local-temp-dir samples/winter/arrivals

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
    click.echo("  Atmospheric specifications directory: " + atmo_dir)
    click.echo("  Specification pattern: " + atmo_pattern)
    if topo_file is not None:
        click.echo("  Terrain file: " + topo_file)
        station_centered = True       

    if station_centered:
        click.echo("    --> note: running station-centered construction via back projection")

    click.echo("  Model output path: " + output_path)

    click.echo('\n' + "infraGA/GeoAc parameters:")
    click.echo("  Source location: " + src_loc)
    click.echo("  Inclination angles (min, max, step): " + inclinations)
    click.echo("  Azimuth angles (min, max, step): " + azimuths)
    click.echo("  Bounces: " + str(bounces))
    click.echo("  Ground elevation: " + str(z_grnd))
    click.echo("  Range max: " + str(rng_max))
    click.echo("  Frequency: " + str(freq))
    click.echo("  Local temporary directory: " + str(local_temp_dir))
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
    
    propagation.run_infraga(atmo_dir, output_path + ".arrivals.dat", pattern=atmo_pattern, cpu_cnt=cpu_cnt, geom="sph", bounces=bounces, 
                    inclinations=inclinations, azimuths=azimuths, freq=freq, z_grnd=z_grnd, rng_max=rng_max, src_loc=src_loc, infraga_path=infraga_path, 
                    local_temp_dir=local_temp_dir, prof_format=prof_format, reverse_winds=station_centered, topo_file=topo_file, verbose=verbose)

    # update source altitude if below ground surface
    src_loc[2] = max(src_loc[2], z_grnd)

    pgm = propagation.PathGeometryModel()
    pgm.build(output_path + ".arrivals.dat", output_path + ".pgm", geom="sph", src_loc=src_loc, rng_width=rng_window, rng_spacing=rng_step, rng_max=rng_max,
                az_bin_cnt=az_bin_cnt, az_bin_wdth=az_bin_width, min_turning_ht=min_turning_ht, station_centered=station_centered)



@click.command('build-tlm', short_help="Build a transmission loss model (TLM)")
@click.option("--atmo-dir", help="Directory containing atmospheric specifications", prompt="Path to directory with atmospheric specifications")
@click.option("--output-path", help="Path and prefix for TLM output", prompt="Path and prefix for TLM output")
@click.option("--atmo-pattern", help="Atmosphere file pattern (default: '*.met')", default="*.met")
@click.option("--ncpaprop-method", help="NCPAprop method ('modess' or 'epape')", default='modess')
@click.option("--ncpaprop-path", help="Path to NCPAprop binaries (if not on path)", default="")
@click.option("--freq", help="Frequency for simulation (default: 0.5 Hz)", default=0.5)
@click.option("--azimuths", help="Azimuth min, max, and step (default: [0, 360, 3]", default = "[0.0, 360.0, 3.0]")
@click.option("--z-grnd", help="Ground elevation for simulations (default: 0.0)", default=0.0)
@click.option("--rng-max", help="Maximum range for simulations (default: 1000.0)", default=1000.0)
@click.option("--rng-resol", help="Range resolution for output (default: 1.0)", default=1.0)
@click.option("--src-loc", help="Source location (lat, lon, alt)", prompt="Enter source location (lat, lon, alt)")
@click.option("--local-temp-dir", help="Local storage for individual NCPAprop results", default=None)
@click.option("--cpu-cnt", help="Number of CPUs for propagation simulations", default=None)
@click.option("--az-bin-cnt", help="Number of azimuth bins in TLM (default: 16)", default=16)
@click.option("--az-bin-width", help="Azimuth bin width in TLM (default: 30 deg)", default=30.0)
@click.option("--rng-lims", help="Range limits in TLM (default: [1, 1000])", default='[1, 1000.0]')
@click.option("--rng-cnt", help="Range intervals in TLM (default: 100)", default=100)
@click.option("--rng-spacing", help="Option for range sampling ('linear' or 'log')", default='linear')
@click.option("--use-coherent-tl", help="Use coherent transmission loss (default: False)", default=False)
@click.option("--station-centered", help="Build a station-centered model (default: False)", default=False)
@click.option("--use-topo", help="Include terrain in simulations (Default: False)", default=None)

def build_tlm(atmo_dir, output_path, atmo_pattern, ncpaprop_method, ncpaprop_path, freq, azimuths, z_grnd,
              rng_max, rng_resol, src_loc, local_temp_dir, cpu_cnt, az_bin_cnt, az_bin_width, rng_lims, rng_cnt, rng_spacing, 
              use_coherent_tl, station_centered, use_topo):
    '''
    \b
    stochprop prop build-tlm
    ------------------------
    \b
    Example Usage:
    \t stochprop prop build-tlm --atmo-dir samples/winter/ --output-path prop/winter/winter --freq 0.2  --cpu-cnt 8 --local-temp-dir samples/winter/tloss

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
    click.echo("  Atmospheric specifications directory: " + atmo_dir)
    click.echo("  Specification pattern: " + atmo_pattern)
    if use_topo:
        click.echo("Including terrain in analysis")
        click.echo("Source location: " + src_loc)
    click.echo("  Model output path: " + output_path)

    if use_topo:
        ncpaprop_method = "epape"
        station_centered = True

    click.echo('\n' + "NCPAprop parameters:")
    if len(ncpaprop_path) > 0:
        click.echo("  NCPAprop path:", ncpaprop_path)
    click.echo("NCPAprop method: " + ncpaprop_method)
    click.echo("  Frequency: " + str(freq))
    click.echo("  Azimuth angles (min, max, step): " + azimuths)
    click.echo("  Ground elevation: " + str(z_grnd))
    click.echo("  Range max, resolution: " + str(rng_max) + ", " + str(rng_resol))
    click.echo("  Local temporary directory: " + str(local_temp_dir))
    if cpu_cnt is not None:
        click.echo("  CPU count: " + str(cpu_cnt))

    click.echo('\n' + "Transmission Loss Model (TLM) parameters:")
    click.echo("  Azimuth bin count: " + str(az_bin_cnt))
    click.echo("  Azimuth bin width: " + str(az_bin_width))
    click.echo("  Range limits: " + str(rng_lims))
    click.echo("  Range count: " + str(rng_cnt))
    click.echo("  Range spacing: " + str(rng_spacing))
    click.echo("  Use coherent transmission loss: " + str(use_coherent_tl))
    if use_topo:
        click.echo("  Use terrain in epape: True")
    click.echo("")

    src_loc = [float(val) for val in src_loc.strip(' ()[]').split(',')]
    azimuths = [float(val) for val in azimuths.strip(' ()[]').split(',')]
    rng_lims = [float(val) for val in rng_lims.strip(' ()[]').split(',')]

    if cpu_cnt is None:
        cpu_cnt = 1
    else:
        cpu_cnt = int(cpu_cnt)

    propagation.run_ncpaprop(ncpaprop_method, atmo_dir, output_path + "_%.3f" %freq + "Hz", pattern=atmo_pattern, 
                             azimuths=azimuths, freq=freq, z_grnd=z_grnd, rng_max=rng_max, rng_resol=rng_resol, 
                             src_loc=src_loc, ncpaprop_path=ncpaprop_path, use_topo=use_topo, reverse_winds=station_centered,
                             local_temp_dir=local_temp_dir, cpu_cnt=cpu_cnt)

    if ncpaprop_method == "modess":
        output_suffix = ".nm"
    else:
        output_suffix = ".pe"
        
    tlm = propagation.TLossModel()
    tlm.build(output_path + "_%.3f" %freq + "Hz" + output_suffix, output_path + "_%.3f" %freq + "Hz.tlm", use_coh=use_coherent_tl, az_bin_cnt=az_bin_cnt,
                az_bin_wdth=az_bin_width, rng_lims=rng_lims, rng_cnt=rng_cnt, rng_smpls=rng_spacing)


@click.command('yld-hob', short_help="Compute statistics for atmospheric explosions")
@click.option("--infraga-config", help="InfraGA config file", prompt="infraGA config file: ")
@click.option("--output-path", help="Path and prefix for output", prompt="Path for output")
@click.option("--yld-lims", help="Yield limits [kg eq. TNT]", default="1.0e3, 10.0e6")
@click.option("--yld-cnt", help="Yield values to consider (log scaling)", default=50)
@click.option("--hob-lims", help="Height-of-burst limits [km]", default="0.0, 50.0")
@click.option("--hob-dz", help="Height-of-burst resolution [km]", default=1.0)
@click.option("--channel-cnt", help="Number of channels in array", default=6)
@click.option("--local-temp-dir", help="Local storage for individual infraGA results", default=None)
def yld_hob(infraga_config, output_path, yld_lims, yld_cnt, hob_lims, hob_dz, channel_cnt, local_temp_dir):
    '''
    \b
    stochprop prop yld-hob
    ---------------------
    \b
    Example Usage:
    \t stochprop prop yld-hob --infraga-config yld-hob_test.cnfg --output-path yld-hob1 --local-temp-dir yld-hob_temp

    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##            stochprop            ##")
    click.echo("##       Propagation Methods       ##")
    click.echo("##      Yield vs. HOB Analysis     ##")
    click.echo("##                                 ##")
    click.echo("######################################")
    click.echo("")  


    click.echo('\n' + "Parameter summary:")
    click.echo("  Output path: " + output_path)
    click.echo("  Local temporary directory: " + str(local_temp_dir))

    click.echo("\n  InfraGA config file: " + infraga_config)

    infraga_cnfg = cnfg.ConfigParser()
    infraga_cnfg.read(infraga_config)

    click.echo("    Atmosphere file: " + infraga_cnfg['GENERAL']['atmo_file'])
    click.echo("    Source lat/lon: " + infraga_cnfg['EIGENRAY']['src_lat'] + ", " + infraga_cnfg['EIGENRAY']['src_lon']) 
    click.echo("    Receiver lat/lon: " + infraga_cnfg['EIGENRAY']['rcvr_lat'] + ", " + infraga_cnfg['EIGENRAY']['rcvr_lon']) 

    click.echo("\n  Grid info:")
    click.echo("    Yield lims [kg eq TNT]: " + yld_lims)
    click.echo("    Yield count: " + str(yld_cnt))

    click.echo("    HOB lims [km]: " + hob_lims)
    click.echo("    HOB resolution [km]: " + str(hob_dz))

    yld_lims = [float(val) for val in yld_lims.strip(' ()[]').split(',')]
    hob_lims = [float(val) for val in hob_lims.strip(' ()[]').split(',')]

    alt_vals = np.arange(hob_lims[0], hob_lims[1] + hob_dz / 2.0, hob_dz)
    yld_vals = np.logspace(np.log10(yld_lims[0]), np.log10(yld_lims[1]), yld_cnt)

    propagation.yield_hob_stats(yld_vals, alt_vals, infraga_config, output_path, channel_cnt, local_temp_dir=local_temp_dir)

