# stochprop/cli/prop.py
#
# Command line interface for building propagation statistics
#
# Philip Blom (pblom@lanl.gov)

import click

from stochprop import propagation

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
@click.option("--clean-up", help="Remove individual results after merge (default: True)", default=True)
@click.option("--cpu-cnt", help="Number of CPUs for propagation simulations", default=None)
@click.option("--rng-window", help="Range window in PGM (default: 50 km)", default=50.0)
@click.option("--rng-step", help="Range resolution in PGM (default: 10 km)", default=10.0)
@click.option("--az-bin-cnt", help="Number of azimuth bins in PGM (default: 16)", default=16)
@click.option("--az-bin-width", help="Azimuth bin width in PGM (default: 30 deg)", default=30.0)
def build_pgm(atmos_dir, atmos_pattern, output_path, src_loc, inclinations, azimuths, bounces, z_grnd, rng_max,
                freq, clean_up, cpu_cnt, rng_window, rng_step, az_bin_cnt, az_bin_width):
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
                    inclinations=inclinations, azimuths=azimuths, freq=freq, z_grnd=z_grnd, rng_max=rng_max, src_loc=src_loc, clean_up=clean_up)

    pgm = propagation.PathGeometryModel()
    pgm.build(output_path + ".arrivals.dat", output_path + ".pgm", geom="sph", src_loc=src_loc, rng_width=rng_window, 
                rng_spacing=rng_step, az_bin_cnt=az_bin_cnt, az_bin_wdth=az_bin_width)


@click.command('build-tlm', short_help="Build a transmission loss model (TLM)")
@click.option("--atmos-dir", help="Directory containing atmospheric specifications", prompt="Path to directory with atmospheric specifications")
@click.option("--atmos-pattern", help="Atmosphere file pattern (default: '*.met')", default="*.met")
@click.option("--output-path", help="Path and prefix for TLM output", prompt="Path and prefix for TLM output")
@click.option("--freq", help="Frequency for simulation (default: 0.5 Hz)", default=0.5)
@click.option("--azimuths", help="Azimuth min, max, and step (default: [0, 360, 6]", default = "[0.0, 360.0, 6.0]")
@click.option("--z-grnd", help="Ground elevation for simulations (default: 0.0)", default=0.0)
@click.option("--rng-max", help="Maximum range for simulations (default: 1000.0)", default=1000.0)
@click.option("--clean-up", help="Remove individual results after merge (default: True)", default=True)
@click.option("--cpu-cnt", help="Number of CPUs for propagation simulations", default=None)
@click.option("--az-bin-cnt", help="Number of azimuth bins in TLM (default: 16)", default=16)
@click.option("--az-bin-width", help="Azimuth bin width in TLM (default: 30 deg)", default=30.0)
@click.option("--rng-lims", help="Range limits in TLM (default: [1, 1000])", default='[1, 1000.0]')
@click.option("--rng-cnt", help="Range intervals in TLM (default: 100)", default=100)
@click.option("--rng-spacing", help="Option for range sampling ('linear' or 'log')", default='linear')
@click.option("--use-coherent-tl", help="Use coherent transmission loss (default: False", default=False)
def build_tlm(atmos_dir, atmos_pattern, output_path, freq, azimuths, z_grnd, rng_max, clean_up, cpu_cnt, az_bin_cnt, 
                az_bin_width, rng_lims, rng_cnt, rng_spacing, use_coherent_tl):
    '''
    \b
    stochprop prop build-tlm 
    ---------------------
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
    click.echo("  Model output path: " + output_path)

    click.echo('\n' + "NCPAprop modess parameters:")
    click.echo("  Frequency: " + str(freq))
    click.echo("  Azimuth angles (min, max, step): " + azimuths)
    click.echo("  Ground elevation: " + str(z_grnd))
    click.echo("  Range max: " + str(rng_max))
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

    propagation.run_modess(atmos_dir, output_path + "_%.3f" %freq + "Hz", pattern=atmos_pattern, azimuths=azimuths, freq=freq, 
                            z_grnd=z_grnd, rng_max=rng_max, clean_up=clean_up, cpu_cnt=cpu_cnt)

    tlm = propagation.TLossModel()
    tlm.build(output_path + "_%.3f" %freq + "Hz.nm", output_path + "_%.3f" %freq + "Hz.tlm", use_coh=use_coherent_tl, az_bin_cnt=az_bin_cnt,
                az_bin_wdth=az_bin_width, rng_lims=rng_lims, rng_cnt=rng_cnt, rng_smpls=rng_spacing)

