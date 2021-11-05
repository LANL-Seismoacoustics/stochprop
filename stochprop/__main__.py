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

import click

import numpy as np

from . import eofs
from . import gravity_waves as grav

@click.group(context_settings={'help_option_names': ['-h', '--help']})
def main():
    '''
    stochprop - Python-based tools for quantifying infrasonic propagation uncertainty via stochastic analyses

    More detailed description...
    '''
    pass

@main.command('build-eofs', short_help="Build EOF info through SVD")
@click.option("--atmo-dir", help="Directory of atmspheric specifications (required)", prompt="Atmospheric specifications: ")
@click.option("--eofs-path", help="EOF output path and prefix (required)", prompt="Output path: ")
@click.option("--atmo-pattern", help="Specification file pattern (default: '*.met')", default='*.met')
@click.option("--atmo-format", help="Specification format (default: 'zTuvdp'", default='zTuvdp')
@click.option("--save-datetime", help="Save date time info (default: False)", default=False)
@click.option("--eof-cnt", help="Number of EOFs to store (default: 100)", default=100)
def build_eofs(atmo_dir, eofs_path, atmo_pattern, atmo_format, save_datetime, eof_cnt):
    '''
    stochprop build-eofs

    --------------------

    EOF construction methods...
    
    '''

    if save_datetime:
        A, z0, datetimes = eofs.build_atmo_matrix(atmo_dir, atmo_pattern, atmo_format=atmo_format, return_datetime=True)
        np.save(eofs_path + ".datetimes", datetimes)
    else:
        A, z0 = eofs.build_atmo_matrix(atmo_dir, atmo_pattern, return_datetime=False)
    eofs.compute_eofs(A, z0, eofs_path, eof_cnt=eof_cnt)
    


@main.command('eof-coeffs', short_help="Compute EOF coefficients")
@click.option("--atmo-dir", help="Directory of atmspheric specifications (required)", prompt="Reference atmospheric specification: ")
def eof_coeffs():
    '''
    stochprop eof-coeffs

    --------------------

    EOF coefficient calculation methods...(not on CLI yet)
    
    '''

    print("Compute EOF coefficients...not yet implemented.")


@main.command('eof-perturbation', short_help="Use EOFs to perturb a specification")
@click.option("--atmo-file", help="Reference atmspheric specification (required)", prompt="Reference atmospheric specification: ")
@click.option("--eofs-path", help="Path to EOF info (required)", prompt="EOF results: ")
@click.option("--out", help="Output prefix (required)", prompt="Output prefix: ")
@click.option("--std-dev", help="Standard deviation (default: 10 m/s)", default=10.0)
@click.option("--eof-max", help="Maximum EOF coefficient to use (default: 100)", default=100)
@click.option("--eof-cnt", help="Number of EOFs to use (default: 50)", default=50)
@click.option("--sample-cnt", help="Number of perturbed samples (default: 25)", default=25)
@click.option("--alt-weight", help="Altitude weighting power (default: 2.0)", default=2.0)
@click.option("--singular-value-weight", help="Sing. value weighting power (default: 0.25)", default=0.25)
def eof_perturbation(atmo_file, eofs_path, out, std_dev, eof_max, eof_cnt, sample_cnt, alt_weight, singular_value_weight):
    '''
    stochprop eof-perturbation

    --------------------

    EOF construction methods...
    
    '''
    
    eofs.perturb_atmo(atmo_file, eofs_path, out, stdev=std_dev, eof_max=eof_max, eof_cnt=eof_cnt, sample_cnt=sample_cnt, alt_wt_pow=alt_weight, sing_val_wt_pow=singular_value_weight)


@main.command('gravity-waves', short_help="Construct gravity wave pertuations")
@click.option("--atmo-file", help="Reference atmspheric specification (required)", prompt="Reference atmospheric specification: ")
@click.option("--out", help="Output prefix (required)", prompt="Output prefix: ")
@click.option("--sample-cnt", help="Number of perturbated samples (default: 25)", default=25)
@click.option("--t0", help="Propagation time from source [hr] (default: 8)", default=8.0)
@click.option("--dx", help="Horizontal wavenumber scale [km] (default: 2.0)", default=2.0)
@click.option("--dz", help="Altitude resolution for integration [km] (default: 0.2)", default=0.2)
@click.option("--nk", help="Horizontal wavenumber resolution (default: 128)", default=128)
@click.option("--nom", help="Frequency resolution (default: 5)", default=5)
@click.option("--random-phase", help="Randomize phase at source [bool] (default: False)", default=False)
@click.option("--z-src", help="Gravity wave source altitude [km] (default: 20.0)", default=20.0)
@click.option("--m-star", help="Gravity wave source spectrum peak [km] (default: 2.5 / (2 pi))", default=2.0 * np.pi / 2.5)
@click.option("--taper-below", help="Taper perturbation below source height [bool] (default: True)", default=True)
@click.option("--cpu-cnt", help="Number of CPUs to use in parallel analysis (default: None)", default=None, type=int)
@click.option("--debug-fig", help="Output for figures to aid in debugging (default: None)", default=None, type=str)
def gravity_waves(atmo_file, out, sample_cnt, t0, dx, dz, nk, nom, random_phase, z_src, m_star, taper_below, cpu_cnt, debug_fig):
    '''
    stochprop gravity-waves

    -----------------------

    Gravity wave perturbation methods based on Drob et al. (2013) method.

    More info...
    
    '''

    click.echo("")
    click.echo("#####################################")
    click.echo("##                                 ##")
    click.echo("##            stochprop            ##")
    click.echo("##    Gravity Wave Perturbations   ##")
    click.echo("##                                 ##")
    click.echo("#####################################")
    click.echo("") 

    grav.perturb_atmo(atmo_file, out, sample_cnt=sample_cnt, t0=t0 * 3600.0, dx=dx, dz=dz, Nk=nk, N_om=nom, 
            random_phase=random_phase, z_src=z_src, m_star=m_star, env_below=taper_below, cpu_cnt=cpu_cnt, fig_out=debug_fig)




if __name__ == '__main__':
    main()

