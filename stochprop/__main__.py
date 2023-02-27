# __main__.py
#
# stochprop CLI control methods
#  
#
# Philip Blom (pblom@lanl.gov)


import click

from . import cli 


@click.group(context_settings={'help_option_names': ['-h', '--help']})
def main():
    '''
    \b
    stochprop
    ---------

    Python-based tools for quantifying infrasonic propagation uncertainty via stochastic analyses

    '''
    pass


@click.group('eof', short_help="Empirical Orthogonal Function (EOF) methods", context_settings={'help_option_names': ['-h', '--help']})
def eof():
    '''
    stochprop eof - Analysis methods using Empirical Orthogonal Function (EOFs) to identify seasonal trends and sample the atmospheric variability
    
    '''
    pass 


@click.group('prop', short_help="Propagation model construction methods", context_settings={'help_option_names': ['-h', '--help']})
def prop():
    '''
    stochprop prop - Construct and interact with stochastic propagation models
    
    '''
    pass 


@click.group('perturb', short_help="Atmospheric specification perturbing methods", context_settings={'help_option_names': ['-h', '--help']})
def perturb():
    '''
    stochprop perturb - Generate perturbed atmospheric specifications from a reference file
    
    '''
    pass 


main.add_command(eof)
main.add_command(prop)
main.add_command(perturb)

eof.add_command(cli.eof_build)
eof.add_command(cli.eof_coeffs)
eof.add_command(cli.eof_seasonality)
eof.add_command(cli.eof_sample)

prop.add_command(cli.build_pgm)
prop.add_command(cli.build_tlm)
prop.add_command(cli.plot_model)
prop.add_command(cli.plot_detection_stats)
prop.add_command(cli.plot_network_performance)

perturb.add_command(cli.eof_perturb)
perturb.add_command(cli.gravity_waves)


if __name__ == '__main__':
    main()

