# __main__.py
#
# stochprop CLI control methods
#  
#
# Philip Blom (pblom@lanl.gov)


import click

from . import stats as stats_cli
from . import prop as prop_cli
from . import plot as plot_cli


@click.group(context_settings={'help_option_names': ['-h', '--help']})
def main():
    '''
    \b
    stochprop
    ---------

    Python-based tools for quantifying infrasonic propagation uncertainty via stochastic analyses

    '''
    pass


@click.group('stats', short_help="Atmosphere statistics methods", context_settings={'help_option_names': ['-h', '--help']})
def stats():
    '''
    stochprop stats - Atmospheric statistics methods using Empirical Orthogonal Function (EOF) and related analysis methods.
    
    '''
    pass 


@click.group('prop', short_help="Propagation model construction methods", context_settings={'help_option_names': ['-h', '--help']})
def prop():
    '''
    stochprop prop - Construct and interact with stochastic propagation models
    
    '''
    pass 


@click.group('plot', short_help="Visualization methods", context_settings={'help_option_names': ['-h', '--help']})
def plot():
    '''
    stochprop plot - Methods to visualize various results and analyses
    
    '''
    pass 


main.add_command(stats)
main.add_command(prop)
main.add_command(plot)

stats.add_command(stats_cli.eof_build)
stats.add_command(stats_cli.eof_coeffs)
stats.add_command(stats_cli.coeff_overlap)
stats.add_command(stats_cli.sample_eofs)
stats.add_command(stats_cli.perturb)

prop.add_command(prop_cli.fit_celerity)
prop.add_command(prop_cli.build_pgm)
prop.add_command(prop_cli.build_tlm)
prop.add_command(prop_cli.yld_hob)

plot.add_command(plot_cli.ess_ratio)
plot.add_command(plot_cli.eofs)
plot.add_command(plot_cli.eof_fit)
plot.add_command(plot_cli.atmo_ensemble)
plot.add_command(plot_cli.coeff_overlap)
plot.add_command(plot_cli.prop_model)
plot.add_command(plot_cli.detection_stats)
plot.add_command(plot_cli.network_performance)


if __name__ == '__main__':
    main()

