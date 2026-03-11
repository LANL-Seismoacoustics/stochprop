# __main__.py
#
# stochprop CLI control methods
#  
#
# Philip Blom (pblom@lanl.gov)

import os 
import click
import webbrowser
import subprocess
import shlex

from importlib.util import find_spec 

from . import stats as stats_cli
from . import prop as prop_cli
from . import plot as plot_cli
from . import utils as utils_cli


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

@click.group('utils', short_help="Utility methods", context_settings={'help_option_names': ['-h', '--help']})
def utils():
    '''
    stochprop utils - Utility and support functionality for stochprop
    
    '''
    pass 


@click.group('utils', short_help="Utility functions", context_settings={'help_option_names': ['-h', '--help']})
def utils():
    '''
    stochprop utils - Utility functions supplementing statistics and other methods
    
    '''
    pass 


#######################
##    Open Manual    ##
#######################
@click.command('doc', short_help="Open stochprop manual")
def open_doc():

    pkg_loc = find_spec('stochprop').submodule_search_locations[0]
    filename = pkg_loc + '/doc/build/html/index.html'

    if not os.path.isfile(filename):      
        print("Compiling manual...")
        subprocess.run(shlex.split("make html -C " + pkg_loc + "/doc/"), shell=False)

    webbrowser.open('file://' + os.path.realpath(filename), new=2)


main.add_command(open_doc)
main.add_command(stats)
main.add_command(prop)
main.add_command(plot)
main.add_command(utils)

stats.add_command(stats_cli.eof_build)
stats.add_command(stats_cli.eof_coeffs)
stats.add_command(stats_cli.coeff_overlap)
stats.add_command(stats_cli.sample_eofs)
stats.add_command(stats_cli.perturb)

# prop.add_command(prop_cli.fit_celerity)
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

utils.add_command(utils_cli.eig_wvfrm2json)
utils.add_command(utils_cli.celerity_gmm)


if __name__ == '__main__':
    main()

