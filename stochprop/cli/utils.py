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

import gzip
import json

import numpy as np

from obspy import UTCDateTime

from infrapy.utils import data_io

from stochprop import utils as sp_utils




@click.command('eig-wvfrm2json', short_help="Create JSON detections from infraGA results")
@click.option("--wvfrm-file", help="InfraGA waveform output", prompt="infraGA waveform output: ")
@click.option("--output-label", help="Output label to write JSON file", prompt="Detection JSON file label: ")
@click.option("--origin-time", help="Origin time of synthetic event", prompt="Event origin time: ")
@click.option("--sta", help="Station code for synthetic detection(s)", prompt="Station code: ")
@click.option("--loc", help="Station code for synthetic detection(s)", default="")
@click.option("--cha", help="Station code for synthetic detection(s)", default="BDF")
@click.option("--c0", help="Sound speed to relate inclination to trace velocity", default=340.0)
@click.option("--ns-option", help="'min', 'mid', or 'max' IMS noise statistic", default="min")
@click.option("--array-dim", help="Array dimension (station size)", default=6)
def eig_wvfrm2json(wvfrm_file, output_label, origin_time, sta, loc, cha, c0, ns_option, array_dim):
    '''
    \b
    stochprop plot eof-fit
    -----------------------
    \b
    Example Usage:
    \t stochprop utils eig-wvfrm2json --wvfrm-file test.wvfrms.dat --output-label "test" --origin-time "2010-01-01T12:10:00" --sta "STA"
    
    '''

    click.echo("")
    click.echo("########################################")
    click.echo("##                                    ##")
    click.echo("##             stochprop              ##")
    click.echo("##           Utility Tools            ##")
    click.echo("##  infraGA eig-wvfrm 2 infraPy JSON  ##")
    click.echo("##                                    ##")
    click.echo("########################################")
    click.echo("")  

    click.echo('\n' + "Data I/O and Configuration Summary:")
    click.echo("  infraGA eigenray + waveform results file: " + wvfrm_file)
    click.echo("  Output JSON detection file: " + output_label + ".dets.json.gz")
    click.echo("  Event origin time: " + origin_time)
    click.echo("  Reference trace ID: SY." + sta + "." + loc + "." + cha)
    click.echo("  Ground sound speed [m/s]: " + str(c0))
    click.echo("  IMS noise option ('min', 'mid', or 'max'): " + ns_option)
    click.echo("  Array dimension: " + str(array_dim))
    click.echo("")


    dets_dict = sp_utils.eig_wvfrm2json(wvfrm_file, UTCDateTime(origin_time), c0, ns_option, array_dim, sta, loc, cha)
    click.echo("Writing " + str(len(dets_dict)) + " detections into JSON file." + '\n')

    det_output = {'det_info' : dets_dict}
    with gzip.open(output_label + ".dets.json.gz", 'wt', encoding='UTF-8') as zipfile:
        json.dump(det_output, zipfile, indent=4, cls=data_io.Infrapy_Encoder)


