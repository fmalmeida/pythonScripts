#!/usr/bin/env python3
# coding: utf-8

########################
### Def help message ###
########################
usage_bacannot2json = """
A script to summarize the main annotation results of fmalmeida/bacannot pipeline as a structured JSON file.

---
Copyright (C) 2022 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain

Usage:
    falmeida-py bacannot2json [ -h|--help ]
    falmeida-py bacannot2json [ --input <bacannot_results> --output <outfile> --print ]

Options:
    -h --help                         Show this screen.
    -i --input=<bacannot_results>     Path to bacannot results folder.
    -o --output=<outfile>             JSON summary output file [Default: bacannot_summary.json].
    -p --print                        Also print resolved JSON to stdout.
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt
import pandas as pd
import json
import simplejson
import os
import yaml
from pathlib import Path
from .utils import find_files
from .general_stats_function import *
from .plasmid_function import *
from .virulence_function import *
from .resistance_function import *

##############################
### fix keys in dictionary ###
##############################
def stringify_keys(d):
    """Convert a dict's keys to strings if they are not."""
    for key in d.keys():

        # check inner dict
        if isinstance(d[key], dict):
            value = stringify_keys(d[key])
        else:
            value = d[key]

        # convert nonstring to string if needed
        if not isinstance(key, str):
            try:
                d[str(key)] = value
            except Exception:
                try:
                    d[repr(key)] = value
                except Exception:
                    raise

            # delete old key
            del d[key]
    return d

###############################################
### based on annotations figure sample name ###
###############################################
def get_samples(filepaths):
    samples = []
    for annotation in filepaths:
        if not 'jbrowse' in str(annotation):
            sample = annotation.split('/')[-2]
            samples.append(sample)
    return samples

#############################
### initialize dictionary ###
#############################
def dict_init(indir):
    # main dictionary for json final output
    bacannot_summary = {}

    # detect available sample annotations
    available_annotations = [ str(x) for x in find_files(start_dir=indir, pattern='annotation') ]

    # detect available samples
    available_samples = [ str(x) for x in get_samples(available_annotations) ]

    # initiate first dictionary level
    for sample in available_samples:
        bacannot_summary[sample] = {}
        bacannot_summary[sample]['sample'] = f"{sample}"
        bacannot_summary[sample]['results_dir'] = f"{indir}/{sample}"
    
    return bacannot_summary

#######################################
### Def main bacannot2json function ###
#######################################
def bacannot2json(indir, outfile, check):

    # initialize
    bacannot_dir     = os.path.abspath( indir )
    bacannot_summary = dict_init( bacannot_dir )

    # check general annotation stats
    general_stats( bacannot_summary )

    # check virulence annotation stats
    virulence_stats( bacannot_summary )

    # check plasmids annotation stats
    plasmids_stats( bacannot_summary )

    # check resistance annotation stats
    resistance_stats( bacannot_summary )

    # save results
    final_results = simplejson.dumps( 
        stringify_keys( bacannot_summary ), 
        sort_keys=True, 
        indent=4, 
        ignore_nan=True 
    )
    with open(outfile, 'w') as file:
        file.write(final_results)

    # keep checking
    if check:
        print( final_results )
    
    # bye bye
    print(f"==> Output generated and saved at:\n\t{outfile}")