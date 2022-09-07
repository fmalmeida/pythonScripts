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
import os
import yaml
from pathlib import Path
from .utils import find_files
from .general_stats_function import *
from .plasmid_function import *
from .virulence_function import *
from .resistance_function import *

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

    # get input absolute path
    bacannot_dir = '/'.join(os.path.abspath(available_annotations[0]).split('/')[:-2])

    # initiate first dictionary level
    for sample in available_samples:
        bacannot_summary[sample] = {}
        bacannot_summary[sample]['results_dir'] = f"{bacannot_dir}/{sample}"
    
    return bacannot_dir,bacannot_summary

#######################################
### Def main bacannot2json function ###
#######################################
def bacannot2json(indir, outfile, check):

    # initialize
    bacannot_dir, bacannot_summary = dict_init( indir )

    # check general annotation stats
    general_stats( bacannot_summary )

    # check virulence annotation stats
    virulence_stats( bacannot_summary )

    # check plasmids annotation stats
    plasmids_stats( bacannot_summary )

    # check resistance annotation stats
    resistance_stats( bacannot_summary )

    # save results
    final_results = json.dumps( bacannot_summary, sort_keys=True, indent=4 )
    with open(outfile, 'w') as file:
        file.write(final_results)

    # keep checking
    if check:
        print( final_results )
    
    # bye bye
    print(f"==> Output generated and saved at:\n\t{outfile}")