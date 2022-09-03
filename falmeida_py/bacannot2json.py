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
    falmeida-py bacannot2json [ --input <bacannot_results> --output <outfile> ]

Options:
    -h --help                         Show this screen.
    -i --input=<bacannot_results>     Path to bacannot results folder.
    -o --output=<outfile>             JSON summary output file [Default: bacannot_summary.json].
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
from .plasmid_function import *
from .virulence_function import *
from .resistance_function import *

###############################################
### based on annotations figure sample name ###
###############################################
def get_samples(filepaths):
    samples = []
    for annotation in filepaths:
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

################################
### check general annotation ###
################################
def general_stats(bacannot_summary):

    # iterate over available samples
    for sample in bacannot_summary:

        # load dir of samples' results
        results_dir = bacannot_summary[sample]['results_dir']

        # load annotation stats
        general_results = yaml.safe_load(
            Path(f"{results_dir}/annotation/{sample}.txt").read_text()
        )

        # load MLST
        mlst_results = pd.read_csv(
            f"{results_dir}/MLST/{sample}_mlst_analysis.txt",
            sep='\t', header=None
        )

        # load refseq_masher
        refseq_masher_results = pd.read_csv(
            f"{results_dir}/refseq_masher/refseq_masher_results.txt",
            sep='\t'
        )
        refseq_masher_results.sort_values(by='distance', ascending=True, inplace=True)

        # save annotation stats
        bacannot_summary[sample]['general_annotation'] = {}
        bacannot_summary[sample]['general_annotation']['mlst']  = mlst_results[2].item().replace('-', 'NaN')
        bacannot_summary[sample]['general_annotation']['cds']   = general_results['CDS']
        bacannot_summary[sample]['general_annotation']['rrna']  = general_results['rRNA']
        bacannot_summary[sample]['general_annotation']['trna']  = general_results['tRNA']
        bacannot_summary[sample]['general_annotation']['tmrna'] = general_results['tmRNA']

        bacannot_summary[sample]['general_annotation']['closest_reference'] = {}
        bacannot_summary[sample]['general_annotation']['closest_reference']['strain'] = refseq_masher_results.head(1)['top_taxonomy_name'].item()
        bacannot_summary[sample]['general_annotation']['closest_reference']['distance'] = refseq_masher_results.head(1)['distance'].item()
        bacannot_summary[sample]['general_annotation']['closest_reference']['accession'] = refseq_masher_results.head(1)['assembly_accession'].item()

#######################################
### Def main bacannot2json function ###
#######################################
def bacannot2json(indir, outfile):

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
    print( final_results )