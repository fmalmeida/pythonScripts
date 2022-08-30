#!/usr/bin/env python3
# coding: utf-8

########################
### Def help message ###
########################
usage_bacannot2json = """
A script meant to subset a genbank annotation file based on alignments against a query (Nucleotide) FASTA file

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
from .utils import find_files
import json
import os
import yaml
from pathlib import Path

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

        # save annotation stats
        bacannot_summary[sample]['general_annotation'] = {}
        bacannot_summary[sample]['general_annotation']['CDS'] = general_results['CDS']
        bacannot_summary[sample]['general_annotation']['rRNA'] = general_results['rRNA']
        bacannot_summary[sample]['general_annotation']['tRNA'] = general_results['tRNA']
        bacannot_summary[sample]['general_annotation']['tmRNA'] = general_results['tmRNA']

##################################
### check virulence annotation ###
##################################
def virulence_stats(bacannot_summary):

    # iterate over available samples
    for sample in bacannot_summary:

        # load dir of samples' results
        results_dir = bacannot_summary[sample]['results_dir']

        # load annotation stats
        if os.path.exists(f"{results_dir}/virulence"):

            # init virulence annotation dictionary
            bacannot_summary[sample]['virulence_annotation'] = {}
            
            if os.path.exists(f"{results_dir}/virulence/vfdb/{sample}_vfdb_blastn_onGenes.summary.txt"):

                # init VFDB annotation dictionary
                bacannot_summary[sample]['virulence_annotation']['VFDB'] = {}

                # load vfdb results
                vfdb_results = pd.read_csv(
                    f"{results_dir}/virulence/vfdb/{sample}_vfdb_blastn_onGenes.summary.txt",
                    sep='\t'
                )

                # number of virulence genes
                total_number_of_genes = len( vfdb_results['SEQUENCE'].unique() )
                bacannot_summary[sample]['virulence_annotation']['VFDB']['total'] = total_number_of_genes

                # identified VFs
                detected_vf_names = [x.split('_(')[0].replace("[", "") for x in vfdb_results['PRODUCT'].unique() ]
                detected_vf_ids   = [x.split('_(')[1].split(')')[0].replace(")", "") for x in vfdb_results['PRODUCT'].unique() ]
                detected_vf_genes = [x.replace(")", "").replace("(", "") for x in vfdb_results['GENE'].unique() ]

                bacannot_summary[sample]['virulence_annotation']['VFDB']['detected_vf_names'] = ';'.join( detected_vf_names )
                bacannot_summary[sample]['virulence_annotation']['VFDB']['detected_vf_ids'] = ';'.join( detected_vf_ids )
                bacannot_summary[sample]['virulence_annotation']['VFDB']['detected_vf_genes'] = ';'.join( detected_vf_genes )

                # check
                print( vfdb_results )

#######################################
### Def main bacannot2json function ###
#######################################
def bacannot2json(indir, outfile):

    # initialize
    bacannot_dir,bacannot_summary = dict_init(indir)

    # check general annotation stats
    general_stats( bacannot_summary )

    # check virulence annotation stats
    virulence_stats ( bacannot_summary )

    # keep checking
    print(
        json.dumps( bacannot_summary, sort_keys=True, indent=4 )
    )