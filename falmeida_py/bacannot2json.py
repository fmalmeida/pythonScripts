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
            
            # vfdb
            if os.path.exists(f"{results_dir}/virulence/vfdb/{sample}_vfdb_blastn_onGenes.summary.txt"):

                # init VFDB annotation dictionary
                bacannot_summary[sample]['virulence_annotation']['VFDB'] = {}

                # load vfdb results
                results = pd.read_csv(
                    f"{results_dir}/virulence/vfdb/{sample}_vfdb_blastn_onGenes.summary.txt",
                    sep='\t'
                )

                # number of virulence genes
                total_number_of_genes = len( results['SEQUENCE'].unique() )
                bacannot_summary[sample]['virulence_annotation']['VFDB']['total'] = total_number_of_genes

                # identified VFs
                detected_vf_names = [x.split('_(')[0].replace("[", "") for x in results['PRODUCT'].unique() ]
                detected_vf_ids   = [x.split('_(')[1].split(')')[0].replace(")", "") for x in results['PRODUCT'].unique() ]
                detected_vf_genes = [x.replace(")", "").replace("(", "") for x in results['GENE'].unique() ]

                bacannot_summary[sample]['virulence_annotation']['VFDB']['detected_vf_names'] = ';'.join( detected_vf_names )
                bacannot_summary[sample]['virulence_annotation']['VFDB']['detected_vf_ids'] = ';'.join( detected_vf_ids )
                bacannot_summary[sample]['virulence_annotation']['VFDB']['detected_vf_genes'] = ';'.join( detected_vf_genes )
            
            # victors
            if os.path.exists(f"{results_dir}/virulence/victors/{sample}_victors_blastp_onGenes.summary.txt"):

                # init victors annotation dictionary
                bacannot_summary[sample]['virulence_annotation']['Victors'] = {}

                # load victors results
                results = pd.read_csv(
                    f"{results_dir}/virulence/victors/{sample}_victors_blastp_onGenes.summary.txt",
                    sep='\t'
                )

                # number of virulence genes
                total_number_of_genes = len( results['SEQUENCE'].unique() )
                bacannot_summary[sample]['virulence_annotation']['Victors']['total'] = total_number_of_genes

                # identified VFs
                detected_vf_ids   = [x.replace("Victors_", "") for x in results['VICTORS_ID'].unique() ]
                detected_vf_genes = [x for x in results['GENE'].unique() ]

                bacannot_summary[sample]['virulence_annotation']['Victors']['detected_vf_ids'] = ';'.join( detected_vf_ids )
                bacannot_summary[sample]['virulence_annotation']['Victors']['detected_vf_genes'] = ';'.join( detected_vf_genes )

########################################
### check plasmids annotations stats ###
########################################
def plasmids_stats(bacannot_summary):
    
    # iterate over available samples
    for sample in bacannot_summary:

        # load dir of samples' results
        results_dir = bacannot_summary[sample]['results_dir']

        # load annotation stats
        if os.path.exists(f"{results_dir}/plasmids"):

            # init plasmids annotation dictionary
            bacannot_summary[sample]['plasmid_annotation'] = {}
            
            # platon
            if os.path.exists(f"{results_dir}/plasmids/platon/{sample}.tsv"):

                # init platon annotation dictionary
                bacannot_summary[sample]['plasmid_annotation']['platon'] = {}

                # load platon results
                results = pd.read_csv(
                    f"{results_dir}/plasmids/platon/{sample}.tsv",
                    sep='\t'
                )

                # number of plasmid annotations
                total_number = len(results.index)
                bacannot_summary[sample]['plasmid_annotation']['platon']['total'] = total_number

                # contigs that are plasmids
                contigs = [x for x in results['ID'].unique()]
                bacannot_summary[sample]['plasmid_annotation']['platon']['contigs'] = ';'.join( contigs )

                # per plasmid info
                for seq in contigs:
                    bacannot_summary[sample]['plasmid_annotation']['platon'][seq] = {}
                    bacannot_summary[sample]['plasmid_annotation']['platon'][seq]['ORFs'] = results.loc[results['ID'] == seq, '# ORFs'].item()
                    bacannot_summary[sample]['plasmid_annotation']['platon'][seq]['Circular'] = results.loc[results['ID'] == seq, 'Circular'].item()
                    bacannot_summary[sample]['plasmid_annotation']['platon'][seq]['AMRs'] = results.loc[results['ID'] == seq, '# AMRs'].item()
                    bacannot_summary[sample]['plasmid_annotation']['platon'][seq]['Replication'] = results.loc[results['ID'] == seq, '# Replication'].item()
                    bacannot_summary[sample]['plasmid_annotation']['platon'][seq]['Mobilization'] = results.loc[results['ID'] == seq, '# Mobilization'].item()
                    bacannot_summary[sample]['plasmid_annotation']['platon'][seq]['Conjugation'] = results.loc[results['ID'] == seq, '# Conjugation'].item()
            
            # plasmidfinder
            if os.path.exists(f"{results_dir}/plasmids/plasmidfinder/results_tab.tsv"):

                # init platon annotation dictionary
                bacannot_summary[sample]['plasmid_annotation']['plasmidfinder'] = {}

                # load platon results
                results = pd.read_csv(
                    f"{results_dir}/plasmids/plasmidfinder/results_tab.tsv",
                    sep='\t'
                )

                # databases
                bacannot_summary[sample]['plasmid_annotation']['plasmidfinder']['Database'] = results['Database'].unique().item()

                # number of plasmid annotations
                total_number = len(results['Contig'].unique())
                bacannot_summary[sample]['plasmid_annotation']['plasmidfinder']['total'] = total_number

                # identified contigs
                contigs = results['Contig'].unique()
                bacannot_summary[sample]['plasmid_annotation']['plasmidfinder']['contigs'] = ';'.join( contigs )

                # important info over contigs
                for seq in contigs:
                    bacannot_summary[sample]['plasmid_annotation']['plasmidfinder'][seq] = {}
                    bacannot_summary[sample]['plasmid_annotation']['plasmidfinder'][seq]['Inc types'] = ';'.join( results.loc[results['Contig'] == seq, 'Plasmid'].unique() )
                    bacannot_summary[sample]['plasmid_annotation']['plasmidfinder'][seq]['Identity'] = ';'.join( [str(x) for x in results.loc[results['Contig'] == seq, 'Identity'].unique() ] )
                    bacannot_summary[sample]['plasmid_annotation']['plasmidfinder'][seq]['Accessions'] = ';'.join( results.loc[results['Contig'] == seq, 'Accession number'].unique() )

#######################################
### Def main bacannot2json function ###
#######################################
def bacannot2json(indir, outfile):

    # initialize
    bacannot_dir,bacannot_summary = dict_init( indir )

    # check general annotation stats
    general_stats( bacannot_summary )

    # check virulence annotation stats
    virulence_stats( bacannot_summary )

    # check plasmids annotation stats
    plasmids_stats( bacannot_summary )

    # save results
    final_results = json.dumps( bacannot_summary, sort_keys=False, indent=4 )
    with open(outfile, 'w') as file:
        file.write(final_results)

    # keep checking
    print( final_results )