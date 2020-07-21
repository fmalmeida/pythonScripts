#!/usr/bin/env python3
# coding: utf-8

## Def help message
"""
A simple script to filter a GFF file based on a list of gene ids

This was created using GFF field names that specifically matches the GFF fields
from the BNUT genome annotation file

If using for another dataset, please check the script to precisely match your GFF fields

---
Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain
Usage:
    filter_gff.py
    filter_gff.py -h|--help
    filter_gff.py -v|--version
    filter_gff.py [--input <gff> --fofn <file>]

Options:
    -h --help                   Show this screen.
    -v --version                Show version information
    --input=<gff>               GFF file for subset
    --fofn=<file>               Gene list for subsetting input GFF
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt
import pandas as pd
import re



##############################################
### Function for filtering 9th column list ###
##############################################

def Filter(string, substr):
    return [str for str in string if
             any(sub in str for sub in substr)]


########################################
### Function for execution of parser ###
########################################

def filter_gff(input_gff, gene_list):

    # Read GFF file
    my_df = pd.read_csv(input_gff, sep = "\t", comment = "#",
                        names=['Chr', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'])
    my_df_smaller = my_df[~my_df['Type'].isin(['gene','mRNA'])]

    # Create test list
    gene_ids = open(gene_list).readlines()
    gene_ids = list(map(str.strip, gene_ids))


    # Parse
    for index, line in my_df.iterrows():

        # Gene Features  ---  Simple check
        if line['Type'] == 'gene':
            id = Filter(line['Attributes'].split(';'), ['ID'])[0].split('=')[1].split('.g')[0]
            if id in gene_ids:
                print(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], sep="\t")

        # mRNA features --- Check based on Parent gene ID
        # Also extract mRNA id
        # and create a pacid field required for Shiu's pipeline

        elif line['Type'] == 'mRNA':
            parent     = Filter(line['Attributes'].split(';'), ['Parent'])[0].split('=')[1].split('.g')[0]
            pacid      = Filter(line['Attributes'].split(';'), ['ID'])[0].split(':')[1]
            pacid_full = Filter(line['Attributes'].split(';'), ['ID'])[0]
            if parent in gene_ids:
                print(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8] + ";pacid=" + pacid, sep="\t")
                subset    = my_df_smaller[my_df_smaller['Attributes'].str.contains(str(pacid_full))]

                # Add pacid field to each subset line
                # And print

                for index, row in subset.iterrows():
                    print(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8] + ";pacid=" + pacid, sep="\t")




## Main
if __name__ == '__main__':
    arguments = docopt(__doc__, version='v1.0 by Felipe Marques de Almeida')

    ## Single GFF
    if arguments['--input'] and arguments['--fofn']:
        print("##gff-version 3")
        filter_gff(input_gff=arguments['--input'], gene_list=arguments['--fofn'])

    ## None
    else:
        print("Missing mandatory arguments")
        print("Please, check out the help message")
        print("")
        print(arguments)
