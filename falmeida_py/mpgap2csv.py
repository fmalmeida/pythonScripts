#!/usr/bin/env python3
# coding: utf-8

########################
### Def help message ###
########################
usage_mpgap2csv = """
A simple to generate metadata .csv for quickly producing tables for papers. Uses statistics calculated with quast and busco, condensed in multiqc file, generated with fmalmeida/MpGAP pipeline.

---
Copyright (C) 2022 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain

Usage:
    falmeida-py mpgap2csv [ -h|--help ] [ --input <indir> --output <outfile> ]

Options:
    -h --help                   Show this screen.
    --input=<indir>             Path to MpGAP outdir.
    --output=<outfile>          File to save results (CSV). [Default: MpGAP_summary.csv]
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt
import pandas as pd
import os
from pathlib import Path
import json
from pprint import pprint

###################################
### Defifining useful functions ###
###################################
def find_multiqc_file(dir):
    matches = []
    for path in Path(dir).rglob('multiqc_data.json'):
        matches.append(os.path.abspath(path.resolve()))
    return matches

def split_and_retrieve_desired_values(item):
    values=str(item).split('/')
    sample=values[-5]
    method=values[-4]
    del values[-5:]
    outdir='/'.join(values)
    
    return sample, method, outdir

def parse_json(data, assembler, field, item):
    
    return data['report_saved_raw_data'][field][assembler][item]

def get_sample_and_assembly_info(files, df, quast_cols, busco_cols, base_columns):
    
    for idx, item in enumerate(files):

        sample, method, outdir = split_and_retrieve_desired_values(item)
        
        with open(item) as json_file:
            data = json.load(json_file)
            assemblies = list(data['report_general_stats_data'][0].keys())

            # quast desires
            quast_selection = []
            for index, val in enumerate(assemblies):
                quast_selection.append( val )
                for selection in quast_cols:
                    quast_selection.append( parse_json(
                        data, val, 'multiqc_quast', selection
                    ) )
            quast_selection = [quast_selection[n:n+10] for n in range(0, len(quast_selection), 10)]

            # busco desires
            busco_selection = []
            for index, val in enumerate(assemblies):
                busco_selection.append( val )
                for selection in busco_cols:
                    busco_selection.append( parse_json(
                        data, val, 'multiqc_busco', selection
                    ) )
            busco_selection = [busco_selection[n:n+7] for n in range(0, len(busco_selection), 7)]
            
            for index, val in enumerate(assemblies):
                final = [sample, method, outdir, item, val] + quast_selection[index][1:] + busco_selection[index][1:]
                df.loc[len(df)] = final
                # print(final)
        
             

####################
## Defining main ###
####################
def mpgap2csv(indir, output):

    base_columns = [ "sample", "method", "outdir", "multiqc_file", "software" ]
    desired_quast_columns = [
        "# contigs", "N50", "Total length", 
        "# total reads", "Properly paired (%)", "Avg. coverage depth",
        "# predicted rRNA genes", "Complete BUSCO (%)", "Partial BUSCO (%)"
        ]
    desired_busco_columns = [
        "complete_single_copy", "complete_duplicated", "fragmented", 
        "missing", "total", "lineage_dataset"
    ]
    multiqc_files_df = pd.DataFrame(columns = list(base_columns + desired_quast_columns + desired_busco_columns))

    get_sample_and_assembly_info(
        find_multiqc_file(indir), 
        multiqc_files_df, 
        desired_quast_columns, 
        desired_busco_columns,
        base_columns
    )

    # save file
    multiqc_files_df.to_csv(output, index=False)
    print(f"=> Saved results in:\n\t{output}")
