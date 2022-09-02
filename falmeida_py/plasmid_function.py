##################################
### Loading Necessary Packages ###
##################################
import pandas as pd
from .utils import find_files
import json
import os
import yaml
from pathlib import Path

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
            bacannot_summary[sample]['plasmid'] = {}
            
            # platon
            if os.path.exists(f"{results_dir}/plasmids/platon/{sample}.tsv"):

                # init platon annotation dictionary
                bacannot_summary[sample]['plasmid']['platon'] = {}

                # load platon results
                results = pd.read_csv(
                    f"{results_dir}/plasmids/platon/{sample}.tsv",
                    sep='\t'
                )

                # number of plasmid annotations
                total_number = len(results.index)
                bacannot_summary[sample]['plasmid']['platon']['total'] = total_number

                # per plasmid info
                for seq in [x for x in results['ID'].unique()]:
                    bacannot_summary[sample]['plasmid']['platon'][seq] = {}
                    bacannot_summary[sample]['plasmid']['platon'][seq]['ORFs'] = results.loc[results['ID'] == seq, '# ORFs'].item()
                    bacannot_summary[sample]['plasmid']['platon'][seq]['Circular'] = results.loc[results['ID'] == seq, 'Circular'].item()
                    bacannot_summary[sample]['plasmid']['platon'][seq]['AMRs'] = results.loc[results['ID'] == seq, '# AMRs'].item()
                    bacannot_summary[sample]['plasmid']['platon'][seq]['Replication'] = results.loc[results['ID'] == seq, '# Replication'].item()
                    bacannot_summary[sample]['plasmid']['platon'][seq]['Mobilization'] = results.loc[results['ID'] == seq, '# Mobilization'].item()
                    bacannot_summary[sample]['plasmid']['platon'][seq]['Conjugation'] = results.loc[results['ID'] == seq, '# Conjugation'].item()
            
            # plasmidfinder
            if os.path.exists(f"{results_dir}/plasmids/plasmidfinder/results_tab.tsv"):

                # init platon annotation dictionary
                bacannot_summary[sample]['plasmid']['plasmidfinder'] = {}

                # load platon results
                results = pd.read_csv(
                    f"{results_dir}/plasmids/plasmidfinder/results_tab.tsv",
                    sep='\t'
                )

                # databases
                bacannot_summary[sample]['plasmid']['plasmidfinder']['meta'] = {}
                bacannot_summary[sample]['plasmid']['plasmidfinder']['meta']['database'] = results['Database'].unique().item()

                # number of plasmid annotations
                total_number = len(results['Contig'].unique())
                bacannot_summary[sample]['plasmid']['plasmidfinder']['total'] = total_number

                # plasmid annotations contigs
                for seq in [ str(x) for x in results['Contig'].unique() ]:
                    bacannot_summary[sample]['plasmid']['plasmidfinder'][seq] = {}
                    bacannot_summary[sample]['plasmid']['plasmidfinder'][seq]['inc_types'] = {}
                    bacannot_summary[sample]['plasmid']['plasmidfinder'][seq]['identity'] = {}
                    bacannot_summary[sample]['plasmid']['plasmidfinder'][seq]['accession'] = {}
                for index, row in results.iterrows():
                    contig = row['Contig']
                    bacannot_summary[sample]['plasmid']['plasmidfinder'][contig]['inc_types'] = row['Plasmid']
                    bacannot_summary[sample]['plasmid']['plasmidfinder'][contig]['identity'] = row['Identity']
                    bacannot_summary[sample]['plasmid']['plasmidfinder'][contig]['accession'] = row['Accession number']