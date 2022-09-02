##################################
### Loading Necessary Packages ###
##################################
import pandas as pd
from .utils import find_files
import json
import os
import yaml
from pathlib import Path

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
            bacannot_summary[sample]['virulence'] = {}
            
            # vfdb
            if os.path.exists(f"{results_dir}/virulence/vfdb/{sample}_vfdb_blastn_onGenes.summary.txt"):

                # init VFDB annotation dictionary
                bacannot_summary[sample]['virulence']['VFDB'] = {}

                # load vfdb results
                results = pd.read_csv(
                    f"{results_dir}/virulence/vfdb/{sample}_vfdb_blastn_onGenes.summary.txt",
                    sep='\t'
                )

                # number of virulence genes
                total_number_of_genes = len( results['SEQUENCE'].unique() )
                bacannot_summary[sample]['virulence']['VFDB']['total'] = total_number_of_genes

                # gene annotations
                for gene in [ str(x) for x in results['SEQUENCE'].unique() ]:
                    row = results.loc[results['SEQUENCE'] == gene]
                    bacannot_summary[sample]['virulence']['VFDB'][gene] = {}
                    bacannot_summary[sample]['virulence']['VFDB'][gene]['virulence_factor'] = row['PRODUCT'].item().split('_(')[0].replace("[", "")
                    bacannot_summary[sample]['virulence']['VFDB'][gene]['id'] = row['PRODUCT'].item().split('_(')[1].split(')')[0].replace(")", "")
                    bacannot_summary[sample]['virulence']['VFDB'][gene]['name'] = row['GENE'].item().replace(")", "").replace("(", "")
            
            # victors
            if os.path.exists(f"{results_dir}/virulence/victors/{sample}_victors_blastp_onGenes.summary.txt"):

                # init victors annotation dictionary
                bacannot_summary[sample]['virulence']['Victors'] = {}

                # load victors results
                results = pd.read_csv(
                    f"{results_dir}/virulence/victors/{sample}_victors_blastp_onGenes.summary.txt",
                    sep='\t'
                )

                # number of virulence genes
                total_number_of_genes = len( results['SEQUENCE'].unique() )
                bacannot_summary[sample]['virulence']['Victors']['total'] = total_number_of_genes

                # gene annotations
                for gene in [ str(x) for x in results['SEQUENCE'].unique() ]:
                    row = results.loc[results['SEQUENCE'] == gene]
                    bacannot_summary[sample]['virulence']['Victors'][gene] = {}
                    bacannot_summary[sample]['virulence']['Victors'][gene]['id'] = row['VICTORS_ID'].item().replace("Victors_", "")
                    bacannot_summary[sample]['virulence']['Victors'][gene]['name'] = row['GENE'].item()