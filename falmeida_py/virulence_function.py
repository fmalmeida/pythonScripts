##################################
### Loading Necessary Packages ###
##################################
import pandas as pd
from .utils import find_files
import json
import os
import yaml
from pathlib import Path
from .utils import load_and_subset_gff

##################################
### check virulence annotation ###
##################################
def virulence_stats(bacannot_summary):

    # iterate over available samples
    for sample in bacannot_summary:

        # load dir of samples' results
        results_dir = bacannot_summary[sample]['results_dir']

        # load gff_file
        gff_file = f"{results_dir}/gffs/{sample}.gff"

        # load annotation stats
        if os.path.exists(f"{results_dir}/virulence"):

            # init virulence annotation dictionary
            bacannot_summary[sample]['virulence'] = {}
            
            # vfdb
            if os.path.exists(f"{results_dir}/virulence/vfdb/{sample}_vfdb_blastn_onGenes.summary.txt") and os.stat(f"{results_dir}/virulence/vfdb/{sample}_vfdb_blastn_onGenes.summary.txt").st_size > 0:

                # init VFDB annotation dictionary
                bacannot_summary[sample]['virulence']['VFDB'] = {}

                # load gff
                gff = load_and_subset_gff(gff_file, 'source', 'VFDB')

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
                    
                    # init values
                    bacannot_summary[sample]['virulence']['VFDB'][gene] = {}
                    row = results.loc[results['SEQUENCE'] == gene]
                    vf_name = row['PRODUCT'].item().split('_(')[0].replace("[", "")
                    vf_id = row['PRODUCT'].item().split('_(')[1].split(')')[0].replace(")", "")
                    vf_fullname = row['PRODUCT'].item().replace("[", "").replace("]", "")
                    gene_name = row['GENE'].item().replace(")", "").replace("(", "")
                    gff_row = gff[gff['attributes'].str.contains(gene)]
                    contig = gff_row['seq'].item()
                    start  = gff_row['start'].item()
                    end    = gff_row['end'].item()

                    # add to dict
                    bacannot_summary[sample]['virulence']['VFDB'][gene]['virulence_factor'] = vf_name
                    bacannot_summary[sample]['virulence']['VFDB'][gene]['product'] = vf_fullname
                    bacannot_summary[sample]['virulence']['VFDB'][gene]['id'] = vf_id
                    bacannot_summary[sample]['virulence']['VFDB'][gene]['gene'] = gene_name
                    bacannot_summary[sample]['virulence']['VFDB'][gene]['chr'] = contig
                    bacannot_summary[sample]['virulence']['VFDB'][gene]['start'] = start
                    bacannot_summary[sample]['virulence']['VFDB'][gene]['end'] = end
            
            # victors
            if os.path.exists(f"{results_dir}/virulence/victors/{sample}_victors_blastp_onGenes.summary.txt") and os.stat(f"{results_dir}/virulence/victors/{sample}_victors_blastp_onGenes.summary.txt").st_size > 0:

                # init victors annotation dictionary
                gff = bacannot_summary[sample]['virulence']['Victors'] = {}

                # load victors results
                results = pd.read_csv(
                    f"{results_dir}/virulence/victors/{sample}_victors_blastp_onGenes.summary.txt",
                    sep='\t'
                )

                # load gff
                gff = load_and_subset_gff(gff_file, 'source', 'Victors')

                # number of virulence genes
                total_number_of_genes = len( results['SEQUENCE'].unique() )
                bacannot_summary[sample]['virulence']['Victors']['total'] = total_number_of_genes

                # gene annotations
                for gene in [ str(x) for x in results['SEQUENCE'].unique() ]:

                    # init values
                    row = results.loc[results['SEQUENCE'] == gene]
                    bacannot_summary[sample]['virulence']['Victors'][gene] = {}
                    vf_id = row['VICTORS_ID'].item().replace("Victors_", "")
                    gene_name = row['GENE'].item()
                    product = row['DESCRIPTION'].item()
                    gff_row = gff[gff['attributes'].str.contains(gene)]
                    contig = gff_row['seq'].item()
                    start  = gff_row['start'].item()
                    end    = gff_row['end'].item()

                    # add values to dict
                    bacannot_summary[sample]['virulence']['Victors'][gene]['id'] = vf_id
                    bacannot_summary[sample]['virulence']['Victors'][gene]['name'] = gene_name
                    bacannot_summary[sample]['virulence']['Victors'][gene]['product'] = product
                    bacannot_summary[sample]['virulence']['Victors'][gene]['contig'] = contig
                    bacannot_summary[sample]['virulence']['Victors'][gene]['start'] = start
                    bacannot_summary[sample]['virulence']['Victors'][gene]['end'] = end