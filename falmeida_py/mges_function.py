##################################
### Loading Necessary Packages ###
##################################
import pandas as pd
import os
from .utils import load_and_subset_gff

####################################
### check MGEs annotations stats ###
####################################
def mges_stats(bacannot_summary):
    
    # iterate over available samples
    for sample in bacannot_summary:

        # load dir of samples' results
        results_dir = bacannot_summary[sample]['results_dir']

        # load gff_file
        gff_file = f"{results_dir}/gffs/{sample}.gff"
            
        # integron_finder
        if os.path.exists(f"{results_dir}/integron_finder/{sample}_integrons.gff") and os.stat(f"{results_dir}/integron_finder/{sample}_integrons.gff").st_size > 0:

            # init MGE annotation dictionary
            if 'MGE' not in bacannot_summary[sample]:
                bacannot_summary[sample]['MGE'] = {}

            # init integron_finder annotation dictionary
            bacannot_summary[sample]['MGE']['integron_finder'] = {}

            # load integron_finder results
            results = pd.read_csv(
                f"{results_dir}/integron_finder/{sample}_integrons.gff",
                sep='\t',
                header=None,
                names=[
                    'chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'atts'
                ]
            )

            # number of integron_finder annotations
            total_number = len(results.index)
            bacannot_summary[sample]['MGE']['integron_finder']['total'] = total_number

            # per integron info
            bacannot_summary[sample]['MGE']['integron_finder'] = {}
            if int(results.shape[0]) > 0:
                for seq in [ str(x) for x in results['chr'].unique() ]:
                    
                    bacannot_summary[sample]['MGE']['integron_finder'][seq] = {}
                    for index, row in results[results['chr'] == seq].reset_index().iterrows():
                        id = row['atts'].split(';')[0].split('=')[-1]
                        bacannot_summary[sample]['MGE']['integron_finder'][seq][id] = {}
                        bacannot_summary[sample]['MGE']['integron_finder'][seq][id]['id'] = id
                        bacannot_summary[sample]['MGE']['integron_finder'][seq][id]['contig'] = row['chr']
                        bacannot_summary[sample]['MGE']['integron_finder'][seq][id]['start'] = row['start']
                        bacannot_summary[sample]['MGE']['integron_finder'][seq][id]['end'] = row['end']
                        bacannot_summary[sample]['MGE']['integron_finder'][seq][id]['type'] = row['atts'].split(';')[1].split('=')[-1]
                        bacannot_summary[sample]['MGE']['integron_finder'][seq][id]['source'] = row['source']
                        bacannot_summary[sample]['MGE']['integron_finder'][seq][id]['product'] = row['type']

        # ICE database
        ice_db_blastp = f"{results_dir}/ICEs/{sample}_iceberg_blastp_onGenes.summary.txt"
        if os.path.exists(ice_db_blastp) and os.stat(ice_db_blastp).st_size > 0:

            # init MGE annotation dictionary
            if 'MGE' not in bacannot_summary[sample]:
                bacannot_summary[sample]['MGE'] = {}
            
            # init iceberg annotation dictionary
            if 'ICE' not in bacannot_summary[sample]['MGE']:
                bacannot_summary[sample]['MGE']['ICEberg'] = {}

            # init iceberg blastp annotation dictionary
            bacannot_summary[sample]['MGE']['ICEberg']['blastp'] = {}

            # load integron_finder results
            results = pd.read_csv(
                ice_db_blastp,
                sep='\t'
            )

            # load gff
            gff = load_and_subset_gff(gff_file, 'source', 'ICE')

            # number of integron_finder annotations
            total_number = len(results.index)
            bacannot_summary[sample]['MGE']['ICEberg']['blastp']['total'] = total_number

            # per gene info
            if int(results.shape[0]) > 0:
                for seq in [ str(x) for x in results['SEQUENCE'].unique() ]:

                    # details missing in output but available in gff
                    gff_row = gff[gff['attributes'].str.contains(seq)]
                    contig = gff_row['seq'].item()
                    start  = gff_row['start'].item()
                    end    = gff_row['end'].item()
                    
                    bacannot_summary[sample]['MGE']['ICEberg']['blastp'][seq] = {}
                    for index, row in results[results['SEQUENCE'] == seq].reset_index().iterrows():
                        bacannot_summary[sample]['MGE']['ICEberg']['blastp'][seq][id] = {}
                        bacannot_summary[sample]['MGE']['ICEberg']['blastp'][seq][id]['id']          = row['ICEBERG_ID']
                        bacannot_summary[sample]['MGE']['ICEberg']['blastp'][seq][id]['contig']      = contig
                        bacannot_summary[sample]['MGE']['ICEberg']['blastp'][seq][id]['start']       = start
                        bacannot_summary[sample]['MGE']['ICEberg']['blastp'][seq][id]['end']         = end
                        bacannot_summary[sample]['MGE']['ICEberg']['blastp'][seq][id]['accession']   = row['ACCESSION']
                        bacannot_summary[sample]['MGE']['ICEberg']['blastp'][seq][id]['product']     = row['PRODUCT']
                        bacannot_summary[sample]['MGE']['ICEberg']['blastp'][seq][id]['description'] = row['DESCRIPTION']
                        bacannot_summary[sample]['MGE']['ICEberg']['blastp'][seq][id]['blast_start'] = row['START']
                        bacannot_summary[sample]['MGE']['ICEberg']['blastp'][seq][id]['blast_end']   = row['END']
                        bacannot_summary[sample]['MGE']['ICEberg']['blastp'][seq][id]['strand']      = row['STRAND']