##################################
### Loading Necessary Packages ###
##################################
import pandas as pd
import os

####################################
### check MGEs annotations stats ###
####################################
def mges_stats(bacannot_summary):
    
    # iterate over available samples
    for sample in bacannot_summary:

        # load dir of samples' results
        results_dir = bacannot_summary[sample]['results_dir']

        # init MGE annotation dictionary
        bacannot_summary[sample]['MGE'] = {}
            
        # integron_finder
        if os.path.exists(f"{results_dir}/integron_finder/{sample}_integrons.gff"):

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