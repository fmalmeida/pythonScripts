##################################
### Loading Necessary Packages ###
##################################
import pandas as pd
from .utils import find_files
import json
import os
import yaml
from pathlib import Path

##########################################
### check resistance annotations stats ###
##########################################
def resistance_stats(bacannot_summary):
    
    # iterate over available samples
    for sample in bacannot_summary:

        # load dir of samples' results
        results_dir = bacannot_summary[sample]['results_dir']

        # load annotation stats
        if os.path.exists(f"{results_dir}/resistance"):

            # init plasmids annotation dictionary
            bacannot_summary[sample]['resistance'] = {}
            
            # amrfinderplus
            if os.path.exists(f"{results_dir}/resistance/AMRFinderPlus/AMRFinder_resistance-only.tsv"):

                # init amrfinderplus annotation dictionary
                bacannot_summary[sample]['resistance']['amrfinderplus'] = {}

                # load amrfinderplus results
                results = pd.read_csv(
                    f"{results_dir}/resistance/AMRFinderPlus/AMRFinder_resistance-only.tsv",
                    sep='\t'
                )
                results.sort_values('Gene symbol', inplace=True)

                # number of annotations
                total_number = len(results['Protein identifier'].unique())
                bacannot_summary[sample]['resistance']['amrfinderplus']['total'] = total_number

                # gene annotations
                for gene in [ str(x) for x in results['Protein identifier'].unique() ]:
                    row = results.loc[results['Protein identifier'] == gene]
                    bacannot_summary[sample]['resistance']['amrfinderplus'][gene] = {}
                    bacannot_summary[sample]['resistance']['amrfinderplus'][gene]['gene'] = row['Gene symbol'].item()
                    bacannot_summary[sample]['resistance']['amrfinderplus'][gene]['subclass'] = row['Subclass'].item()
                    bacannot_summary[sample]['resistance']['amrfinderplus'][gene]['identity'] = row['% Identity to reference sequence'].item()
            

            # resfinder
            if os.path.exists(f"{results_dir}/resistance/resfinder/ResFinder_results_tab.txt"):

                # init resfinder annotation dictionary
                bacannot_summary[sample]['resistance']['resfinder'] = {}

                # load resfinder results
                results = pd.read_csv(
                    f"{results_dir}/resistance/resfinder/ResFinder_results_tab.txt",
                    sep='\t'
                )
                results.drop_duplicates(inplace=True)

                # number of annotations
                bacannot_summary[sample]['resistance']['amrfinderplus']['total'] = len(results.index)

                # gene annotations
                for gene in [ str(x) for x in results['Resistance gene'].unique() ]:

                    # init
                    bacannot_summary[sample]['resistance']['resfinder'][gene] = {}
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['chr']   = list()
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['start'] = list()
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['end']   = list()
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['Identity']    = list()
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['phenotype']   = list()
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['accession']   = list()

                # parse
                for index, row in results.iterrows():

                    gene = row['Resistance gene']

                    bacannot_summary[sample]['resistance']['resfinder'][gene]['chr'].append(
                        row['Contig']
                    )
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['start'].append(
                        row['Position in contig'].split('..')[0]
                    )
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['end'].append(
                        row['Position in contig'].split('..')[1]
                    )
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['Identity'].append(
                        row['Identity']
                    )
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['phenotype'].append(
                        row['Phenotype']
                    )
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['accession'].append(
                        row['Accession no.']
                    )
            
            # rgi
            if os.path.exists(f"{results_dir}/resistance/RGI/RGI_{sample}.txt"):

                # init amrfinderplus annotation dictionary
                bacannot_summary[sample]['resistance']['rgi'] = {}

                # load amrfinderplus results
                results = pd.read_csv(
                    f"{results_dir}/resistance/RGI/RGI_{sample}.txt",
                    sep='\t'
                )
                results.drop_duplicates(inplace=True)

                # number of annotations
                total_number = len(results['ORF_ID'].unique())
                bacannot_summary[sample]['resistance']['rgi']['total'] = total_number

                # gene annotations
                for gene in [ str(x) for x in results['ORF_ID'].unique() ]:
                    row = results.loc[results['ORF_ID'] == gene]
                    name_split_list = gene.split(' ')
                    gene = name_split_list[0]
                    name = ' '.join(name_split_list[1:])
                    bacannot_summary[sample]['resistance']['rgi'][gene] = {}
                    bacannot_summary[sample]['resistance']['rgi'][gene]['name'] = name
                    bacannot_summary[sample]['resistance']['rgi'][gene]['gene'] = row['Best_Hit_ARO'].item()
                    bacannot_summary[sample]['resistance']['rgi'][gene]['cut_off'] = row['Cut_Off'].item()
                    bacannot_summary[sample]['resistance']['rgi'][gene]['resistance_mechanism'] = row['Resistance Mechanism'].item()
                    bacannot_summary[sample]['resistance']['rgi'][gene]['gene_family'] = row['AMR Gene Family'].item()
                    bacannot_summary[sample]['resistance']['rgi'][gene]['subclass'] = row['Drug Class'].item()
                    bacannot_summary[sample]['resistance']['rgi'][gene]['identity'] = row['Best_Identities'].item()
                    bacannot_summary[sample]['resistance']['rgi'][gene]['accession'] = row['Model_ID'].item()