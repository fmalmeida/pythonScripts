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

##########################################
### check resistance annotations stats ###
##########################################
def resistance_stats(bacannot_summary):
    
    # iterate over available samples
    for sample in bacannot_summary:

        # load dir of samples' results
        results_dir = bacannot_summary[sample]['results_dir']

        # load gff_file
        gff_file = f"{results_dir}/gffs/{sample}.gff"

        # load annotation stats
        if os.path.exists(f"{results_dir}/resistance"):

            # init plasmids annotation dictionary
            bacannot_summary[sample]['resistance'] = {}
            
            #####################
            ### amrfinderplus ###
            #####################
            if os.path.exists(f"{results_dir}/resistance/AMRFinderPlus/AMRFinder_resistance-only.tsv"):

                # init amrfinderplus annotation dictionary
                bacannot_summary[sample]['resistance']['amrfinderplus'] = {}

                # load amrfinderplus results
                results = pd.read_csv(
                    f"{results_dir}/resistance/AMRFinderPlus/AMRFinder_resistance-only.tsv",
                    sep='\t'
                )

                # load gff
                gff = load_and_subset_gff(gff_file, 'source', 'AMRFinderPlus')

                # number of annotations
                total_number = len(results['Protein identifier'].unique())
                bacannot_summary[sample]['resistance']['amrfinderplus']['total'] = total_number

                # gene annotations
                for gene in [ str(x) for x in results['Protein identifier'].unique() ]:

                    # init values
                    row = results.loc[results['Protein identifier'] == gene]
                    bacannot_summary[sample]['resistance']['amrfinderplus'][gene] = {}
                    gene_name = row['Gene symbol'].item()
                    drug_class = row['Subclass'].item()
                    identity = row['% Identity to reference sequence'].item()

                    gff_row = gff[gff['attributes'].str.contains(f"ID={gene}")]
                    contig = gff_row['seq'].item()
                    start  = gff_row['start'].item()
                    end    = gff_row['end'].item()

                    # add values to dict
                    bacannot_summary[sample]['resistance']['amrfinderplus'][gene]['gene'] = gene_name
                    bacannot_summary[sample]['resistance']['amrfinderplus'][gene]['subclass'] = drug_class
                    bacannot_summary[sample]['resistance']['amrfinderplus'][gene]['identity'] = identity
                    bacannot_summary[sample]['resistance']['amrfinderplus'][gene]['contig'] = contig
                    bacannot_summary[sample]['resistance']['amrfinderplus'][gene]['start'] = start
                    bacannot_summary[sample]['resistance']['amrfinderplus'][gene]['end'] = end
            
            #################
            ### resfinder ###
            #################
            #
            # TODO: Include genomic coordinates info
            #
            if os.path.exists(f"{results_dir}/resistance/resfinder/ResFinder_results_tab.txt"):

                # init resfinder annotation dictionary
                bacannot_summary[sample]['resistance']['resfinder'] = {}

                # load gff
                gff = load_and_subset_gff(gff_file, 'source', 'Resfinder')

                # number of annotations
                bacannot_summary[sample]['resistance']['resfinder']['total'] = len(gff.index)

                # since resfinder output does not has locus_tag information
                # for this module, we will be using the resolved gff from bacannot
                for index, row in gff.iterrows():

                    # init attributes as dict
                    attributes = dict()
                    for keyvaluepair in row['attributes'].split(';'):
                        items = keyvaluepair.split('=')
                        key   = items[0]
                        value = '='.join(items[1:])
                        attributes[key] = value
                    gene = attributes['ID']

                    # parse
                    bacannot_summary[sample]['resistance']['resfinder'][gene] = {}
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['start'] = row['start']
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['end']   = row['end']
                    # bacannot_summary[sample]['resistance']['resfinder'][gene]['Identity'] = row['Identity']
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['name'] = attributes['Resfinder_gene']
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['phenotype'] = attributes['Resfinder_phenotype']
                    bacannot_summary[sample]['resistance']['resfinder'][gene]['accession'] = attributes['Resfinder_reference']
            
            ###########
            ### rgi ###
            ###########
            if os.path.exists(f"{results_dir}/resistance/RGI/RGI_{sample}.txt"):

                # init amrfinderplus annotation dictionary
                bacannot_summary[sample]['resistance']['rgi'] = {}

                # load amrfinderplus results
                results = pd.read_csv(
                    f"{results_dir}/resistance/RGI/RGI_{sample}.txt",
                    sep='\t'
                )
                results.drop_duplicates(inplace=True)

                # load gff
                gff = load_and_subset_gff(gff_file, 'source', 'CARD')

                # number of annotations
                total_number = len(results['ORF_ID'].unique())
                bacannot_summary[sample]['resistance']['rgi']['total'] = total_number

                # gene annotations
                for gene in [ str(x) for x in results['ORF_ID'].unique() ]:

                    # init values
                    row = results.loc[results['ORF_ID'] == gene]
                    name_split_list = gene.split(' ')
                    gene = name_split_list[0]
                    name = ' '.join(name_split_list[1:])
                    bacannot_summary[sample]['resistance']['rgi'][gene] = {}
                    gff_row = gff[gff['attributes'].str.contains(f"ID={gene}")]
                    contig = gff_row['seq'].item()
                    start  = gff_row['start'].item()
                    end    = gff_row['end'].item()

                    # add values to dict
                    
                    
                    
                    bacannot_summary[sample]['resistance']['rgi'][gene]['name'] = name
                    bacannot_summary[sample]['resistance']['rgi'][gene]['gene'] = row['Best_Hit_ARO'].item()
                    bacannot_summary[sample]['resistance']['rgi'][gene]['cut_off'] = row['Cut_Off'].item()
                    bacannot_summary[sample]['resistance']['rgi'][gene]['resistance_mechanism'] = row['Resistance Mechanism'].item()
                    bacannot_summary[sample]['resistance']['rgi'][gene]['gene_family'] = row['AMR Gene Family'].item()
                    bacannot_summary[sample]['resistance']['rgi'][gene]['subclass'] = row['Drug Class'].item()
                    bacannot_summary[sample]['resistance']['rgi'][gene]['identity'] = row['Best_Identities'].item()
                    bacannot_summary[sample]['resistance']['rgi'][gene]['accession'] = row['Model_ID'].item()
                    bacannot_summary[sample]['resistance']['rgi'][gene]['contig'] = contig
                    bacannot_summary[sample]['resistance']['rgi'][gene]['start'] = start
                    bacannot_summary[sample]['resistance']['rgi'][gene]['end'] = end