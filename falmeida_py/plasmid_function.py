##################################
### Loading Necessary Packages ###
##################################
import pandas as pd
import os

########################################
### check plasmids annotations stats ###
########################################
def plasmids_stats(bacannot_summary):
    
    # iterate over available samples
    for sample in bacannot_summary:

        # load dir of samples' results
        results_dir = bacannot_summary[sample]['results_dir']

        # init plasmids annotation dictionary
        bacannot_summary[sample]['plasmid'] = {}            
            
        # platon
        if os.path.exists(f"{results_dir}/plasmids/platon/{sample}.tsv") and os.stat(f"{results_dir}/plasmids/platon/{sample}.tsv").st_size > 0:

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
            if int(results.shape[0]) > 0:
                for seq in [x for x in results['ID'].unique()]:
                    bacannot_summary[sample]['plasmid']['platon'][seq] = {}
                    bacannot_summary[sample]['plasmid']['platon'][seq]['Length'] = results.loc[results['ID'] == seq, 'Length'].item()
                    bacannot_summary[sample]['plasmid']['platon'][seq]['ORFs'] = results.loc[results['ID'] == seq, '# ORFs'].item()
                    bacannot_summary[sample]['plasmid']['platon'][seq]['Circular'] = results.loc[results['ID'] == seq, 'Circular'].item()
                    bacannot_summary[sample]['plasmid']['platon'][seq]['AMRs'] = results.loc[results['ID'] == seq, '# AMRs'].item()
                    bacannot_summary[sample]['plasmid']['platon'][seq]['Replication'] = results.loc[results['ID'] == seq, '# Replication'].item()
                    bacannot_summary[sample]['plasmid']['platon'][seq]['Mobilization'] = results.loc[results['ID'] == seq, '# Mobilization'].item()
                    bacannot_summary[sample]['plasmid']['platon'][seq]['Conjugation'] = results.loc[results['ID'] == seq, '# Conjugation'].item()
        
        # plasmidfinder
        if os.path.exists(f"{results_dir}/plasmids/plasmidfinder/results_tab.tsv") and os.stat(f"{results_dir}/plasmids/plasmidfinder/results_tab.tsv").st_size > 0:

            # init platon annotation dictionary
            bacannot_summary[sample]['plasmid']['plasmidfinder'] = {}

            # load platon results
            results = pd.read_csv(
                f"{results_dir}/plasmids/plasmidfinder/results_tab.tsv",
                sep='\t'
            )

            if not results.empty:

                # databases
                print( results['Database'].unique() )
                bacannot_summary[sample]['plasmid']['plasmidfinder']['meta'] = {}
                db_arr = results['Database'].unique()
                bacannot_summary[sample]['plasmid']['plasmidfinder']['meta']['database'] = db_arr.tolist() if len(db_arr) > 1 else db_arr.item()

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
                    contig = str(row['Contig'])
                    bacannot_summary[sample]['plasmid']['plasmidfinder'][contig]['inc_types'] = row['Plasmid']
                    bacannot_summary[sample]['plasmid']['plasmidfinder'][contig]['identity'] = row['Identity']
                    bacannot_summary[sample]['plasmid']['plasmidfinder'][contig]['accession'] = row['Accession number']
        
        # mob suite
        if os.path.exists(f"{results_dir}/plasmids/mob_suite/{sample}_mobtyper_results.txt") and os.stat(f"{results_dir}/plasmids/mob_suite/{sample}_mobtyper_results.txt").st_size > 0:

            # init integron_finder annotation dictionary
            bacannot_summary[sample]['plasmid']['mob_suite'] = {}

            # load integron_finder results
            results = pd.read_csv(
                f"{results_dir}/plasmids/mob_suite/{sample}_mobtyper_results.txt",
                sep='\t',
                header='infer',
                # sample_id	num_contigs	size	gc	md5	rep_type(s)	rep_type_accession(s)	relaxase_type(s)	relaxase_type_accession(s)	mpf_type	mpf_type_accession(s)	orit_type(s)	orit_accession(s)	predicted_mobility	mash_nearest_neighbor	mash_neighbor_distance	mash_neighbor_identification	primary_cluster_id	secondary_cluster_id	predicted_host_range_overall_rank	predicted_host_range_overall_name	observed_host_range_ncbi_rank	observed_host_range_ncbi_name	reported_host_range_lit_rank	reported_host_range_lit_name	associated_pmid(s)
            )

            # number of plasmid types annotated annotations
            # total_number = len(results.index) - 1 # always counts chromosome
            # bacannot_summary[sample]['plasmid']['mob_suite']['total'] = total_number

            # per integron info
            bacannot_summary[sample]['plasmid']['mob_suite'] = {}
            if int(results.shape[0]) > 0:
                for seq in [ str(x) for x in results['sample_id'].unique() ]:
                    
                    bacannot_summary[sample]['plasmid']['mob_suite'][seq] = {}
                    for index, row in results[results['sample_id'].astype(str) == seq].reset_index().iterrows():
                        id = row['sample_id'] # they are the same for this result
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id] = {}
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['size'] = row['size']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['rep_type'] = row['rep_type(s)'].replace(',', '; ')
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['rep_type_accession'] = row['rep_type_accession(s)'].replace(',', '; ')
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['relaxase_type'] = row['relaxase_type(s)']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['relaxase_type_accession(s)'] = row['relaxase_type_accession(s)']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['mpf_type'] = row['mpf_type']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['mpf_type_accession'] = row['mpf_type_accession(s)']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['orit_type'] = row['orit_type(s)']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['orit_accession'] = row['orit_accession(s)']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['mash_nearest_neighbor'] = row['mash_nearest_neighbor']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['mash_neighbor_distance'] = row['mash_neighbor_distance']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['mash_neighbor_identification'] = row['mash_neighbor_identification']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['primary_cluster_id'] = row['primary_cluster_id']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['secondary_cluster_id'] = row['secondary_cluster_id']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['predicted_host_range_overall_rank'] = row['predicted_host_range_overall_rank']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['predicted_host_range_overall_name'] = row['predicted_host_range_overall_name']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['observed_host_range_ncbi_rank'] = row['observed_host_range_ncbi_rank']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['observed_host_range_ncbi_name'] = row['observed_host_range_ncbi_name']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['reported_host_range_lit_rank'] = row['reported_host_range_lit_rank']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['reported_host_range_lit_name'] = row['reported_host_range_lit_name']
                        bacannot_summary[sample]['plasmid']['mob_suite'][seq][id]['associated_pmid'] = row['associated_pmid(s)']