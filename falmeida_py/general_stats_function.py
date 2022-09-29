##################################
### Loading Necessary Packages ###
##################################
import pandas as pd
import yaml
from pathlib import Path

################################
### check general annotation ###
################################
def general_stats(bacannot_summary):

    # iterate over available samples
    for sample in bacannot_summary:

        # load dir of samples' results
        results_dir = bacannot_summary[sample]['results_dir']

        # load annotation stats
        general_results = yaml.safe_load(
            Path(f"{results_dir}/annotation/{sample}.txt").read_text()
        )

        # load MLST
        mlst_results = pd.read_csv(
            f"{results_dir}/MLST/{sample}_mlst_analysis.txt",
            sep='\t', header=None
        )

        # load refseq_masher
        refseq_masher_results = pd.read_csv(
            f"{results_dir}/refseq_masher/refseq_masher_results.txt",
            sep='\t'
        )
        refseq_masher_results.sort_values(by='distance', ascending=True, inplace=True)

        # save annotation stats
        bacannot_summary[sample]['general_annotation'] = {}
        bacannot_summary[sample]['general_annotation']['mlst']  = str(mlst_results[2].item()).replace('-', 'null')
        bacannot_summary[sample]['general_annotation']['cds']   = general_results['CDS']
        bacannot_summary[sample]['general_annotation']['rrna']  = general_results['rRNA']
        bacannot_summary[sample]['general_annotation']['trna']  = general_results['tRNA']
        bacannot_summary[sample]['general_annotation']['tmrna'] = general_results['tmRNA']

        bacannot_summary[sample]['general_annotation']['closest_reference'] = {}
        bacannot_summary[sample]['general_annotation']['closest_reference']['strain'] = refseq_masher_results.head(1)['top_taxonomy_name'].item()
        bacannot_summary[sample]['general_annotation']['closest_reference']['distance'] = refseq_masher_results.head(1)['distance'].item()
        bacannot_summary[sample]['general_annotation']['closest_reference']['accession'] = refseq_masher_results.head(1)['assembly_accession'].item()