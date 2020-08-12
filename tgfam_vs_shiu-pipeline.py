#!/usr/bin/env python3
# coding: utf-8

## Def help message
"""
This script was created in order to easily and automatically identify overlapping
regions (sequences) between TGFam-Finder and Shiu's Pseudogene pipeline annotations.
This pipeline is meant to standardize the annotation between all target gene families.

---
Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain
Usage:
    tgfam_vs_shiu-pipeline.py.py
    tgfam_vs_shiu-pipeline.py.py -h|--help
    tgfam_vs_shiu-pipeline.py.py -v|--version
    tgfam_vs_shiu-pipeline.py.py summary [--tgfam <dir> --shiu <gff> --mutations <int> --bedtools <path>]
    tgfam_vs_shiu-pipeline.py.py compare [--tgfam <dir> --shiu <gff> --mutations <int> --bedtools <path> --tgfam_fraction <int> --outdir <dir>]

Options:
    -h --help                   Show this screen.
    -v --version                Show version information
    --bedtools=<path>           Path to bedtools executable, v2.29.2 is suggested [default: bedtools]
    --tgfam_fraction=<int>      Minimum overlap required as a fraction of TGFam genes in order to report a intersection [default: 0.1]
    --shiu=<gff>                GFF file containing Pseudogenes predicted with Shiu's pipeline [default: ./gffs/Pseudogenes/bnut.shiu.Pseudogenes.gff]
    --mutations=<int>           The total of disabling mutations to accept while searching for "relaxed Pseudogenes" [default: 3]
    --tgfam=<dir>               Directory containing TGFam-Finder GFF files [default: ./gffs/TGFam]
    --outdir=<dir>              Path to output directory [default: ./_results]
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt
import pandas as pd
import re
import os
import subprocess

###########################################
### Function for filtering python lists ###
###########################################
def filter(string, substr):
    return [str for str in string if
             any(sub in str for sub in substr)]

###############################################################
### Function for reading bedtools intersection as dataframe ###
###############################################################
def intersectBed_2_DF(inputFile):
    df = pd.read_csv(inputFile, sep = "\t", comment = "#",
                        names=['Shiu.Chr', 'Shiu.Source', 'Shiu.Type', 'Shiu.Start', 'Shiu.End', 'Shiu.Score', 'Shiu.Strand', 'Shiu.Phase', 'Shiu.Attributes',
                               'TGFAM.Chr', 'TGFAM.Source', 'TGFAM.Type', 'TGFAM.Start', 'TGFAM.End', 'TGFAM.Score', 'TGFAM.Strand', 'TGFAM.Phase', 'TGFAM.Attributes', 'Overlap'])
    return df

############################################################
### Function for running bedtools intersect between GFFs ###
############################################################
def intersectBed(shiu, tgfam, bedtools, overlap_fraction, out):
    os.system(f"{bedtools} intersect -wo -a {shiu} -b {tgfam} -F {overlap_fraction} | awk '{{ if ($12 == \"gene\") print }}' > {out}")

##############################################
### Filter a GFF -- Based on a list of IDS ###
##############################################

def filter_gff(input_gff, entry_list, output, mode="normal"):
    # mode = reverse sets the function to reversely match a list of values

    # Read GFF file
    gff_df = pd.read_csv(input_gff, sep = "\t", comment = "#",
                        names=['Chr', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'])

    # Parse -- It works simple as this because all the lines from both GFFs
    # contain identical IDs, even though they are parent and children features
    # Different of JGI annotation which uses a different ID for parents and children
    for index, line in gff_df.iterrows():

        # Grab ID
        id = filter(line['Attributes'].split(';'), ['ID'])[0].split('=')[1]

        # Select desired entries
        with open(f"{output}", 'a+') as gff:

            # Check if user wants to reversely match a list
            if mode == "reverse":
                if id not in entry_list:
                    print(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8],
                          sep="\t", file=gff)
            else:
                if id in entry_list:
                    print(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8],
                          sep="\t", file=gff)

##################################
### GFF summarization pipeline ###
##################################
def summary_gffs(shiuGFF, tgfamDIR, bedtools, overlap_fraction, relaxation):
    # Remove result if existent
    if os.path.exists("Intersection_summary.md"):
        os.remove("Intersection_summary.md")

    # Recreate  summary
    for filename in os.listdir(tgfamDIR):
        if filename.endswith(".gff") or filename.endswith(".gff3"):

            # Get current family name
            gene_family = filename.split(".")[0]
            tgfamGFF = os.path.join(tgfamDIR, filename)

            # Check overlaps - Only in gene features (TGFam)
            intersectBed(shiu=shiuGFF, tgfam=tgfamGFF, bedtools=bedtools, out="tmp.intersect", overlap_fraction=f"{overlap_fraction}")

            # Parse intersections
            intersect_df = intersectBed_2_DF("tmp.intersect")

            # Lists and start
            overlaps_ids      = [] # List to store TGFam all gene ids with intersection
            relaxed_genes_ids = [] # List to store TGFam all relaxed genes ids with intersection
            intact_genes_ids  = [] # List to store TGFam all intact genes ids with intersection

            # Loop intersections
            for index, line in intersect_df.iterrows():
                tgfam_id = filter(line["TGFAM.Attributes"].split(';'), ['ID'])[0].split("=")[1] # Select TGFam gene id from att column
                shiuNotes     = filter(line['Shiu.Attributes'].split(';'), ['Note'])[0].split('=')[1] # Select Shiu's notes in GFFs
                notesEvidence = shiuNotes.split("_")[3].split(",") # Disabling mutations given as evidence
                disabling_mutations = [int(i) for i in notesEvidence] # Convert to integer
                total = sum(disabling_mutations) # Sum the number of disabling mutations

                # The following construction is used to avoid counting a TGFam gene with more than one intersection
                # multiple times. It checks if the gene ID is already present in the list of gene ids before appending.
                # The final list will be used to count the amount of (unique) TGFam genes with intersections.
                if total == 0:
                    if tgfam_id not in overlaps_ids:
                        overlaps_ids.append(tgfam_id)
                    if tgfam_id not in intact_genes_ids:
                        intact_genes_ids.append(tgfam_id)

                elif total >= int(relaxation):
                    if tgfam_id not in overlaps_ids:
                        overlaps_ids.append(tgfam_id)
                    if tgfam_id not in relaxed_genes_ids:
                        relaxed_genes_ids.append(tgfam_id)
                else:
                    if tgfam_id not in overlaps_ids:
                        overlaps_ids.append(tgfam_id)

            # Count the number of ids in each list
            count         = len(overlaps_ids)
            intact_genes  = len(intact_genes_ids)
            relaxed_genes = len(relaxed_genes_ids)



            # Done -- Write summary
            with open(f"Intersection_summary.md", 'a+') as f:
                print(f"""
# Summary of {gene_family} (TGFam) gene family

Overview:
    A total of {count} TGFam genes have intersection with one or more pseudogenes annotated with Shiu's pipeline.
        * Using a minimum of {overlap_fraction} overlapping fraction of TGFam genes (bedtools intersect -F parameter).

Pseudogene evidence:

    * Amount of genes (from overlaps) without a single disabling mutation predicted by Shiu's pipeline: {intact_genes}
    * Amount of genes (from overlaps) with 1 to {relaxation} disabling mutations predicted by Shiu's pipeline: {relaxed_genes}
                """, file=f)

            # Remove temporary files
            os.remove("tmp.intersect")

            # print(os.path.join(tgfamDIR, filename))
            continue
        else:
            continue

###############################
### GFF comparison pipeline ###
###############################
def gff_compare(shiuGFF, tgfamDIR, bedtools, outdir, relaxation, overlap_fraction):

    # Remove result if existent
    if os.path.exists(f"{outdir}"):
        os.system(f"rm -r {outdir}")

    # Create output directory
    os.system(f"mkdir -p {outdir}")

    # Create list to store the Pseudogenes that have been
    # promoted to true or relaxed genes in all the gene
    # families annotated
    promoted_Pseudogenes = []

    # Iterate over each TGFam GFF file
    for filename in os.listdir(tgfamDIR):
        if filename.endswith(".gff") or filename.endswith(".gff3"):

            # Get current family name
            gene_family = filename.split(".")[0]
            tgfamGFF = os.path.join(tgfamDIR, filename)

            # Create a results dir for the given gene family
            current_dir = f"{outdir}/{gene_family}"
            os.system(f"mkdir -p {current_dir}")

            # Debug file -- remove if existent
            if os.path.exists(f"./{current_dir}/gff_comparison.log"):
                os.remove(f"./{current_dir}/gff_comparison.log")
            f_debug = open(f"./{current_dir}/gff_comparison.log", "a+")

            # First Step
            # Detect intersections with bedtools
            log = f"""
# STEP 1

Find intersection between TGFam and Pseudogene annotations using bedtools intersect (default)
    + CMD: {bedtools} intersect -wo -a {shiuGFF} -b {tgfamGFF} -F {overlap_fraction} | awk '{{ if ($12 == \"gene\") print }}' > intersected_genes.txt
    + Note: TGFam GFF contains information of genes, mRNAs and CDS. Thus, we used awk to select only overlapping gene features.
    + Bedtools intersect result at: intersected_genes.txt
            """
            f_debug.write(log) # Write to file
            intersectBed(shiu=shiuGFF, tgfam=tgfamGFF, bedtools=bedtools, overlap_fraction=f"{overlap_fraction}", out=f"{current_dir}/intersected_genes.txt")

            # Second Step
            # Separate gene ids by the number of disabling mutations
            log = f"""
# STEP 2

Detecting which of the TGFam genes are true genes, "relaxed" genes or Pseudogenes.
    + Using a threshold of {relaxation} disabling mutations for "relaxed" genes

This is done by parsing the bedtools intersection file over a pandas dataframe and
then evaluating the number of disabling mutations predicted by Shiu's pipeline for
each entry.

TGFam gene ids and Pseudogene ids will be stored at:
    + True genes (thus the promoted Pseudogenes): gene_ids/accepted_as_genes_ids.txt
    + Relaxed genes (thus the Pseudogenes in observation): gene_ids/relaxed_genes_ids.txt
    + True Pseudogenes (thus the TGFam genes that will be relegated): gene_ids/accepted_as_Pseudogenes_ids.txt
            """
            f_debug.write(log) # Write to file

            # Initiate lists
            true_genes        = {} # Dict for storage of TGFam gene ids of "true" genes and consequently of promoted Pseudogenes
            relaxed_genes     = {} # Dict for storage of TGFam gene ids (and Pseudogenes) of genes with few disabling mutations
            true_pseudogenes  = {} # Dict for storage of TGFam gene ids of actually Pseudogenes (and consequently true Pseudogenes)

            # Outputs

            # Parse and loop in the intersections
            intersect_df = intersectBed_2_DF(f"{current_dir}/intersected_genes.txt")
            for index, line in intersect_df.iterrows():

                # Grab TGFam ID and Pseudogene ID
                tgfam_id = filter(line["TGFAM.Attributes"].split(';'), ['ID'])[0].split("=")[1] # Select TGFam gene id from att column
                shiu_id = filter(line["Shiu.Attributes"].split(';'), ['ID'])[0].split("=")[1] # Select Pseudogene id from att column

                # Check Pseudogenes evidences
                shiuNotes     = filter(line['Shiu.Attributes'].split(';'), ['Note'])[0].split('=')[1] # Select Shiu's notes in GFFs
                notesEvidence = shiuNotes.split("_")[3].split(",") # Disabling mutations given as evidence
                disabling_mutations = [int(i) for i in notesEvidence] # Convert to integer
                total = sum(disabling_mutations) # Sum the number of disabling mutations

                if total == 0:
                    true_genes[f"{tgfam_id}"] = f"{shiu_id}"
                elif total >= int(relaxation):
                    relaxed_genes[f"{tgfam_id}"] = f"{shiu_id}"
                else:
                    true_pseudogenes[f"{tgfam_id}"] = f"{shiu_id}"

            # Save ids to lists
            os.system(f"mkdir -p {current_dir}/gene_ids")

            ## true genes
            with open(f"{current_dir}/gene_ids/accepted_as_genes_ids.txt", "w") as outfile:
                print("TGFam ID\tPseudogene ID", file=outfile)
                for key, value in true_genes.items():
                    print(f"{key}\t{value}", file=outfile)

            ## relaxed genes
            with open(f"{current_dir}/gene_ids/relaxed_genes_ids.txt", "w") as outfile:
                print("TGFam ID\tPseudogene ID", file=outfile)
                for key, value in relaxed_genes.items():
                    print(f"{key}\t{value}", file=outfile)

            ## Pseudogenes
            with open(f"{current_dir}/gene_ids/accepted_as_Pseudogenes_ids.txt", "w") as outfile:
                print("TGFam ID\tPseudogene ID", file=outfile)
                for key, value in true_pseudogenes.items():
                    print(f"{key}\t{value}", file=outfile)

            ## Save promoted Pseudogenes in its MAIN list
            for i in [*{**true_genes, **relaxed_genes}.values()]:
                if i not in promoted_Pseudogenes:
                    promoted_Pseudogenes.append(i)

            # Third Step
            # Separate TGFam GFFs per family into: true genes and relaxed genes
            # Relegated genes will be excluded from TGFam results and maintained only in shiu's Pseudogenes results
            log = f"""
# STEP 3

Filtering TGFam results in order to remove "relegated genes", and to create
to separate files for "true genes" and "relaxed genes" of TGFam annotations.
Everything based on the intersection of Pseudogenes and TGFam annotations.

Explanations:

* "True genes"      -- True genes are all the genes from TGFam that we will accept as genes.
                       it contains all the genes that do not had a single base overlap with
                       any of the predicted Pseudogenes in addition to the genes that had any
                       overlap with Pseudogenes without disabling mutations. These former
                       Pseudogenes will be removed from Shiu's final results

* "Relaxed genes"   -- Relaxed genes are all the genes from TGFam that had any overlap with
                       Pseudogenes with 1 to {relaxation} disabling mutations. These are genes
                       that we maintain in observation due to the possibility of problems with
                       the sequencing technology. These former Pseudogenes will be removed from
                       Shiu's final results

* "Relegated genes" -- Relegated genes are all the genes from TGFam that had any overlap with
                       Pseudogenes with more than {relaxation} disabling mutations. These will
                       be accepted as true Pseudogenes, maintaining Shiu's annotation and
                       removing them from TGFam's final results.
            """
            f_debug.write(log) # Write to file

            ## Create outdir
            os.system(f"mkdir -p {current_dir}/tgfam_gffs")

            ## True genes (not relegated and not relaxed)
            in_list = [*{**relaxed_genes, **true_pseudogenes}.keys()]
            filter_gff(input_gff=f"{tgfamGFF}", entry_list=in_list,
                       output=f"{current_dir}/tgfam_gffs/true_genes_tgfam.gff", mode="reverse")

            ## Relaxed genes
            in_list = [*relaxed_genes.keys()]
            filter_gff(input_gff=f"{tgfamGFF}", entry_list=in_list,
                       output=f"{current_dir}/tgfam_gffs/relaxed_genes_tgfam.gff", mode="normal")

            ## Relegated genes
            in_list = [*true_pseudogenes.keys()]
            filter_gff(input_gff=f"{tgfamGFF}", entry_list=in_list,
                       output=f"{current_dir}/tgfam_gffs/relegated_genes_tgfam.gff", mode="normal")


    # Final Step
    # Create a final GFF for Shiu's Pseudogenes
    log = f"""
# STEP 4

Produce a final GFF final containing the remaining Pseudogenes that were accepted as Pseudogenes.

For this file we will remove the Pseudogenes that have been promoted to genes (true and relaxed)
based the intersection with TGFam gene families annotations

Output at: 00_final_Pseudogenes/bnut_Pseudogenes.gff
    """
    f_debug.write(log) # Write to file

    ## Create outdir
    os.system(f"mkdir -p {outdir}/00_final_Pseudogenes")

    ## Save Pseudogenes
    ## Those that have not been promoted
    filter_gff(input_gff=f"{shiuGFF}", entry_list=promoted_Pseudogenes,
               output=f"tmp.gff", mode="reverse")
    os.system(f"sort -V -k 1,4 tmp.gff > {outdir}/00_final_Pseudogenes/bnut_Pseudogenes.gff ; rm tmp.gff")


############
### Main ###
############

if __name__ == '__main__':
    arguments = docopt(__doc__, version='v1.0 by Felipe Marques de Almeida')

    #############################
    ### Checking dependencies ###
    #############################

    # Check bedtools
    try:
        subprocess.run([f"{arguments['--bedtools']}", "--help"], check = True, stdout=subprocess.DEVNULL)
    except:
        print("\n----------------------------------------")
        print("Bedtools is not available")
        print("Please, make sure it is available on PATH")
        print("or set it with --bedtools parameter (v2.29.2)")
        print("----------------------------------------\n")
        exit(1)
    # Check shiu file
    try:
        f = open(arguments['--shiu'])
        f.close()
    except FileNotFoundError:
        print("\n----------------------------------------")
        print("Input (shiu gff) file does not exist")
        print("Please, check out the help message below")
        print("----------------------------------------\n")
        os.system("./tgfam_vs_shiu-pipeline.py.py -h")
        exit(1)

    # Check tgfam gffs
    try:
        len(os.listdir(arguments['--tgfam'])) == 0
    except:
        print("\n----------------------------------------")
        print("TGFam directory does not exist or is empty")
        print("Please, check out the help message below")
        print("----------------------------------------\n")
        os.system("./tgfam_vs_shiu-pipeline.py.py -h")
        exit(1)

    ########################
    ### Summary pipeline ###
    ########################
    if arguments['summary'] and arguments['--tgfam'] and arguments['--shiu']:

        # Run
        summary_gffs(shiuGFF=arguments['--shiu'], tgfamDIR=arguments['--tgfam'],
                     bedtools=arguments['--bedtools'], relaxation=arguments['--mutations'],
                     overlap_fraction=arguments['--tgfam_fraction'])
        print("Finished!\n\nA summary of the intersections between TGFam and Shiu's Pseudogene annotation was written in: ./Intersection_summary.md\n")


    ###########################
    ### Comparison pipeline ###
    ###########################
    elif arguments['compare'] and arguments['--shiu'] and arguments['--tgfam']:

        # Run
        print(f"""

    > BEGIN

-----------------------------------------------------------------------
Starting the comparison between TGFam and Shiu's Pseudogene annotations

Please be patient, this may take a while.

Outputs will be in {arguments['--outdir']}.
-----------------------------------------------------------------------
""")
        gff_compare(shiuGFF=arguments['--shiu'], tgfamDIR=arguments['--tgfam'],
                    bedtools=arguments['--bedtools'], relaxation=arguments['--mutations'],
                    outdir=arguments['--outdir'], overlap_fraction=arguments['--tgfam_fraction'])
        print(f"""

    > END

-----------------------------------------------------------------------
The execution has finished
-----------------------------------------------------------------------
""")

    ## None
    else:
        print("----------------------------------------")
        print("Missing mandatory arguments")
        print("Please, check out the help message below")
        print("----------------------------------------\n")
        os.system("./tgfam_vs_shiu-pipeline.py.py -h")
