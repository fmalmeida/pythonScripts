#!/usr/bin/env python3
# coding: utf-8

## Def help message
"""
A simple script for the search of intersections between transposons and pseudogene annotations.

---
Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain
Usage:
    TEs_vs_shiu-pipeline.py
    TEs_vs_shiu-pipeline.py -h|--help
    TEs_vs_shiu-pipeline.py -v|--version
    TEs_vs_shiu-pipeline.py compare [--transposons <gff> --shiu <gff> --bedtools <path> --shiu_fraction <int> --outdir <dir>]
    TEs_vs_shiu-pipeline.py plot [--parsed_gff <gff> --transposons <gff>]

Options:
    -h --help                   Show this screen.
    -v --version                Show version information
    --bedtools=<path>           Path to bedtools executable, v2.29.2 is suggested [default: bedtools]
    --shiu_fraction=<int>       Minimum overlap required as a fraction of TGFam genes in order to report a intersection [default: 0.2]
    --shiu=<gff>                GFF file containing Pseudogenes predicted with Shiu's pipeline [default: ./bnut_pseudogenes.gff]
    --transposons=<gff>         GFF containing transposon annotation [default: ./repet_bnut_final.gff]
    --parsed_gff=<gff>          GFF produced with the compare pipeline. It already is parsed and contains the TE families intersected with each predicted pseudogene [default: ./_results/bnut_pseudogenes_TEs.gff]
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
from plotnine import *

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
                               'TEs.Chr', 'TEs.Source', 'TEs.Type', 'TEs.Start', 'TEs.End', 'TEs.Score', 'TEs.Strand', 'TEs.Phase', 'TEs.Attributes', 'Overlap'])
    return df

############################################################
### Function for running bedtools intersect between GFFs ###
############################################################
def intersectBed(shiu, transposon, bedtools, overlap_fraction, out):
    os.system(f"{bedtools} intersect -wo -a {shiu} -b {transposon} -f {overlap_fraction} > {out}")

###############################
### GFF comparison pipeline ###
###############################
def gff_compare(shiuGFF, repetGFF, bedtools, outdir, overlap_fraction):

    # Remove result if existent
    if os.path.exists(f"{outdir}"):
        os.system(f"rm -rf {outdir}")

    # Create output directory
    os.system(f"mkdir -p {outdir}")

    # Debug file -- remove if existent
    if os.path.exists(f"./{outdir}/gff_comparison.log"):
        os.remove(f"./{outdir}/gff_comparison.log")
    f_debug = open(f"./{outdir}/gff_comparison.log", "a+")

    # First Step
    # Detect intersections with bedtools
    log = f"""
# STEP 1

Find intersection between REPET and Pseudogene annotations using bedtools intersect (default)
    + CMD: {bedtools} intersect -wo -a {shiuGFF} -b {repetGFF} -f {overlap_fraction} > pseudogenes_intersection.txt
    + Bedtools intersect results at: pseudogenes_intersection.txt
    """
    f_debug.write(log) # Write to file
    intersectBed(shiu=shiuGFF, transposon=repetGFF, bedtools=bedtools, overlap_fraction=f"{overlap_fraction}", out=f"{outdir}/pseudogenes_intersection.txt")

    # Second Step
    # Filter out duplicated intersections
    log = f"""
# STEP 2

Find and filter out duplicated intersections based on pseudogene IDs. From duplicated
intersections we will keep the result with the biggest overlap.

For that, we'll use a pandas funtion called: 'drop_duplicated'.

* The filtered results will be at: pseudogenes_intersection.filtered.txt
    """
    f_debug.write(log) # Write to file

    # Load results
    bedDF = intersectBed_2_DF(f"{outdir}/pseudogenes_intersection.txt")
    # Apply filter
    bedDF = bedDF.sort_values(by='Overlap', ascending=False)
    bedDF = bedDF.drop_duplicates(subset='Shiu.Attributes', keep="first")
    bedDF = bedDF.sort_values(by=['Shiu.Chr', 'Shiu.Start', 'Overlap'], ascending=True)
    # Save
    bedDF.to_csv(f"{outdir}/pseudogenes_intersection.filtered.txt", sep="\t", index=False, header=False)

    # Third step
    # Save, per ID, the transposon family for its intersection
    log = f"""
# STEP 3

Detect and save (per ID) the transposon family of each intersection with pseudogenes.

* Results will be at: transposable_pseudogenes_family.txt
    """
    f_debug.write(log) # Write to file

    # Create empty dict
    te_families = {}

    # Iterate over results
    for index, line in bedDF.iterrows():

        # Grab pseudogene id
        shiu_id   = filter(line["Shiu.Attributes"].split(';'), ['ID'])[0].split("=")[1]
        # Grab TE classification
        te_family = line["TEs.Type"]
        # Save
        if shiu_id in te_families:
            te_families[f"{shiu_id}"].append(f"{te_family}")
        else:
            te_families[f"{shiu_id}"] = [ f"{te_family}" ]

    # Write to file
    with open(f"{outdir}/transposable_pseudogenes_family.txt", "w") as outfile:
        print("Pseudogene ID\TE family", file=outfile)
        for key, value in te_families.items():
            value = ",".join(value)
            print(f"{key}\t{value}", file=outfile)

    # Fourth step
    # Re-create a GFF with TE families
    log = f"""
# STEP 4

Create a new GFF that contains the information of TE families with intersections with
the predicted pseudogenes. All will be put in the attributes column.

* Results will be at: bnut_pseudogenes_TEs.gff
    """
    f_debug.write(log) # Write to file

    # Read GFF file
    gff_df = pd.read_csv(f"{shiuGFF}", sep = "\t", comment = "#",
                        names=['Chr', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'])
    gff_df = gff_df.sort_values(by=['Chr', 'Start', 'Attributes'], ascending=True)

    # Final file
    # Select desired entries
    with open(f"{outdir}/bnut_pseudogenes_TEs.gff", 'w') as gff:
        print("## gff version 3", file=gff)

        # Iterate
        for index, line in gff_df.iterrows():

            # Grab ID
            shiu_id   = filter(line["Attributes"].split(';'), ['ID'])[0].split("=")[1]
            te_family = ",".join(te_families.get(f"{shiu_id}", ""))

            # Does it have a TE hit?
            ## Yes
            if shiu_id in [*{**te_families}.keys()]:
                 print(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7],
                 f"{line[8]};Putative_TE=Yes;TE_family={te_family}", sep="\t", file=gff)
            ## No
            else:
                print(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7],
                f"{line[8]};Putative_TE=No", sep="\t", file=gff)
    # Final step
    # Draw an image showing the TE families with
    # intersections to predicted pseudogenes
    log = f"""
# STEP 5

Draw plots of intersected TE families. Must use the 'plot' function from pipeline.

* Results will be at: *{{png,svg}}
    """
    f_debug.write(log) # Write to file

##############################################
### Plot the TE families with intersection ###
##############################################
def te_plot(finalGFF, repetGFF, outdir):

    # Remove if existent
    if os.path.exists(f"./{outdir}/TE_families_plot.png"):
        os.system(f"rm -rf {outdir}/*.png")
    if os.path.exists(f"./{outdir}/TE_families_plot.svg"):
        os.system(f"rm -rf {outdir}/*.svg")

    # Final step
    # Draw an image showing the TE families with
    # intersections to predicted pseudogenes

    # Read GFF file
    gff_df = pd.read_csv(f"{finalGFF}", sep = "\t", comment = "#",
                        names=['Chr', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'])

    # For plot
    data = []

    # Iterate
    for index, line in gff_df.iterrows():

        # It has a intersection with TE?
        putative_te = filter(line["Attributes"].split(';'), ['Putative_TE'])[0].split("=")[1]
        if putative_te == "Yes":
            # Grab TE family
            te_family = filter(line["Attributes"].split(';'), ['TE_family'])[0].split("=")[1]

            # Append fields
            data.append([line["Chr"], te_family])

    # Select
    intersect_df = pd.DataFrame(data, columns=['Chr', 'Class'])

    # Remove unwanted TE families
    plot_df = intersect_df[~intersect_df.Class.isin(['CONFUSED', 'FILTERED', 'UNK', 'SSR'])]

    # Count
    plot_df['freq'] = plot_df.groupby('Class')['Class'].transform('count')
    plot_df = plot_df.drop_duplicates(subset='Class', keep="first")

    # Sort
    plot_df = plot_df.sort_values(by='Class', ascending=True)
    plot_df['Class'] = pd.Categorical(plot_df['Class'], categories=plot_df['Class'], ordered=True)

    # Plot
    p = (ggplot(plot_df, aes(x='Class', y='freq'))
    + geom_col(aes(fill='Class'), size=20) # Percent stack bar
    #+ coord_flip()        # flipping the x- and y-axes
    + labs(title='TE families with intersections to pseudogenes', x = "", y='Count', fill='TE families')
    + theme_classic()
    + theme(axis_text_x = element_text(angle=90),
            legend_position = "none")
    ) # customizing labels
    ggsave(plot = p, filename = f"{outdir}/TE_families_plot.png")
    ggsave(plot = p, filename = f"{outdir}/TE_families_plot.svg")



    # Normalized plot
    # Percentage
    # Read REPET GFF file
    repet_df = pd.read_csv(f"{repetGFF}", sep = "\t", comment = "#",
                        names=['Chr', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'])

    # Count the total number of TEs from each family
    ## TOTAL
    tes_total = pd.DataFrame(list(repet_df['Type']), columns=['Class'])
    tes_total['freq'] = tes_total.groupby('Class')['Class'].transform('count')
    tes_total = tes_total.drop_duplicates(subset='Class', keep="first")

    ## INTERSECTS
    tes_intersected = pd.DataFrame(list(intersect_df['Class']), columns=['Class'])
    tes_intersected['freq'] = tes_intersected.groupby('Class')['Class'].transform('count')
    tes_intersected = tes_intersected.drop_duplicates(subset='Class', keep="first")


    # Merge
    plot_df = pd.merge(left = tes_total, right = tes_intersected,
                       on = "Class")
    # Normalize
    plot_df['Norm'] = plot_df['freq_y']/plot_df['freq_x']

    # Sort
    plot_df = plot_df.sort_values(by='Class', ascending=True)
    plot_df['Class'] = pd.Categorical(plot_df['Class'], categories=plot_df['Class'], ordered=True)

    # Remove unwanted TE families
    plot_df = plot_df[~plot_df.Class.isin(['CONFUSED', 'FILTERED', 'UNK', 'SSR'])]

    # Plot
    p = (ggplot(plot_df, aes(x='Class', y='Norm'))
    + geom_col(aes(fill='Class'), size=20) # Percent stack bar
    #+ coord_flip()        # flipping the x- and y-axes
    + labs(title='TE families with intersections to pseudogenes (Normalized)', x = "", y='Count (Intersections / Total)', fill='TE families')
    + theme_classic()
    + theme(axis_text_x = element_text(angle=90),
            legend_position = "none")
    ) # customizing labels
    ggsave(plot = p, filename = f"{outdir}/TE_families_plot_normalized.png")
    ggsave(plot = p, filename = f"{outdir}/TE_families_plot_normalized.svg")

############
### Main ###
############

if __name__ == '__main__':
    arguments = docopt(__doc__, version='v1.0 by Felipe Marques de Almeida')

    ###########################
    ### Comparison pipeline ###
    ###########################
    if arguments['compare'] and arguments['--shiu'] and arguments['--transposons']:

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
            os.system("./TEs_vs_shiu-pipeline.py -h")
            exit(1)

        # Check REPET file
        try:
            f = open(arguments['--transposons'])
            f.close()
        except FileNotFoundError:
            print("\n----------------------------------------")
            print("Input (repet gff) file does not exist")
            print("Please, check out the help message below")
            print("----------------------------------------\n")
            os.system("./TEs_vs_shiu-pipeline.py -h")
            exit(1)

        # Run
        print(f"""

    > BEGIN

-----------------------------------------------------------------------
Starting the comparison between transposons and pseudogenes annotations

Please be patient, this may take a while.

Outputs will be in {arguments['--outdir']}.
-----------------------------------------------------------------------
""")
        gff_compare(shiuGFF=arguments['--shiu'], repetGFF=arguments['--transposons'],
                    bedtools=arguments['--bedtools'], outdir=arguments['--outdir'],
                    overlap_fraction=arguments['--shiu_fraction'])
        print(f"""

    > END

-----------------------------------------------------------------------
The execution has finished
-----------------------------------------------------------------------
""")

    #########################
    ### plotting pipeline ###
    #########################
    elif arguments['plot'] and arguments['--parsed_gff']:

        # Check REPET file
        try:
            f = open(arguments['--transposons'])
            f.close()
        except FileNotFoundError:
            print("\n----------------------------------------")
            print("Input (repet gff) file does not exist")
            print("Please, check out the help message below")
            print("----------------------------------------\n")
            os.system("./TEs_vs_shiu-pipeline.py -h")
            exit(1)

        # Run
        print(f"""

    > BEGIN

---------------------------------------------
Generating the final plot, please be patient.

Outputs will be in {arguments['--outdir']}.
---------------------------------------------
""")
        te_plot(finalGFF=arguments['--parsed_gff'], outdir=arguments['--outdir'],
                repetGFF=arguments['--transposons'])
        print(f"""

    > END

---------------------------------------------
The plotting has finished
---------------------------------------------
""")

    ## None
    else:
        print("----------------------------------------")
        print("Missing mandatory arguments")
        print("Please, check out the help message below")
        print("----------------------------------------\n")
        os.system("./TEs_vs_shiu-pipeline.py -h")
