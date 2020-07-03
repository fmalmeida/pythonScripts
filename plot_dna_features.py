#!/usr/bin/env python
# coding: utf-8

## Def help message
"""
A simple script to plot DNA features using the DNA features python package
Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain

Usage:
    plot_dna_features.py
    plot_dna_features.py -h|--help
    plot_dna_features.py -v|--version
    plot_dna_features.py (--input <gff> | --fofn <file>) (--start <start_base> --end <end_base> --contig <contig_name>) [--feature <feature_type> --title <title> --label <label> --color <color> --output <png_out> --width <width> --height <height>]

Options:
    -h --help                   Show this screen.
    -v --version                Show version information
    --input=<gff>               Used to plot dna features from a single GFF file
    --fofn=<file>               Used to plot dna features from multiple GFF files. Contents must be in csv format with 3 columns:
                                gff,custom_label,color (HEX format). Features from each GFF will have the color set in the 3rd column,
                                labeled as -> 'custom_label: gene id'.
    --start=<start_base>        Starting position for plotting.
    --end=<end_base>            Ending position for plotting.
    --contig=<contig_name>      Name of the contig which you want to plot.
    --title=<title>             Plot title [default: Gene Plot].
    --label=<label>             Custom label for plotting. Legends will be in the following format: 'Custom label: gene id' [default: Gene].
    --feature=<feature_type>    Type of the GFF feature (3rd column) which you want to plot [default: gene].
    --color=<color>             HEX entry for desired plotting color [default: #ccccff].
    --width=<width>             Plot width ratio [default: 20].
    --height=<height>           Plot height ratio [default: 5].
    --output=<png_out>          Output PNG filename [default: ./out.png].
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt
from dna_features_viewer import *
import Bio.SeqIO
import pprint
from BCBio import GFF
from BCBio.GFF import GFFExaminer
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt



######################################################
### Function for execution with a single GFF input ###
######################################################

def single_gff(infile, start, end, contig, feature, coloring, custom_label, outfile, plot_title, plot_width, plot_height):
    # Subset GFF based on chr and feature type
    limit_info = dict(
            gff_id   = [contig],
            gff_type = [feature]
    )

    # Load GFF and its sequences
    gff = GFF.parse(infile, limit_info=limit_info)

    # Create empty features and legened list
    features = []

    ## Populate features list
    ## Filtering by location
    start_nt = int(start)
    end_nt   = int(end)
    length   = end_nt - start_nt

    for rec in gff:
        for i in range(0, len(rec.features)):
            if ( int(rec.features[i].location.start) >= int(start_nt) and int(rec.features[i].location.end) <= int(end_nt) ):

                if (str(rec.features[i].location.strand) == "+"):
                    strand=+1
                else:
                    strand=-1

                # Create input string
                ## Label in the gene plot
                # input= GraphicFeature(start=int(rec.features[i].location.start), end=int(rec.features[i].location.end),
                #                       strand=int(strand), label="{0}: {1}".format(custom_label, str(rec.features[i].qualifiers['Name'][0])), color=coloring)

                ## Label not in the gene plot
                input= GraphicFeature(start=int(rec.features[i].location.start), end=int(rec.features[i].location.end),
                                      strand=int(strand), label=str(rec.features[i].qualifiers['Name'][0]), color=coloring)

                # Append
                features.append(input)


    # Draw plot
    record = GraphicRecord(sequence_length=length, features=features, first_index=start_nt)
    ax, _ = record.plot(figure_width=int(plot_width), figure_height=int(plot_height))
    ## Label not in the gene plot (using separate legend box)
    legend = ax.legend(handles=[mpatches.Patch(facecolor=coloring, label="{0}".format(custom_label), linewidth = 0.5, edgecolor = 'black')],
                       loc = 1, title=plot_title, fontsize = 'medium', fancybox = True)
    ax.figure.savefig(outfile, bbox_inches='tight')



#######################################################
### Function for execution with multiple GFF inputs ###
#######################################################

def multiple_gff(input_fofn, start, end, contig, feature, outfile, plot_title, plot_width, plot_height):

    # Open list of filenames containing GFFs
    file = open(input_fofn, 'r')
    content = file.readlines()

    # Create empty features and legened list
    features = []
    legend_entries = []

    # Begin Parsing
    for line in content:
        data = line.strip().split(",", 3)
        infile    = data[0]
        labeling  = data[1]
        coloring  = data[2]

        # Subset GFF based on chr and feature type
        limit_info = dict(
                gff_id   = [contig],
                gff_type = [feature]
        )

        # Load GFF and its sequences
        gff = GFF.parse(infile, limit_info=limit_info)

        ## Populate features list
        ## Filtering by location
        start_nt = int(start)
        end_nt   = int(end)
        length   = end_nt - start_nt

        for rec in gff:
            for i in range(0, len(rec.features)):
                if ( int(rec.features[i].location.start) >= int(start_nt) and int(rec.features[i].location.end) <= int(end_nt) ):

                    if (str(rec.features[i].location.strand) == "+"):
                        strand=+1
                    else:
                        strand=-1

                    # Create input string
                    ## Label in the gene plot
                    # input= GraphicFeature(start=int(rec.features[i].location.start), end=int(rec.features[i].location.end),
                    #                       strand=int(strand), label="{0}: {1}".format(labeling, str(rec.features[i].qualifiers['Name'][0])), color=coloring)

                    ## Label not in the gene plot
                    input= GraphicFeature(start=int(rec.features[i].location.start), end=int(rec.features[i].location.end),
                                          strand=int(strand), label=str(rec.features[i].qualifiers['Name'][0]), color=coloring)

                    # Append DNA features plot
                    features.append(input)

        # Append to legend
        legend_entries.append(
            mpatches.Patch(facecolor=coloring, label="{0}".format(labeling), linewidth = 0.5, edgecolor = 'black')
        )


    # Draw plot
    record = GraphicRecord(sequence_length=length, features=features, first_index=start_nt)
    ax, _ = record.plot(figure_width=int(plot_width), figure_height=int(plot_height))
    ## Label not in the gene plot (using separate legend box)
    legend = ax.legend(handles=legend_entries, loc = 1, title=plot_title, fontsize = 'medium', fancybox = True)
    ax.figure.savefig(outfile, bbox_inches='tight')





## Main
if __name__ == '__main__':
    arguments = docopt(__doc__, version='v1.0 by Felipe Marques de Almeida')

    ## Single GFF
    if arguments['--input'] and arguments['--start'] and arguments['--end'] and arguments['--contig']:
        print("Executing the pipeline for a single GFF input")
        single_gff(infile=arguments['--input'], start=arguments['--start'], end=arguments['--end'],
                   contig=arguments['--contig'], feature=arguments['--feature'], coloring=arguments['--color'],
                   custom_label=arguments['--label'], outfile=arguments['--output'], plot_title=arguments['--title'],
                   plot_width=arguments['--width'], plot_height=arguments['--height'])
        print("Done, checkout the results in {}".format(arguments['--output']))

    ## Multiple GFFs
    elif arguments['--fofn'] and arguments['--start'] and arguments['--end'] and arguments['--contig']:
        print("Executing the pipeline for multiple GFF inputs")
        multiple_gff(input_fofn=arguments['--fofn'], start=arguments['--start'], end=arguments['--end'],
                     contig=arguments['--contig'], feature=arguments['--feature'], outfile=arguments['--output'],
                     plot_title=arguments['--title'], plot_width=arguments['--width'], plot_height=arguments['--height'])
        print("Done, checkout the results in {}".format(arguments['--output']))

    ## None
    else:
        print("Missing mandatory arguments")
        print("Please, check out the help message")
        print("")
        print(arguments)
