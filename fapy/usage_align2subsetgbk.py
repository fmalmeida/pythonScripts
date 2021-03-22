#!/usr/bin/env python3
# coding: utf-8

########################
### Def help message ###
########################
usage_align2subsetgbk = """
A script meant to subset a genbank annotation file based on alignments
against a query FASTA file

---
Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain

Usage:
    fa-py align2subsetgbk [ -h|--help ]
    fa-py align2subsetgbk [ --gbk <in_gbk> --fasta <fasta> --out <out_gbk> --minid <int> --mincov <int> --culling_limit <int> --extension <int> ]

Options:
    -h --help                      Show this screen.
    -g --gbk=<in_gbk>              Gbk file for subset
    -f --fasta=<fasta>             FASTA file for querying the gbk
    -o --out=<out_gbk>             Gbk filtered output file [Default: out.gbk].
    --extension=<int>              Base pair length to extend the flank regions in the alignment [Default: 0].
    --minid=<int>                  Min. Identity percentage for gene annotation [Default: 80].
    --mincov=<int>                 Min. Covereage for gene annotation [Default: 80].
    --culling_limit=<int>          Blast culling_limit for best hit only [Default: 1].
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt
from Bio import SeqIO
from io import StringIO
import pandas as pd
import os
import sys
from .blasts import *

##########################################
### Function to convert gbk into fasta ###
##########################################
def gbk2fasta(gbk):
    f=open("tmp_gbk.fa", "a")
    for seq_record in SeqIO.parse(gbk, 'genbank'):
        print(f">{seq_record.id}\n{seq_record.seq}\n", file=f)

############################################
### Function to filter gbk based on hits ###
############################################
def filtergbk(gbk, out, extension):

    f=open(f"{out}", "w")
    # Read blast results
    blast_res = pd.read_csv('out.blast', sep = '\t')
    print(blast_res)

    # Get locus_tags
    contigs = sorted(set(blast_res["sseqid"].tolist()))

    # Subset
    filtered = []
    for contig in contigs:
        for seq_record in SeqIO.parse(gbk, 'genbank'):
            if str(contig) == str(seq_record.id):
                small_df = blast_res[blast_res['sseqid'].isin([seq_record.id])]
                for index, row in small_df.iterrows():
                    for features in seq_record.features:
                        # plus strand
                        if int(row["sstart"]) <= int(row["send"]):
                            if int(features.location.start) >= int(row["sstart"] - extension) and int(features.location.start) <= int(row["send"] + extension) and features.type != "source":
                                filtered.append(features)
                        # minus strand
                        elif int(row["sstart"]) >= int(row["send"]):
                            if int(features.location.start) >= int(row["send"] - extension) and int(features.location.start) <= int(row["sstart"] + extension) and features.type != "source":
                                filtered.append(features)
            else:
                pass

        # Print results
        seq_record.features = filtered
        if len(seq_record.features) > 0:
            SeqIO.write(seq_record, f, 'gb')
