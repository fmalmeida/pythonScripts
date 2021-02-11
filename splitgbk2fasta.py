#!/usr/bin/env python3
# coding: utf-8

########################
### Def help message ###
########################
"""
A script meant to subset a genbank annotation file based on alignments
against a query FASTA file

---
Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain

Usage:
    splitgbk2fasta.py
    splitgbk2fasta.py -h|--help
    splitgbk2fasta.py -v|--version
    splitgbk2fasta.py ( --gbk <gbk> --fasta <fasta> --out <gbk> ) [ --minid <int> --mincov <int> --culling_limit <int> ]

Options:
    -h --help                      Show this screen.
    -v --version                   Show version information
    -g --gbk=<gbk>                 Gbk file for subset
    -f --fasta=<fasta>             FASTA file for querying the gbk
    -o --out=<gbk>                 Gbk filtered output file
    --minid=<int>                  Min. Identity percentage for gene annotation [default: 80]
    --mincov=<int>                 Min. Covereage for gene annotation [default: 80]
    --culling_limit=<int>          Blast culling_limit for best hit only [default: 1]
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

##########################################
### Function to convert gbk into fasta ###
##########################################
def gbk2fasta(gbk):
    f=open("tmp_gbk.fa", "a")
    for seq_record in SeqIO.parse(gbk, 'genbank'):
        for seq_feature in seq_record.features:
            if seq_feature.type=="CDS":
                if 'translation' in seq_feature.qualifiers:
                    print(f">{seq_feature.qualifiers['locus_tag'][0]}", file=f)
                    print(f"{seq_feature.qualifiers['translation'][0]}", file=f)
                else:
                    start    = seq_feature.location.nofuzzy_start
                    end      = seq_feature.location.nofuzzy_end
                    tl_table = seq_feature.qualifiers['transl_table'][0]
                    if str(seq_feature.strand) == '-1':
                        my_seq = seq_record.seq[start:end].reverse_complement().translate(table=tl_table)
                    else:
                        my_seq = seq_record.seq[start:end].translate(table=tl_table)
                    seq_feature.qualifiers['translation'] = str(my_seq)
                    print(f">{seq_feature.qualifiers['locus_tag'][0]}", file=f)
                    print(f"{seq_feature.qualifiers['translation'][0]}", file=f)

##############################################################
### Function to run blast and detect locus_tags that match ###
##############################################################
def blastgbk(fasta, culling, minid, mincov):
    # Outfmt
    outfmt="6 qseqid qstart qend qlen sseqid sstart send slen evalue length pident gaps gapopen stitle"

    # Run blast
    os.system(f"makeblastdb -dbtype nucl -in {fasta} -out query_db")
    os.system(f"tblastn -db query_db -query tmp_gbk.fa -outfmt \"{outfmt}\" -culling_limit {culling} | \
        awk -v minid={minid} -v mincov={mincov} '{{ if ($11 >= minid && (($10 - $12) / $4 * 100) >= mincov) {{print $0}}  }}' > out.blast")

############################################
### Function to filter gbk based on hits ###
############################################
def filtergbk(gbk, out):

    f=open(f"{out}", "w")
    # Read blast results
    blast_res = pd.read_csv('out.blast', sep = '\t', header=None)

    # Get locus_tags
    sel_locus = blast_res.iloc[:, 0].tolist()

    # Subset
    for seq_record in SeqIO.parse(gbk, 'genbank'):

        cds = [feat for feat in seq_record.features if feat.type == 'CDS']
        filtered = [ft for ft in cds if ft.qualifiers['locus_tag'][0] in sel_locus]
        seq_record.features = [i for i in filtered]
        if len(seq_record.features) > 0:
            SeqIO.write(seq_record, f, 'gb')

############
### Main ###
############
if __name__ == '__main__':
    arguments = docopt(__doc__, version='v1.0 by Felipe Marques de Almeida')

    ## Run pipeline
    if arguments['--gbk'] and arguments['--fasta']:

        # Run
        gbk2fasta(gbk=arguments['--gbk'])
        blastgbk(fasta=arguments['--fasta'], culling=arguments['--culling_limit'],
                 minid=arguments['--minid'], mincov=arguments['--mincov'])
        filtergbk(gbk=arguments['--gbk'], out=arguments['--out'])

        # Clean dir
        os.system(f"rm tmp_gbk.fa out.blast query_db.n*")

    ## None
    else:
        print("Missing mandatory arguments")
        print("Please, check out the help message")
        print("")
        print(arguments)
