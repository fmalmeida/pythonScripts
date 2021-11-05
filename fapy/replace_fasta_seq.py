#!/usr/bin/env python3
#!/usr/bin/env python3
# coding: utf-8

########################
### Def help message ###
########################
usage_replace_fasta_seq = """
A script meant to replace a string in a FASTA file using defitions in a BED file.

---
Copyright (C) 2021 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain

Usage:
    fa-py replace_fasta_seq [ -h|--help ]
    fa-py replace_fasta_seq [ --fasta <fasta> --bed <bed> --out <out_fasta>  --sep <sep> ]

Options:
    -h --help                      Show this screen.
    -f --fasta=<fasta>             FASTA file to replace sequences in.
    -o --out=<out_fasta>           Output FASTA file [Default: out.fasta].
    -b --bed=<bed>                 BED file with replacement definitions.
    -s --sep=<sep>                 Separator for BED file [Default: '\t'].

Comments:
    The BED file MUST be a 4 column file with the following format:
        contig start end sub_seq
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt
from Bio import SeqIO

##########################
### load fasta as dict ###
##########################
def fasta_as_dict(fasta):
    fasta_dict = {}
    for seq_record in SeqIO.parse(fasta, "fasta"):
        fasta_dict[seq_record.id] = seq_record
    return fasta_dict

#############################################
### replace sequences using values in bed ###
#############################################
def replace_seq_in_dict(bed, dict, sep):
    with open(bed) as f:
        for line in f:
            contig, start, end, sub_seq = line.strip().split(sep)
            dict[contig].seq = dict[contig].seq[:int(start) - 1] + sub_seq + dict[contig].seq[int(end):]

####################
### output fasta ###
####################
def output_fasta(dict, out):
    with open(out, 'w') as f:
        for record in dict.values():
            print(">" + record.id, file=f)
            print(record.seq, file=f)

#####################
### main function ###
#####################
def replace_fasta_seq(input, output, bed, sep):
    fasta_dict = fasta_as_dict(input)
    replace_seq_in_dict(bed, fasta_dict, sep)
    output_fasta(fasta_dict, output)