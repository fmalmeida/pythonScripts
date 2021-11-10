#!/usr/bin/env python
# coding: utf-8

## Def help message
usage_splitgbk = """
A very simple script to split multisequence genbank files into separate files
Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain

usage:
    falmeida-py splitgbk
    falmeida-py splitgbk [ -h|--help ]
    falmeida-py splitgbk [ --gbk <file> ] [ -o|--outdir <outdir> ]

options:
    -h --help                               Show this screen.
    -g --gbk=<file>                         Input genbank file to split into multiple individual files.
    -o --outdir=<outdir>                    Directory in which to write the splitted files [Default: ./].
"""

##################################
### Loading Necessary Packages ###
##################################
from Bio import SeqIO
import os

####################
### GBK splitter ###
####################
def splitgbk(gbk, outdir):
    # just parse the dir
    if type(outdir) == list:
        outdir = outdir[0]
    # exec biopython
    for rec in SeqIO.parse(gbk, "genbank"):
        SeqIO.write([rec], open(os.path.basename(os.path.normpath(outdir)) + "/" + rec.id + ".gbk", "w"), "genbank")
    # finish
    print(f"Done!\nIndividual files have been written at: {outdir}")
