#!/usr/bin/env python3
# coding: utf-8

## Def help message
usage_gbk2fasta = """
A simple script to automatize the execution and filtering of blast alignments.

---
Copyright (C) 2021 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain
Usage:
    falmeida-py gbk2fasta
    falmeida-py gbk2fasta -h|--help
    falmeida-py gbk2fasta -v|--version
    falmeida-py gbk2fasta ( --gbk <genbank> ) [ --out <fasta> --fofn <file> --type <type> ]

Options:
    -h --help                   Show this screen.
    -v --version                Show version information
    -g --gbk <genbank>          Genbank file to be converted to fasta
    -o --out <fasta>            Output fasta file [Default: stdout]
    -f --fofn <file>            File with the list of genes to be extracted. One gene per line, using the 'locus_tag' field.
    -t --type <type>            Type of sequence to output genes: nucl or prot [Default: prot]
"""

##################################
### Loading Necessary Packages ###
##################################
from Bio import SeqIO

#########################################
### load list of genes as python list ###
#########################################
def load_gene_list(genes_list):
    """
    Loads a list of genes from a file.
    """
    genes = []
    with open(genes_list, 'r') as f:
        for line in f:
            genes.append(line.strip())
    return genes

#################################################
### prints a fasta from gbk biopython feature ###
#################################################
def print_fasta(feature, genes_list):
    if genes_list == None:
        print(f">{feature.qualifiers['locus_tag'][0]} {feature.qualifiers['product'][0]}")
        print(feature.qualifiers['translation'][0])
    else:
        if feature.qualifiers['locus_tag'][0] in genes_list:
            print(f">{feature.qualifiers['locus_tag'][0]} {feature.qualifiers['product'][0]}")
            print(feature.qualifiers['translation'][0])

###########################################
### loads and converts genbank to fasta ###
###########################################
def convertgbk(genbank, genes_list):
    """
    Loads a genbank file and converts it to fasta.
    """
    if genes_list:
        genes = load_gene_list(genes_list)
    else:
        genes = None
    for record in SeqIO.parse(genbank, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                print_fasta(feature, genes)