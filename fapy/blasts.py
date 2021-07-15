#!/usr/bin/env python3
# coding: utf-8

## Def help message
usage_blasts = """
A simple script to automatize the execution and filtering of blast alignments.

---
Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain
Usage:
    fa-py blasts
    fa-py blasts -h|--help
    fa-py blasts -v|--version
    fa-py blasts ( --query <fasta> --subject <subject> ) [ --task <task> --minid <int> --mincov <int> --culling_limit <int> --out <string> --threads <int> --2way ]

Options:
    -h --help                   Show this screen.
    -v --version                Show version information
    --task=<task>               Select which task to run (blastn, blastp, tblastn or blastx) [Default: blastn].
    --2way                      Sets the pipeline to filter alignments by coverage in a 2way manner.
                                Which means an alignment must cover at least n from the query and
                                subject lengths. Otherwise it just needs to cover n from query seq.
                                This method is good when comparing query genes to subject genes.
    --query=<fasta>             Query fasta file.
    --subject=<subject>         Subject fasta file.
    --minid=<int>               Min. Identity percentage for gene annotation [Default: 80]
    --mincov=<int>              Min. Covereage for gene annotation [Default: 80]
    --culling_limit=<int>       Blast culling_limit for best hit only [Default: 1]
    --out=<string>              File for saving blast outputs [Default: out.blast]
    --threads=<int>             Number of threads to be used [Default: 1]
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt
import pandas as pd
import os
import sys

######################
### BLAST FUNCTION ###
######################
def blast(task, query, subject, culling, minid, mincov, out, threads, twoway):

    # Outfmt
    outfmt="6 qseqid qstart qend qlen sseqid sstart send slen evalue length pident gaps gapopen bitscore"

    # format header
    os.system(f"echo \"qseqid\tqstart\tqend\tqlen\tsseqid\tsstart\tsend\tslen\tevalue\tlength\tpident\tgaps\tgapopen\tbitscore\" > {out}")

    # format db
    if task == "blastn" or task == "tblastn":
        db_type = "nucl"
    elif task == "blastx" or task == "blastp":
        db_type = "prot"
    os.system(f"makeblastdb -in {subject} -parse_seqids -out ./FA-PY-SUBJECT-DB -dbtype {db_type} 1> /dev/null")

    # run blast
    if twoway:
        os.system(f"{task} -query {query} -db ./FA-PY-SUBJECT-DB -outfmt \"{outfmt}\" -num_threads {threads} -culling_limit {culling} | \
        awk -v minid={minid} -v mincov={mincov} '{{ if ($11 >= minid && (($10 - $12) / $8 * 100) >= mincov && (($10 - $12) / $4 * 100) >= mincov) {{print $0}}  }}' >> {out} ")
    else:
        os.system(f"{task} -query {query} -db ./FA-PY-SUBJECT-DB -outfmt \"{outfmt}\" -num_threads {threads} -culling_limit {culling} | \
        awk -v minid={minid} -v mincov={mincov} '{{ if ($11 >= minid && (($10 - $12) / $4 * 100) >= mincov) {{print $0}}  }}' >> {out} ")

    # clear work dir
    os.system("rm -rf ./FA-PY-SUBJECT-DB*")

########################
### Summary function ###
########################
def summary(output):

    # Outfmt
    columns="QUERY\tQUERY_START\tQUERY_END\tQUERY_STRAND\t%QUERY_COV\tSUBJECT\tSUBJECT_START\tSUBJECT_END\tSUBJECT_STRAND\t%SUBJECT_COV\t%IDENTITY\tGAPS"

    # Summary
    blast = pd.read_csv(output, sep="\t")
    print(columns)
    for index, line in blast.iterrows():
        # Query strand
        if (line['qstart'] > line['qend']):
            strand="-"
        else:
            strand="+"
        # Subject strand
        if (line['sstart'] > line['send']):
            sstrand='-'
        else:
            sstrand='+'
        # Parse headers
        subject=line["sseqid"]
        # Query coverage
        qcov=round((100 * (line["length"] - line["gaps"]) / line["qlen"]), 2)
        # Subject coverage
        scov=round((100 * (line["length"] - line["gaps"]) / line["slen"]), 2)
        # Identity
        id=round(line["pident"], 2)
        # Gaps
        gaps=str(line["gapopen"]) + "/" + str(line["gaps"])

        # Print
        print(line["qseqid"], line["qstart"], line["qend"], strand, qcov,
              subject, line["sstart"], line["send"], sstrand, scov, id, gaps, sep = "\t")
