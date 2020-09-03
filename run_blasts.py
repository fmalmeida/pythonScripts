#!/usr/bin/env python3
# coding: utf-8

## Def help message
"""
Script created exclusevely for executing blast alignments inside Bacannot pipeline.

---
Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain
Usage:
    run_blasts.py
    run_blasts.py -h|--help
    run_blasts.py -v|--version
    run_blasts.py blastn  [--query <fasta> --db <db> --minid <int> --mincov <int> --culling_limit <int> --out <string> --threads <int>]
    run_blasts.py blastp  [--query <fasta> --db <db> --minid <int> --mincov <int> --culling_limit <int> --out <string> --threads <int>]
    run_blasts.py blastx  [--query <fasta> --db <db> --minid <int> --mincov <int> --culling_limit <int> --out <string> --threads <int>]
    run_blasts.py tblastn [--query <fasta> --db <db> --minid <int> --mincov <int> --culling_limit <int> --out <string> --threads <int>]

Options:
    -h --help                   Show this screen.
    -v --version                Show version information
    --query=<fasta>             Query genome or genes to be searched.
    --db=<db>                   Blast or Diamond database to be used.
    --minid=<int>               Min. Identity percentage for gene annotation [default: 80]
    --mincov=<int>              Min. Covereage for gene annotation [default: 80]
    --culling_limit=<int>       Blast | Diamond culling_limit for best hit only [default: 1]
    --out=<string>              File for saving blast outputs [default: out.blast]
    --threads=<int>             Number of threads to be used [default: 1]
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt
import pandas as pd
import os
import sys

###########################################
### Function for filtering python lists ###
###########################################
def filter(string, substr):
    return [str for str in string if
             any(sub in str for sub in substr)]

#######################
### BLASTN function ###
#######################
def blastn(query, db, culling, minid, mincov, out, threads):

    # Outfmt
    outfmt="6 qseqid qstart qend qlen sseqid sstart send slen evalue length pident gaps gapopen stitle"

    # Run blastn
    os.system(f"echo \"qseqid\tqstart\tqend\tqlen\tsseqid\tsstart\tsend\tslen\tevalue\tlength\tpident\tgaps\tgapopen\tstitle\" > {out}")
    os.system(f"blastn -query {query} -db {db} -outfmt \"{outfmt}\" -num_threads {threads} -culling_limit {culling} -perc_identity {minid} | \
    awk -v minid={minid} -v mincov={mincov} '{{ if ($11 >= minid && (($10 - $12) / $8 * 100) >= mincov) && (($10 - $12) / $4 * 100) >= mincov) {{print $0}}  }}' >> {out} ")

########################
### TBLASTN function ###
########################
def tblastn(query, db, culling, minid, mincov, out, threads):

    # Outfmt
    outfmt="6 qseqid qstart qend qlen sseqid sstart send slen evalue length pident gaps gapopen stitle"

    # Run blastn
    os.system(f"echo \"qseqid\tqstart\tqend\tqlen\tsseqid\tsstart\tsend\tslen\tevalue\tlength\tpident\tgaps\tgapopen\tstitle\" > {out}")
    os.system(f"tblastn -query {query} -db {db} -outfmt \"{outfmt}\" -num_threads {threads} -culling_limit {culling} -perc_identity {minid} | \
    awk -v minid={minid} -v mincov={mincov} '{{ if ($11 >= minid && (($10 - $12) / $8 * 100) >= mincov) && (($10 - $12) / $4 * 100) >= mincov) {{print $0}}  }}' >> {out} ")

#######################
### BLASTX function ###
#######################
def blastx(query, db, culling, minid, mincov, out, threads):

    # Outfmt
    outfmt="6 qseqid qstart qend qlen sseqid sstart send slen evalue length pident gaps gapopen stitle"

    # Run blastn
    os.system(f"echo \"qseqid\tqstart\tqend\tqlen\tsseqid\tsstart\tsend\tslen\tevalue\tlength\tpident\tgaps\tgapopen\tstitle\" > {out}")
    os.system(f"diamond blastx --query {query} --db {db} --outfmt {outfmt} --max-target-seqs {culling} \
    --threads {threads} --id {minid} --subject-cover {mincov} --query-cover {mincov} >> {out} ")

#######################
### BLASTP function ###
#######################
def blastp(query, db, culling, minid, mincov, out, threads):

    # Outfmt
    outfmt="6 qseqid qstart qend qlen sseqid sstart send slen evalue length pident gaps gapopen stitle"

    # Run blastn
    os.system(f"echo \"qseqid\tqstart\tqend\tqlen\tsseqid\tsstart\tsend\tslen\tevalue\tlength\tpident\tgaps\tgapopen\tstitle\" > {out}")
    os.system(f"diamond blastp --query {query} --db {db} --outfmt {outfmt} --max-target-seqs {culling} \
    --threads {threads} --id {minid} --subject-cover {mincov} --query-cover {mincov} >> {out} ")

########################
### Summary function ###
########################
def summary(output):

    # Outfmt
    columns="SEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tDESCRIPTION"

    # Summary
    blast = pd.read_csv(output, sep="\t")
    print(columns)
    for index, line in blast.iterrows():
        # Query strand
        if (line['qstart'] > line['qend']):
            strand="-"
        else:
            strand="+"
        # Subsject strand
        if (line['sstart'] > line['send']):
            cov_map=str(line["send"]) + '-' + str(line["sstart"]) + "/" + str(line["slen"])
        else:
            cov_map=str(line["sstart"]) + '-' + str(line["send"]) + "/" + str(line["slen"])
        # Parse headers
        db=line["stitle"].split('~~~')[0]
        gene=line["stitle"].split('~~~')[1]
        acc=line["stitle"].split('~~~')[2]
        prodc=line["stitle"].split('~~~')[3].split(" ", 1)[0]
        desc=line["stitle"].split('~~~')[3].split(" ", 1)[1]
        # Subject coverage
        cov=round((100 * (line["length"] - line["gaps"]) / line["slen"]), 2)
        # Identity
        id=round(line["pident"], 2)
        # Gaps
        gaps=str(line["gapopen"]) + "/" + str(line["gaps"])

        # Print
        print(line["qseqid"], line["qstart"], line["qend"], strand, gene,
              cov_map, gaps, cov, id, db, acc, prodc, desc, sep = "\t")




############
### MAIN ###
############

if __name__ == '__main__':
    arguments = docopt(__doc__, version='v1.0 by Felipe Marques de Almeida')

    ##############
    ### BLASTN ###
    ##############
    if arguments['blastn'] and arguments['--query'] and arguments['--db']:
        blastn(query=arguments['--query'], db=arguments['--db'],
               culling=arguments['--culling_limit'], minid=arguments['--minid'],
               mincov=arguments['--mincov'], out=arguments['--out'],
               threads=arguments['--threads'])
        summary(output=arguments['--out'])

    ##############
    ### BLASTX ###
    ##############
    elif arguments['blastx'] and arguments['--query'] and arguments['--db']:
        blastx(query=arguments['--query'], db=arguments['--db'],
               culling=arguments['--culling_limit'], minid=arguments['--minid'],
               mincov=arguments['--mincov'], out=arguments['--out'],
               threads=arguments['--threads'])
        summary(output=arguments['--out'])

    ##############
    ### BLASTP ###
    ##############
    elif arguments['blastp'] and arguments['--query'] and arguments['--db']:
        blastp(query=arguments['--query'], db=arguments['--db'],
               culling=arguments['--culling_limit'], minid=arguments['--minid'],
               mincov=arguments['--mincov'], out=arguments['--out'],
               threads=arguments['--threads'])
        summary(output=arguments['--out'])

    ###############
    ### TBLASTN ###
    ###############
    elif arguments['tblastn'] and arguments['--query'] and arguments['--db']:
        tblastn(query=arguments['--query'], db=arguments['--db'],
               culling=arguments['--culling_limit'], minid=arguments['--minid'],
               mincov=arguments['--mincov'], out=arguments['--out'],
               threads=arguments['--threads'])
        summary(output=arguments['--out'])

    ############
    ### None ###
    ############
    else:
        print("Missing mandatory arguments")
        print("Please, check out the help message")
        print("")
        print(arguments)
