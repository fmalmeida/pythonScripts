#!/usr/bin/env python3
license="""
Copyright 2021 Felipe Almeida (almeidafmarques@gmail.com)
https://github.com/fmalmeida/pythonScripts

This file is part of my custom python scripts (falmeida-py) package, which is free: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. This package is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with falmeida-py package.
If not, see <http://www.gnu.org/licenses/>.
"""

## Def main help
usage="""
falmeida-py: a package to the simple distribution of my custom scripts.

Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain

usage:
    falmeida-py [ -h|--help ] [ -v|--version ] [ --license ]
    falmeida-py <command> [ -h|--help ] [ <args>... ]

options:
    -h --help                  Show this screen
    -v --version               Show version information
    --license                  Show LEGAL LICENSE information

commands:
    tsv2markdown               Command for rapid convertion of tsv or csv to markdown tables.
    splitgbk                   Command to split multisequence genbank files into individual files.
    align2subsetgbk            Command to subset genbank files based on alignments to a FASTA file.
    gbk2fasta                  Command to convert genbank files to fasta files.
    blasts                     Command to execute automatized blast commands.
    replace_fasta_seq          Command to replace strings in a FASTA using defitinitions from a BED file
    mpgap2csv                  Command to summarize main mpgap multiqc assembly statistics into a CSV file
    bacannot2json              Command to summarize main bacannot annotation results into JSON file

Use: `falmeida-py <commmand> -h` to get more help and see examples.
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt

########################
### Import functions ###
########################
from .version import *
from .tsv2markdown import *
from .splitgbk import *
from .gbk2fasta import *
from .blasts import *
from .align2subsetgbk import *
from .replace_fasta_seq import *
from .bacannot2json import usage_bacannot2json,bacannot2json
from .mpgap2csv import usage_mpgap2csv,mpgap2csv

## Defining main
def main():
    # Parse docopt
    __version__ = get_version()
    arguments = docopt(usage, version=__version__, help=False, options_first=True)

    ############################
    ### tsv2markdown command ###
    ############################
    if arguments['<command>'] == 'tsv2markdown':

        # Parse docopt
        args = docopt(usage_tsv2markdown, version=__version__, help=False)

        # run script
        if args['--help']:
            print(usage_tsv2markdown.strip())

        elif args['--tsv']:
            file2mw(args['--tsv'], '\t', args['--header'])

        elif args['--csv']:
            file2mw(args['--csv'], ',', args['--header'])

        else:
            print(usage_tsv2markdown.strip())

    #########################
    ### Split gbk command ###
    #########################
    elif arguments['<command>'] == 'splitgbk':
        # Parse docopt
        args = docopt(usage_splitgbk, version=__version__, help=False)

        # Run
        if args['--help']:
            print(usage_splitgbk.strip())

        elif args['--gbk']:
            splitgbk(args['--gbk'], args['--outdir'])

        else:
            print(usage_splitgbk.strip())

    ######################
    ### Blast commands ###
    ######################
    elif arguments['<command>'] == 'blasts':
        # Parse docopt
        args = docopt(usage_blasts, version=__version__, help=False)

        # Run
        if args['--help']:
            print(usage_blasts.strip())

        elif args['--query'] and args['--subject']:

            ## check if task is correct
            if args['--task'].lower() in ['blastn', 'tblastn', 'blastp', 'blastx']:
                ## run blast
                blast(task=args['--task'].lower(), query=args['--query'],
                subject=args['--subject'], culling=args['--culling_limit'],
                minid=args['--minid'], mincov=args['--mincov'], out=args['--out'],
                threads=args['--threads'], twoway=args['--2way'])
                ## summary
                summary(output=args['--out'])
            else:
                print(f"PROBLEM!\nI could not understand the task \"{args['--task'].lower()}\", please select one of: blastn, tblastn, blastp or blastx.")

        else:
            print(usage_blasts.strip())

    ######################################
    ### Subset gbk with fasta commands ###
    ######################################
    elif arguments['<command>'] == 'align2subsetgbk':
        # Parse docopt
        args  = docopt(usage_align2subsetgbk, version=__version__, help=False)

        # Run
        if args['--help']:
            print(usage_align2subsetgbk.strip())

        elif args['--gbk'] and args['--fasta']:

            # Run
            print(f"Processing file: {args['--gbk']}!")
            gbk2fasta(gbk=args['--gbk'])
            blast(task='blastn', query=args['--fasta'], subject='tmp_gbk.fa',
             culling=args['--culling_limit'], minid=args['--minid'], mincov=args['--mincov'],
             out='out.blast', threads=1, twoway=None)
            filtergbk(gbk=args['--gbk'], out=args['--out'], extension=int(args['--extension']))

            # Clean dir
            os.system(f"rm -rf tmp_gbk.fa out.blast")

        else:
            print(usage_align2subsetgbk.strip())
    
    ###################################
    ### Replace fasta seq using bed ###
    ###################################
    elif arguments['<command>'] == 'replace_fasta_seq':
        # Parse docopt
        args = docopt(usage_replace_fasta_seq, version=__version__, help=False)

        # Run
        if args['--help']:
            print(usage_replace_fasta_seq.strip())

        elif args['--fasta'] and args['--bed']:

            # Run
            print(f"Processing file: {args['--fasta']}!")
            replace_fasta_seq(input=args['--fasta'], bed=args['--bed'], output=args['--out'])
            print(f"Done!")

        else:
            print(usage_replace_fasta_seq.strip())
    
    ###########################
    ### Convert gbk 2 fasta ###
    ###########################
    elif arguments['<command>'] == 'gbk2fasta':
        # Parse docopt
        args = docopt(usage_gbk2fasta, version=__version__, help=False)

        # run script
        if args['--help']:
            print(usage_gbk2fasta.strip())

        elif args['--gbk']:
            convertgbk(genbank=args['--gbk'], genes_list=args['--fofn'])

        else:
            print(usage_gbk2fasta.strip())
    
    #############################
    ### bacannot2json command ###
    #############################
    if arguments['<command>'] == 'bacannot2json':

        # Parse docopt
        args = docopt(usage_bacannot2json, version=__version__, help=False)

        # run script
        if args['--help']:
            print(usage_bacannot2json.strip())
        
        elif args['--input']:
            bacannot2json(args['--input'], args['--output'])

        else:
            print(usage_bacannot2json.strip())

    #########################
    ### mpgap2csv command ###
    #########################
    elif arguments['<command>'] == 'mpgap2csv':

        # Parse docopt
        args = docopt(usage_mpgap2csv, version=__version__, help=False)

        # run script
        if args['--help']:
            print(usage_mpgap2csv.strip())
        
        elif args['--input']:
            mpgap2csv(args['--input'], args['--output'])

        else:
            print(usage_mpgap2csv.strip())

    #####################
    ### Check license ###
    #####################
    elif arguments['--license']:
        print(license.strip())

    #######################################
    ### Without commands nor parameters ###
    #######################################
    else:
        print(usage.strip())

## Calling main
if __name__ == '__main__':
    main()
