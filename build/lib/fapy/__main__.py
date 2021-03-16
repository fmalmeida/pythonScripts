#!/usr/bin/env python3
license="""
Copyright 2021 Felipe Almeida (almeidafmarques@gmail.com)
https://github.com/fmalmeida/pythonScripts

This file is part of my custom python scripts (fa-py) package, which is free: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. This package is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with fa-py package.
If not, see <http://www.gnu.org/licenses/>.
"""

## Def main help
usage="""
fa-py: a package to the simple distribution of my custom scripts.

Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain

usage:
    fa-py [ -h|--help ] [ -v|--version ] [ --license ]
    fa-py <command> [ -h|--help ] [ <args>... ]

options:
    -h --help                                               Show this screen
    -v --version                                            Show version information
    --license                                               Show LEGAL LICENSE information

commands:
    tsv2markdown                                            Command for rapid convertion of tsv or csv to markdown tables.
    splitgbk                                                Command to split multisequence genbank files into individual files.
    blasts                                                  Command to execute automatized blast commands.

Use: `fa-py <commmand> -h` to get more help and see examples.
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt
from subprocess import call
import sys

########################
### Import functions ###
########################
from .version import *
from .tsv2markdown import *
from .splitgbk import *
from .blasts import *

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
            file2mw(arguments['--csv'], ',', args['--header'])

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
