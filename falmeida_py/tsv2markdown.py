#!/usr/bin/env python
# coding: utf-8

## Def help message
usage_tsv2markdown = """
A simple script to convert tsv (or csv) files to markdown tables using tabulate!
Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain

usage:
    falmeida-py tsv2markdown
    falmeida-py tsv2markdown [ -h|--help ]
    falmeida-py tsv2markdown [ --tsv <file> --csv <file> --header <list> ]

options:
    -h --help                               Show this screen.
    --tsv=<file>                            Input tsv file to print as markdown table
    --csv=<file>                            Input csv file to print as markdown table
    --header=<list>                         If file does not have a header, set a
                                            custom header. E.g. --header "Planet,R (km),mass (x 10^29 kg)".
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt
from tabulate import tabulate

###########################
### Read header as list ###
###########################
def header2list(header):
    return header.split(',')

###################################
### Convert with header in file ###
###################################
def file2mw(filep, fsep, header):

    # read data.frame
    with open(filep) as file_in:
        lines = []
        for line in file_in:
            words = line.split(fsep)
            words_striped = [word.strip() for word in words]
            lines.append(words_striped)

    # generate tabulate
    if header:
        print(tabulate(lines, headers=header2list(header), tablefmt="github"))
    else:
        print(tabulate(lines, headers="firstrow", tablefmt="github"))
