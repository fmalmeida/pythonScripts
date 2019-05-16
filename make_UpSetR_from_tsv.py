#!/home/falmeida/falmeidaEXT/miniconda/bin/python

## Help
"""make_UpSetR_from_tsv.py

Usage:
    make_UpSetR_from_tsv.py -i <OrthoGroups.tsv> -o <output.txt>

Options:
    -h --help    Show this screen
    -i, --input  OrthoGroups.tsv input File
    -o, --output Output File

"""
from docopt import docopt

if __name__ == '__main__':
    args = docopt(__doc__)

## Import pandas with package
import pandas as pd

## Load our OrthoGroups.tsv with pandas
orthogroups = pd.read_csv(args['<OrthoGroups.tsv>'], sep="\t")

## Fill empty with 0
orthogroups.fillna(0, inplace=True)

## Get DataFrame Column Names
cols = list(orthogroups)

## Fill those not empty with 1
orthogroups.loc[orthogroups[cols[1]] != 0, cols[1]] = 1
orthogroups.loc[orthogroups[cols[2]] != 0, cols[2]] = 1

## Write Results to File
orthogroups.to_csv(args['--output'], sep='\t')
