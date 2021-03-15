#!/usr/bin/python3

## Help
"""make_UpSetR_from_tsv.py

Usage:
    parse_cardDB.py -i <Input.tsv> --aro_index <aro_index.tsv> \
--aro_categories_index <aro_categories_index.tsv> -o <output.txt>

Options:
    -h --help    Show this screen
    -i, --input  {prefix}_resistance.tsv input File
    --aro_index  card aro_index file
    --aro_categories_index  card aro_categories_index file
    -o, --output Output File

"""
from docopt import docopt

if __name__ == '__main__':
    args = docopt(__doc__)

## Import pandas with package
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

## Load CARD db metadata
aro_index = pd.read_csv(args['<aro_index.tsv>'], sep="\t").drop(columns=['ARO Name']).sort_values(by=['ARO Accession'])
aro_categories_index = pd.read_csv(args['<aro_categories_index.tsv>'], sep="\t")

### Merge them all
card_metadata = pd.merge(aro_index, aro_categories_index, on='Protein Accession', how='inner')

## Load our Annotated Features and Database with pandas
resistance_features = pd.read_csv(args['<Input.tsv>'], sep="\t")
resistance_features.rename(columns={'ARO_Accession': 'ARO Accession'}, inplace=True)
resistance_features['ARO Accession'] = [x.upper() for x in resistance_features['ARO Accession']]

## Merge Description with Our data
complete = pd.merge(card_metadata, resistance_features, on='ARO Accession', how='inner')

## Select wanted columns
col = ['seqname', 'Prokka_product', 'ARO Accession', 'Protein Accession', 'Drug Class', 'Resistance Mechanism', 'AMR Gene Family']
final = complete[col]

## Write Results to File
final.to_csv(args['--output'], sep='\t')
