#!/usr/bin/python3

## Help
"""parse_resfamsDB.py

Usage:
    make_UpSetR_from_tsv.py -i <Input.tsv> -x <Resfams.xlsx> -o <output.txt>

Options:
    -h --help    Show this screen
    -i, --input  {prefix}_resistance.tsv input File
    -x, --xlsx   Resfams metadata file
    -o, --output Output File

"""
from docopt import docopt

if __name__ == '__main__':
    args = docopt(__doc__)

## Import pandas with package
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

## Load our Annotated Features and Database with pandas
resistance_features = pd.read_csv(args['<Input.tsv>'], sep="\t")
resfams_metadata = pd.read_excel(args['<Resfams.xlsx>'], sheet_name=0, header=0, index_col=False, keep_default_na=True)

## Filter our results that are from resfams
our_resfams_filtered = resistance_features[resistance_features.Prokka_inference.str.contains('resfams', case=False)]
## Split the Prokka Prokka_inference column to get resfams id
clist = our_resfams_filtered.Prokka_inference.str.split(",")
our_data = []
for i in clist:
    our_data.append(i[1].split(":")[2])

## To upper
our_data = [x.upper() for x in our_data]

## Adding Resfams ID to dataframe
our_resfams_filtered['ResfamID'] = our_data

## Slice Resfams Metadata related to our Annotation
resfams_metadata_filtered = resfams_metadata[resfams_metadata['ResfamID'].isin(our_data)]
result = pd.merge(our_resfams_filtered, resfams_metadata, on='ResfamID', how='inner')

## Select wanted columns
col = ['seqname', 'Prokka_ID', 'ResfamID', 'Description', 'CARD_ARO_Updated', 'Antibiotic Classification (Resfam Only)', 'Mechanism Classification']
final = result[col]

## Write Results to File
final.to_csv(args['--output'], sep='\t')
