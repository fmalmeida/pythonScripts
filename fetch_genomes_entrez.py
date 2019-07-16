#!/usr/bin/python3

## Help
"""fetch_genomes_entrez.py

Usage:
    fetch_genomes_entrez.py [--show] [--help] [(--query <query-to-search> --db \
<database-to-search> --email <your-email>) --collection <ncbi-collection> --output <output-dir>]

Options:
    -h --help         Show help screen
    --query           Query for search
    --db              NCBI database name
    --output          Path to Output directory [default: ./].
    --collection      From which collection the genome must be downloaded RefSeq or Genbank [default: RefSeq].
    --email           Provide your email address: my_email@gmail.com
    --show            Shows all available Entrez databases to search

This pipeline was developed to facilitate the download of a massive
amount of genomes from NCBI based on a normal NCBI database search.

For a proper use, the user must know what he wants to ask, from which
NCBI database. For instance, one must ask 'Novosphingobium' on the
NCBI assembly database. This will retrieve all genomes related to this search.

Also, the user needs to choose from which assembly he wants the data: RefSeq,
or Genbank.

"""
from docopt import docopt

if __name__ == '__main__':
    args = docopt(__doc__)

# Import Package
from Bio import Entrez
import os

# Fix Broken defaults
if args['<output-dir>'] == None:
    args['<output-dir>'] = './'

if args['<ncbi-collection>'] == None:
    args['<ncbi-collection>'] = 'RefSeq'

# Load Parameters
email      = args['<your-email>']
query      = args['<query-to-search>']
database   = args['<database-to-search>']
outDir = args['<output-dir>']
collection = args['<ncbi-collection>']
if collection == 'RefSeq':
    from_ftp = 'FtpPath_RefSeq'
elif collection == 'Genbank':
    from_ftp = 'FtpPath_GenBank'
else:
    from_ftp = 'FtpPath_RefSeq'

# Initiate Pipeline

## Print Entrez databases
if args['--show'] == True:
    Entrez.email = "Your.Name.Here@example.org"
    handle = Entrez.einfo()
    records = Entrez.read(handle)
    handle.close

    for record in records:
        databases=records[record]
        print(record)
        print()
        print(databases)

## Execute Entrez Search
### Example Result Dict
### DictElement({'Count': '1', 'RetMax': '1', 'RetStart': '0', 'IdList': ['709711'], \
### 'TranslationSet': [DictElement({'From': 'Novosphingobium rosa', 'To': '"Novosphingobium rosa"[Organism]'},\
###  attributes={})], 'TranslationStack': [DictElement({'Term': '"Novosphingobium rosa"[Organism]', \
### 'Field': 'Organism', 'Count': '1', 'Explode': 'Y'}, attributes={}), 'GROUP'], \
### 'QueryTranslation': '"Novosphingobium rosa"[Organism]'}, attributes={})

elif args['--query'] == True:
    # First Part - Get IDs of results
    print('Now performing the search of: ', query, ' in ', database, ' database!')
    Entrez.email = email
    handle = Entrez.esearch(db=database, term=query)
    record = Entrez.read(handle)
    handle.close
    ids=record['IdList']
    # Second Part - Retrieve Document Summary of result IDs
    handle = Entrez.efetch(db=database, id=ids, rettype="docsum")
    result = Entrez.read(handle)
    # Third Part - Get ftp path of results
    ftp_path = result['DocumentSummarySet']['DocumentSummary'][0][from_ftp]
    # Fourth Part - Download Result Genomes
    name = "echo " + ftp_path + " | awk -F\"/\" '{print $0\"/\"$NF\"_genomic.fna.gz\"}'"
    ftp_file = os.popen(name).read()
    cmd = "wget -P " + str(outDir) + " " + str(ftp_file)
    os.system(cmd)
    handle.close
else:
    print('Unknown command. Please check the help!.')
    raise DocoptExit()
