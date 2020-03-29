#!/usr/bin/env python3

'''MongoDB_parse_JSON.py last modified 2020-03-29

MongoDB_parse_JSON.py -i [input.json] -n [Collection Name]

A pipeline to parse JSON files into MongoDB Documents.
Documents are always inserted into a Database named: Annotation
'''

#
import sys
import argparse
import time
from collections import namedtuple
from pymongo import MongoClient
import json
#

def run_mongoDB(input, collection_name):
    # Create connection
    client = MongoClient()

    # Create Database
    db = client['Annotation']

    # Create Collection
    collection = db[collection_name]

    # Input Document
    with open(input) as f:
        file_data = json.load(f)

    collection.insert_many(file_data)


def main(argv):
    if not len(argv):
        argv.append("-h")

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument('-i', '--input',  help="JSON file to be inserted in MongoDB Collection")
    parser.add_argument('-n', '--name',   help="MongoDB Collection name for input file")
    args = parser.parse_args(argv)

    # Execute the program
    run_mongoDB(args.input, args.name)

    # Print Goodbye
    print("\nFile %s added into MongoDB Database (Annotation) under collection %s\n"
    % (args.input, args.name))


if __name__ == "__main__":
    main(sys.argv[1:])
