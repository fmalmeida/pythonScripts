#!/usr/bin/env python3
# coding: utf-8

## Def help message
"""
Ask mongo: This command allows you to ask/interrogate the mongo dbs in your machine

usage:
    ask-mongo.py [ -h|--help ]
    ask-mongo.py [ --dbpath <db_path> ] [ --list_dbs ]
    ask-mongo.py [ --dbpath <db_path> ] [ --db <db_name> --list_collections ]
    ask-mongo.py [ --dbpath <db_path> ] [ --db <db_name> --collection <collection_name> ] [ --overview ] [ --subfield <subfield> --key <key> --val <val> ]

options:

    -h, --help                                               Show this screen
    --dbpath=<db_path>                                       Where you have saved your dbs? Data directory of your mongo dbs. [Default: /data/db]
    --list_dbs                                               List available mongo dbs in your system
    --db=<db_name>                                           Mongo DB to be queryed
    --list_collections                                       Check the available collections in a given database
    --collection=<collection_name>                           Collection name (from mongo db) to be queryed
    --overview                                               Documents in a given collection (from a mongo db)
    --subfield=<subfield>                                    Is your key/val inside a subfield (a nested document)? Give the name.
    --key=<key>                                              Query a collection for documents based on which key?
    --val=<val>                                              What to search in this key?

example:
"""

##################################
### Loading Necessary Packages ###
##################################
import sys
import os
import re
from docopt import docopt
import urllib.request, urllib.parse, urllib.error
import json
import pymongo
from pymongo import MongoClient
import pathlib
import pprint
import random
from bson.objectid import ObjectId

########################
### Useful functions ###
########################
def start_mongod(db_path):

    try:
        # Function to start mongo shell if not started yet
        os.system(f"mongod --dbpath {db_path} --syslog --fork &> /dev/null")
    except:
        pass

#########################
### Check your mongos ###
#########################
def check_mongos():

    # Create connection
    client = MongoClient()

    # Get dbs
    dbs = client.list_database_names()

    # Print databases
    print(f"\nThe available mondo dbs found in your system are: {dbs}\n")

    # Close client
    client.close()

#######################################
### Check collections in a database ###
#######################################
def check_db(db_name):

    # Create connection
    client = MongoClient()

    # Open Database
    db = client[db_name]

    # Check available collections
    cols = db.list_collection_names()
    print(f"\nAll the available collections found in the {db_name} database are given in the list: {cols}\n")

#############################
### Overview a collection ###
#############################
def collection_overview(db_name, collection_name):

    # Create connection
    client = MongoClient()

    # Open Database
    db = client[db_name]

    # Open collection
    collection = db[collection_name]

    # Count
    n_docs = collection.count_documents({})

    # Get keys
    keys = []
    for doc in collection.find({}):
        for key in iter(doc.keys()):
            keys.append(key)
    keys = list(dict.fromkeys(keys)) # remove duplicates

    # Get one as example
    ids = []
    for doc in collection.find({}):
        ids.append(doc['_id'])
    example = collection.find_one({'_id': ObjectId(random.choice(ids))})

    # Give overview
    print(f"When analysing the collection {collection_name} in the database {db_name} we have found a total of {n_docs} documents (entries).")
    print(f"These are the searcheable keys found in these documents: {keys}.")
    print(f"They are searcheable and parseable based on its values.")
    print(f"Here it is an example of one document of the collection:\n")
    pprint.pprint(example)

####################
### Get Document ###
####################
def get_doc(db_name, collection_name, key, val, subfield):

    # Create connection
    client = MongoClient()

    # Open Database
    db = client[db_name]

    # Open collection
    collection = db[collection_name]

    # Find my document
    if key == '_id':
        results = collection.find({ key: ObjectId(val) })
    elif subfield != None:
        results = collection.find({ f"{subfield}.{key}" : val })
    else:
        results = collection.find({ key: val })

    # Print
    for doc in results:
        pprint.pprint(doc)

################
### Def main ###
################
if __name__ == '__main__':
    args = docopt(__doc__, version='v1.0 by Felipe Marques de Almeida')

    # Start db (if not started)
    start_mongod(db_path=args['--dbpath'])

    ## Help
    if args['--help']:
        print(args.strip())

    ## List dbs
    elif args['--list_dbs']:
        check_mongos()

    ## List collections
    elif args['--list_collections'] and args['--db']:
        check_db(db_name=args['--db'])

    ## Overview collection
    elif args['--collection'] and args['--db']:
        if args['--overview']:
            collection_overview(db_name=args['--db'], collection_name=args['--collection'])
        elif args['--key'] and args['--val']:
            get_doc(db_name=args['--db'], collection_name=args['--collection'], key=args['--key'], val=args['--val'],
                    subfield=args['--subfield'])

    ## None
    else:
        print("Missing argument!\n")
        print(__doc__.strip())
