#!/usr/bin/env python3
# coding: utf-8

## Def help message
"""
Ask mongo: This command allows you to ask/interrogate the mongo dbs in your machine

usage:
    ask-mongo.py [ -h|--help ] [ --list_dbs ]
    ask-mongo.py [ --db <db_name> --list_collections ]
    ask-mongo.py [ --db <db_name> --collection <collection_name> ] [ --overview ]

options:

    -h, --help                                               Show this screen
    --list_dbs                                               List available mongo dbs in your system
    --db=<db_name>                                           Mongo DB to be queryed
    --list_collections                                       Check the available collections in a given database
    --collection=<collection_name>                           Collection name (from mongo db) to be queryed
    --overview                                               Documents in a given collection (from a mongo db)

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

#########################
### Check your mongos ###
#########################
def check_mongos():

    # Start message
    print("""
    The available mondo dbs found in your system are:\n
    """)

    # Print databases
    os.system('mongo --quiet --eval  "printjson(db.adminCommand(\'listDatabases\'))"')

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

################
### Def main ###
################
if __name__ == '__main__':
    args = docopt(__doc__, version='v1.0 by Felipe Marques de Almeida')

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

    ## None
    else:
        print("Missing argument!\n")
        print(__doc__.strip())
