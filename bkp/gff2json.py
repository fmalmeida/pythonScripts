#!/usr/bin/env python3

'''gff2json.py last modified 2020-03-29

gff2json.py -i [input.gff] -o [output.json]

A pipeline to convert GFF files into JSON format that can be
used in MongoDB servers.
'''

#
import sys
import argparse
import time
from collections import namedtuple
import gzip
import urllib.request, urllib.parse, urllib.error
import json
#

# Initialized GeneInfo named tuple. Note: namedtuple is immutable
gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

# Function to parse Attributes column
def parseGFFAttributes(attributeString):
    """Parse the GFF3 attribute column and return a dict"""
    if attributeString == ".": return {}
    ret = {}
    for attribute in attributeString.split(";"):
        try:
            key, value = attribute.split("=")
            ret[key] = value.strip()
        except ValueError:
            print(attribute)
            break

    return ret

# Function to parse all file
def parseGFF3(filename, output):
    """
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.

    Supports transparent gzip decompression.
    """
    gff_dict = {}
    final    = []
    #Parse with transparent decompression
    openFunc = gzip.open if filename.endswith(".gz") else open
    with openFunc(filename) as infile:
        for line in infile:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            #If this fails, the file format is not standard-compatible
            assert len(parts) == len(gffInfoFields)
            #Separate Values
            seq        = parts[0]
            source     = parts[1]
            feature    = parts[2]
            start      = parts[3]
            end        = parts[4]
            score      = parts[5]
            strand     = parts[6]
            phase      = parts[7]
            attributes = parseGFFAttributes(parts[8])
            gff_dict = {
                "CDS" : attributes['ID'],
                "seqid"      : seq,
                "source"    : source,
                "type"      : feature,
                "start"     : start,
                "end"       : end,
                "score"     : score,
                "strand"    : strand,
                "phase"     : phase,
                "attributes": attributes
            }
            final.append(gff_dict)

    with open(output, 'a+') as fp:
        json.dump(final, fp, indent=2)


def main(argv):
    if not len(argv):
        argv.append("-h")

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument('-i', '--input',  help="GFF file to be converted into JSON")
    parser.add_argument('-o', '--output', help="JSON file output name")
    args = parser.parse_args(argv)

    # Execute the program
    parseGFF3(args.input, args.output)

    # Print Goodbye
    print("Finished!")


if __name__ == "__main__":
    main(sys.argv[1:])
