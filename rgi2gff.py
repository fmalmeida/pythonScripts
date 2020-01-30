#!/usr/bin/env python3

'''rgi2gff.py last modified 2018-07-09

rgi2gff.py -b RGI_out.txt > output.gff3

RGI tabular output must be first parsed with:

sed -i 's/ # /#/g' and sed -i 's/ /_/g'
'''

#
import sys
import argparse
import time
from collections import defaultdict
#

# Returns s truncated at the n'th (3rd by default) occurrence of the delimiter, d.
def trunc_at(s, d, n=3):
    return d.join(s.split(d, n)[:n])


def write_line(outlist, wayout):
    outline = "\t".join(outlist)
    print(outline, file=wayout)


def main(argv, wayout):
    if not len(argv):
        argv.append("-h")

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument('-f', '--file', help="RGI txt output file")
    args = parser.parse_args(argv)

    # counter for number of lines, and strand flips
    linecounter, writecounter = 0, 0
    plusstrand, minusstrand = 0, 0

    hitDictCounter = defaultdict(int)
    print("Starting RGI parsing on %s" % (args.file), time.asctime(), file=sys.stderr)
    for line in open(args.file, 'r'):
        linecounter += 1
        ORF_ID, Contig, Start, Stop, Orientation, Cut_Off, Pass_Bitscore, Best_Hit_Bitscore, \
        Best_Hit_ARO, Best_Identities, ARO, Model_type, SNPs_in_Best_Hit_ARO, Other_SNPs, \
        Drug_Class, Resistance_Mechanism, AMR_Gene_Family, Predicted_DNA, Predicted_Protein, \
        CARD_Protein_Sequence, Percentage_Length_Reference_Sequence, ID, Model_id= line.rstrip().split("\t")

        # Reformat Contig ID
        Contig = trunc_at(Contig, "_", n=2)
        # Reformat Drug_Classes
        Drug_Class = Drug_Class.replace(";_", "&")
        # Reformat Resistance_Mechanism
        Resistance_Mechanism = Resistance_Mechanism.replace(";_", "&")
        # Define attributes
        attributes = "Additional_database=CARD_RGI;CARD_gene_name={0};CARD_gene_family={1};CARD_ARO={2};CARD_target_drugs={3};CARD_resistance_mechanism={4}".format(
        Best_Hit_ARO, AMR_Gene_Family, ARO, Drug_Class, Resistance_Mechanism)

        #print("Starting base on {0}".format(Start), file=sys.stderr)

        # convert strings of start and end to integers for calculations
        iqend = int(Start)
        iqstart = int(Stop)
        # as start must always be less or equal to end, reverse them for opposite strand hits
        if iqstart <= iqend:
            strand = "+"
            outlist = [Contig, "CARD_RGI", "resistance", Start,
                       Stop, Best_Hit_Bitscore	, strand, ".", attributes]
            plusstrand += 1
        else:
            strand = "-"
            outlist = [Contig, "CARD_RGI", "resistance", Start,
                       Stop, Best_Hit_Bitscore, strand, ".", attributes]
            minusstrand += 1

        writecounter += 1

        # Print results
        write_line(outlist, wayout)

    # Logs
    print("Parsed %d lines" % (linecounter), time.asctime(), file=sys.stderr)
    print("Found %d forward and %d reverse hits" % (
        plusstrand, minusstrand), time.asctime(), file=sys.stderr)
    print("Wrote %d matches" % (writecounter), time.asctime(), file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1:], sys.stdout)
