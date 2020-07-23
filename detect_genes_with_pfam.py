#!/usr/bin/env python3
# coding: utf-8

## Def help message
"""
A script to detect the presence of pfam domains of interest in target protein sequences
Dependencies:
    * It requires the HMMER package to be available in $PATH

---
Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain
Usage:
    ./detect_genes_with_pfam.py
    ./detect_genes_with_pfam.py -h|--help
    ./detect_genes_with_pfam.py -v|--version
    ./detect_genes_with_pfam.py pfam-index  [--pfam <Pfam-A.hmm>]
    ./detect_genes_with_pfam.py pfam-detect [--prots <pep fasta> --prefix <prefix> --pfam <Pfam-A.hmm> --pfam_list <target pfam ids>]

Options:
    -h --help                                     Show this screen.
    -v --version                                  Show version information
    pfam-index                                    Indexes pfam database with hmmpress if not done yet
    pfam-detect                                   Executes the detection pipeline
    --prots=<pep fasta>                           Path to the target protein sequences
    --prefix=<prefix>                             Prefix for writing output files [default: out]
    --pfam=<Pfam-A.hmm>                           Path to pfam database (hmm format)
    --pfam_list=<target pfam ids>                 Path to the file containing the list of pfam domains of interest. It accepts the pfam accession and domain name. It is preferred to
                                                  use the domain names, since some pfam accessions in the hmm file have extra decimal hex values that may conflit with the information
                                                  available on the webpage.
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt
import os

####################################
### Function for Pfam indexation ###
####################################
def pfam_index(pfamHMM):
    os.system("hmmpress {0}".format(pfamHMM))

###################################
### Function for Pfam detection ###
###################################
def pfam_detection(pfamHMM, pfamLIST, pepFASTA, prefix):
    os.system("hmmfetch -f {0} {1} | hmmsearch -o pfam.search.tmp --noali --tblout {3}_pfam_hits.txt - {2}".format(pfamHMM, pfamLIST, pepFASTA, prefix))
    os.system("grep -v \"#\" {0}_pfam_hits.txt | cut -d \"-\" -f 1 | sed 's/ *$//g' > tmp.target.list".format(prefix))
    os.system("awk \'NR==FNR{{ids[$0]; next}} ($1 in ids){{ printf \">\" $0 }}\' tmp.target.list RS='>' {0} > {1}_target.fa".format(pepFASTA, prefix))
    os.system("rm tmp.target.list pfam.search.tmp")



############
### Main ###
############
if __name__ == '__main__':
    arguments = docopt(__doc__, version='v1.0 by Felipe Marques de Almeida')

    ## Pfam indexation
    if arguments['pfam-index'] and arguments['--pfam']:
        print("Initiating Pfam hmm indexation")
        print("\nCommand used:\n\thmmpress {0}".format(arguments['--pfam']))
        print("\nOutput:")
        pfam_index(pfamHMM=arguments['--pfam'])

    # Pfam detection
    elif arguments['pfam-detect'] and arguments['--pfam'] and arguments['--prots'] and arguments['--pfam_list']:
        print("Initiating Pfam domain detection")
        print("\nCommands used:")
        print("\thmmfetch -f {0} {1} | hmmsearch -o pfam.search.tmp --noali --tblout pfam_hits.txt - {2}".format(arguments['--pfam'], arguments['--pfam_list'], arguments['--prots']))
        print("\tgrep -v \"#\" {0}_pfam_hits.txt | cut -d \"-\" -f 1 | sed 's/ *$//g' > tmp.target.list".format(arguments['--prefix']))
        print("\tawk \'NR==FNR{{ids[$0]; next}} ($1 in ids){{ printf \">\" $0 }}\' tmp.target.list RS='>' {0} > {1}".format(arguments['--prots'], prefix=arguments['--prefix']))
        print("\trm tmp.target.list pfam.search.tmp")
        print("\nGenes that contain the pfam domains of interest are in:\n\t{0}".format(prefix=arguments['--prefix']))
        pfam_detection(pfamHMM=arguments['--pfam'], pfamLIST=arguments['--pfam_list'], pepFASTA=arguments['--prots'], prefix=arguments['--prefix'])

    ## None
    else:
        print("Missing mandatory arguments")
        print("Please, check out the help message")
        print("")
        print(arguments)
