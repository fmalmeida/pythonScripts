.. _align2subsetgbk:

align2subsetgbk
===============

This script is a little more tricky. What does it do?

1. It takes a genbank file and a (Nucleotide) fasta file
2. It will align the fasta to the (contig) sequences in the genbank file
3. Then, it will filter the alignments based on the given thresholds and take the its coordinates
4. Finally, using these coordinates it will filter the genbank file and output only a subset of the annotation features that is found inside these alignments

CLI help message
----------------

.. literalinclude:: ./align2subsetgbk_help.txt
   :language: stdout

Usage
-----

The usage is super simple:

.. code-block:: none

   falmeida-py align2subsetgbk \
       --gbk sample.gbk \
       --fasta genomic_islands.fna \
       --extension 50 \
       --culling_limit 2
   
   Processing file: annotation/sample.gbk!
                           qseqid  qstart   qend   qlen       sseqid   sstart     send     slen  evalue  length  pident  gaps  gapopen  bitscore
   0  NC_016845.1_1287077_1326610       1  39533  39533  NC_016845.1  1287078  1326610  5333942     0.0   39533   100.0     0        0     73004
   1  NC_016845.1_1778390_1811349       1  32959  32959  NC_016845.1  1778391  1811349  5333942     0.0   32959   100.0     0        0     60864
   2  NC_016845.1_2286181_2305760       1  19579  19579  NC_016845.1  2286182  2305760  5333942     0.0   19579   100.0     0        0     36156
   3  NC_016845.1_4049987_4083106       1  33119  33119  NC_016845.1  4049988  4083106  5333942     0.0   33119   100.0     0        0     61160
   4  NC_016845.1_4821157_4858652       1  37495  37495  NC_016845.1  4821158  4858652  5333942     0.0   37495   100.0     0        0     69241


.. note::

   This command will align the genomic islands to the genbank to find the annotation features that are found inside the regions that align to these genomic islands.

   In the example, the alignments coordinates will be "extendend" 50nt in both directions and, the two best alignments will be used for each sequence (from fasta).