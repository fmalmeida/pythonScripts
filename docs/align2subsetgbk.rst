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


.. note::

   This command will align the genomic islands to the genbank to find the annotation features that are (possibly) inside the genomic islands.

   In the example, the alignments coordinates will be "extendend" 50nt in both directions and, the two best alignments will be used for each sequence (from fasta).