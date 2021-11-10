.. _splitgbk:

splitgbk
========

This is a super sample script that performs only one task:

* Given a genbank file with multiple sequences, it splits the file into multiple files, each with one sequence and its related annotation features features.

CLI help message
----------------

.. literalinclude:: ./splitgbk_help.txt
   :language: stdout

Usage
-----

The usage is super simple:

.. code-block:: none

   falmeida-py tsv2markdown --gbk multi_contig.gbk --outdir single_contig_gbks