.. _tsv2markdown:

tsv2markdown
============

Often I found myself having to add the contents of TSV (or CSV) files in markdown documents such as jupyter notebooks or rmarkdown documents. The process was very tedious, so I create a super simple script to output the contents so it could just be pasted into these documents.

CLI help message
----------------

.. literalinclude:: ./tsv2markdown_help.txt
   :language: stdout

Usage
-----

The usage is super simple. For example, given the following CSV:

.. code-block:: none

   sample1,brazil,river,2020
   sample2,argentina,sewer,2020
   sample3,brazil,sewer,2020

We could rapidply create a markdown table entry with:

.. code-block:: bash

   falmeida-py tsv2markdown --csv sample.csv --header "sample id,country,collection source,year"

   | sample id   | country   | collection source   |   year |
   |-------------|-----------|---------------------|--------|
   | sample1     | brazil    | river               |   2020 |
   | sample2     | argentina | sewer               |   2020 |
   | sample3     | brazil    | sewer               |   2020 |

.. note::

   If the first line of the document (TSV or CSV) is a header, will to not need to use the ``--header`` option, because it will use the frist line as the header.