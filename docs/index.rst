.. _index:

falmeida-py (Felipe Almeida python scripts)
===========================================

This is just a simple repository of python scripts that I use in a daily basis. Some of them are used in my projects, some are just for fun. Thus, in order to provide a better way to use my scripts throughout machines, I decided to create a package with them.

Installation
------------

Conda package
"""""""""""""

Users are advised to create separate conda environments

.. code-block:: bash

   # Get the conda package
   mamba create -n falmeida-py -c anaconda -c conda-forge -c bioconda -c falmeida falmeida-py

.. tip::
   
   Users are advised to use `mamba <https://github.com/mamba-org/mamba>`_ since it is faster.

Available features
------------------

Available features can be visualised with the command line help message:

.. literalinclude:: ./help_message.txt
   :language: stdout


.. note::

   All the commands will have a page describing its usage.

.. toctree::
   :hidden:
   :caption: Home

   self

.. toctree::
   :hidden:
   :caption: Reference book

   tsv2markdown
   splitgbk
