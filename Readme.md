# fa-py (Felipe Almeida python scripts)

This repository has been turned into an installable python package in order to facilitate the distribution of my custom python scripts and, in turn, make them easier to execute.

## Installation

Installation is super easy and perhaps not required:

```bash
# Download
git clone https://github.com/fmalmeida/pythonScripts.git
cd pythonScripts

# Run without installing
python3 fa-py-runner.py -h

# Or install with pip and run anywhere
pip install setup.py
fa-py -h

# help
fa-py: a package to the simple distribution of my custom scripts.

Copyright (C) 2020 Felipe Marques de Almeida (almeidafmarques@gmail.com)
License: Public Domain

usage:
    fa-py [ -h|--help ] [ -v|--version ] [ --license ]
    fa-py <command> [ -h|--help ] [ <args>... ]

options:
    -h --help                                               Show this screen
    -v --version                                            Show version information
    --license                                               Show LEGAL LICENSE information

commands:
    tsv2markdown                                            Command for rapid convertion of tsv or csv to markdown tables.
    splitgbk                                                Command to split multisequence genbank files into individual files.
    align2subsetgbk                                               Command to subset genbank files based on alignments to a FASTA file.
    blasts                                                  Command to execute automatized blast commands.

Use: `fa-py <commmand> -h` to get more help and see examples.
```

## Old scripts

All my python scripts as single scripts, that may or may not be included in the package are available in the `bkp` folder!

## License

This repository has no warranty and is free to use, modify and share under GNU GENERAL PUBLIC LICENSE version 3.

## Contact

Felipe Almeida <almeidafmarques@gmail.com>
