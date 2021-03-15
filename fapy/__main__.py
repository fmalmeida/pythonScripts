#!/usr/bin/env python3
license="""
Copyright 2020 Felipe Almeida (almeidafmarques@gmail.com)
https://github.com/fmalmeida/pythonScripts

This file is part of my custom python scripts (fa-py) package, which is free: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. This package is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with fa-py package.
If not, see <http://www.gnu.org/licenses/>.
"""

## Def main help
usage="""
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
    include

Use: `fa-py <commmand> -h` to get more help and see examples.
"""

##################################
### Loading Necessary Packages ###
##################################
from docopt import docopt
from subprocess import call
import sys

########################
### Import functions ###
########################
from .version import *

## Defining main
def main():
    # Parse docopt
    __version__ = get_version()
    arguments = docopt(usage, version=__version__, help=False, options_first=True)

    ############################
    ### GFF overview command ###
    ############################
    if arguments['<command>'] == 'overview':

        # Parse docopt
        args_overview = docopt(usage_overview, version=__version__, help=False)

        if args_overview['overview'] and args_overview['--help']:
            print(usage_overview.strip())

        elif args_overview['overview'] and args_overview['--input']:
            overview(args_overview['--input'])

        else:
            print(usage_overview.strip())

    #####################
    ### Check license ###
    #####################
    elif arguments['--license']:
        print(license.strip())

    #######################################
    ### Without commands nor parameters ###
    #######################################
    else:
        print(usage.strip())

## Calling main
if __name__ == '__main__':
    main()
