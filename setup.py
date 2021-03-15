#!/usr/bin/env python3
"""
This is the fa-py installation script. Assuming you're in the same directory, it can be run
like this: `python3 setup.py install`, or (probably better) like this: `pip3 install .`

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

from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


# Get the program version from another file.
__version__ = '0.0.0'
exec(open('fapy/version.py').read())

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(name='fa-py',
      version=__version__,
      description='fa-py: a package to the simple distribution of my custom scripts.',
      long_description=readme(),
      long_description_content_type='text/markdown',
      url='https://github.com/fmalmeida/pythonScripts',
      author='Felipe Almeida',
      author_email='almeidafmarques@gmail.com',
      license='GPLv3',
      packages=['fapy'],
      install_requires=required,
      entry_points={"console_scripts": ['fa-py = fapy.__main__:main']},
      include_package_data=True,
      zip_safe=False,
      python_requires='>=3.6')
