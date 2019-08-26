#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# This file is part of the BHC Complexity Toolkit.
#
# The BHC Complexity Toolkit is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The BHC Complexity Toolkit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the BHC Complexity Toolkit.  If not, 
# see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------
# Copyright 2019, Mark D. Flood
#
# Author: Mark D. Flood
# Last revision: 22-Jun-2019
# -----------------------------------------------------------------------------

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='bhccpx',
    version='0.3.0',
    description='Package for analyzing the structure and complexity of BHCs',
    long_description=readme,
    author='Mark D. Flood',
    license=license,
    packages=find_packages(exclude=('tests', 'docs', 'data', 'out',)),
    install_requires=['graphviz>=0.11',
                      'networkx>=2.2',
                      'numpy>=1.16',
                      'pandas>=0.24',
                      'progressbar']
)
