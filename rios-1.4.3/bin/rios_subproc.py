#!/usr/bin/env python
"""
Main program for RIOS subprocesses. 

"""
# This file is part of RIOS - Raster I/O Simplification
# Copyright (C) 2012  Sam Gillingham, Neil Flood
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import print_function

from rios.parallel import subproc

import sys

if __name__ == "__main__":
    nArgs = len(sys.argv) - 1

    # This is for the case of no commandline args    
    if sys.version_info[0] > 2:
        # For Python 3, use binary buffer objects
        # otherwise unpickling fails
        inf = sys.stdin.buffer
        outf = sys.stdout.buffer
    else:
        inf = sys.stdin
        outf = sys.stdout
    
    # These are for the cases with some commandline args
    inFileName = None
    if nArgs >= 1:
        inFileName = sys.argv[1]
        inf = open(inFileName, 'rb')
    if nArgs == 2:
        outf = open(sys.argv[2], 'wb')

    subproc.runJob(inf, outf, inFileName)
