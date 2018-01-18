"""
Helper function for unpickling a user function and its
parameters and running it. 

Normally called from rios_subproc.py, but some parallel
methods (MPI, multiprocessing etc) manage their own forking
of the Python process so we must have this available as an 
importable function.
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

import os
try:
    import cPickle as pickle
except ImportError:
    import pickle

def runJob(inf, outf, inFileName=None):
    """
    Run one job reading the pickle from inf
    and writing the output to outf.

    if inFileName is not None then then it
    will be deleted once it is read from.
    """

    # Read the pickled input
    (fn, jobInfo) = pickle.load(inf)
    
    # If using a disk file for input, close it and remove it now that we have read it
    if inFileName is not None:
        inf.close()
        os.remove(inFileName)

    params = jobInfo.getFunctionParams()   
    
    # Execute the function, with the given input data
    fn(*params)

    outputs = jobInfo.getFunctionResult(params)
    
    # Pickle and write out the output
    pickle.dump(outputs, outf)
