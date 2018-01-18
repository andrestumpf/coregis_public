#!/usr/bin/env python
"""
Main program for RIOS MPI subprocesses. 
This uses the MPI calls to receive and send data
from the main process.

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

from mpi4py import MPI
from rios.parallel import subproc

import sys
from io import BytesIO

if __name__ == "__main__":

    # get the handle to the parent
    comm = MPI.Comm.Get_parent()

    # keep going until we are told to exit
    while True:
        # get the data from the parent
        data = comm.recv(source=0)

        status, blockData = data
        if not status:
            break

        # wrap it with a BytesIO object so it can be treated
        # like a file by subproc.runJob() - this keeps compatibility
        # with the other sub job types
        inf = BytesIO(blockData)
        inf.seek(0)

        # output file object
        outf = BytesIO()

        # do the processing
        subproc.runJob(inf, outf)

        # send the result back
        outdata = outf.getvalue()
        comm.send(outdata, dest=0)

