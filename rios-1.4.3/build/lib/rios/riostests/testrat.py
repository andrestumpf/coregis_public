#!/usr/bin/env python
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

"""
Test the rios.rat functionality
"""
import os
import sys
import numpy

from rios import rat
from rios import calcstats
from rios import cuiprogress
from rios.riostests import riostestutils

from rios import ratapplier
from rios import applier

TESTNAME = 'TESTRAT'

def run():
    """
    Run tests of the rios.rat functions
    """
    riostestutils.reportStart(TESTNAME)
    allOK = True

    imgfile = 'test.img'
    ratValues = makeTestFile(imgfile)
    nValues = len(ratValues)
    
    columnList = [
        ("Int32", numpy.int32),
        ("Float32", numpy.float32),
        ("Unicode", numpy.dtype('U10'))
    ]
    # Only test old string type for python 2
    if sys.version_info.major < 3:
        columnList.append(("String", numpy.dtype('S10')))
    
    allOk = True
    for (colName, arrayDtype) in columnList:
        # Write the array into the file, with the given datatype
        ratValues_type = ratValues.astype(arrayDtype)
        rat.writeColumn(imgfile, colName, ratValues_type)
        
        # Read it back, and check that the values are the same
        ratValues_fromFile = rat.readColumn(imgfile, colName)[:nValues].astype(ratValues.dtype)
        if not (ratValues_fromFile == ratValues).all():
            riostestutils.report(TESTNAME, "Value mis-match for column %s"%(colName))
            allOk = False
    
    if os.path.exists(imgfile):
        os.remove(imgfile)

    if allOK:
        riostestutils.report(TESTNAME, "Passed")

    return allOk



def makeTestFile(imgfile, withRat=True):
    # Make a test image with a simple RAT
    nRows = 100
    nCols = 1
    ds = riostestutils.createTestFile(imgfile, numRows=nRows, numCols=nCols)
    imgArray = numpy.ones((nRows, nCols), dtype=numpy.uint8)
    imgArray[1:10, 0] = numpy.arange(1, 10)
    imgArray[50:, 0] = 0
    band = ds.GetRasterBand(1)
    band.WriteArray(imgArray)

    band.SetMetadataItem('LAYER_TYPE', 'thematic')
    del ds
    
    # Note that the RAT has a row for lots of values which have no corresponding pixel
    ratValues = (numpy.mgrid[0:nRows]+10).astype(numpy.int32)
    ratValues[0] = 500

    return ratValues

    
if __name__ == "__main__":
    run()
