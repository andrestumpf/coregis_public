#!/usr/bin/env python
"""
Test the calculation of Raster Attribute Table statistics by rios.fileinfo
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

import numpy

from rios import calcstats
from rios import fileinfo
from rios import rat
from rios import cuiprogress
from . import riostestutils

TESTNAME = 'TESTRATSTATS'

def run():
    """
    Run a test of RAT statistics calculation
    """
    riostestutils.reportStart(TESTNAME)
    
    allOK = True

    imgfile = 'test.img'
    nRows = 100
    nCols = 1
    ds = riostestutils.createTestFile(imgfile, numRows=nRows, numCols=nCols)
    imgArray = numpy.ones((nRows, nCols), dtype=numpy.uint8)
    imgArray[1:10, 0] = numpy.arange(1, 10)
    imgArray[50:, 0] = 0
    band = ds.GetRasterBand(1)
    band.WriteArray(imgArray)
    
    nullDN = 0
    calcstats.calcStats(ds, ignore=nullDN, progress=cuiprogress.SilentProgress())
    columnName = 'Value'
    # Note that the RAT has a row for lots of values which have no corresponding pixel
    ratValues = (numpy.mgrid[0:nRows]**2).astype(numpy.float64)
    ratValues[0] = 500
    rat.writeColumnToBand(band, columnName, ratValues)
    band.SetMetadataItem('LAYER_TYPE', 'thematic')
    del ds
    
    ratStats = fileinfo.RatStats(imgfile, columnlist=[columnName])
    
    # Construct an image of the values, by translating pixel values into RAT values
    ratValImg = numpy.zeros(imgArray.shape, dtype=numpy.float64)
    for dn in numpy.unique(imgArray):
        if dn != 0:
            mask = (imgArray == dn)
            ratValImg[mask] = ratValues[dn]
    ratValImgNonNull = ratValImg[imgArray!=0]

    # Now find the "true" values of the various stats for this image (i.e. this is the
    # histogramweighted=True case, which I think will be the most common one)
    trueMean = ratValImgNonNull.mean()
    trueStddev = ratValImgNonNull.std()
    trueMin = ratValImgNonNull.min()
    trueMax = ratValImgNonNull.max()
    
    tolerance = 0.000001
    if not equalTol(ratStats.Value.mean, trueMean, tolerance):
        riostestutils.report(TESTNAME, "Mismatched means: %s, %s" % 
            (repr(ratStats.Value.mean), repr(trueMean)))
        allOK = False
    if not equalTol(ratStats.Value.stddev, trueStddev, tolerance):
        riostestutils.report(TESTNAME, "Mismatched stddevs: %s, %s" % 
            (repr(ratStats.Value.stddev), repr(trueStddev)))
        allOK = False
    if not equalTol(ratStats.Value.min, trueMin, tolerance):
        riostestutils.report(TESTNAME, "Mismatched mins: %s, %s" % 
            (repr(ratStats.Value.min), repr(trueMin)))
        allOK = False
    if not equalTol(ratStats.Value.max, trueMax, tolerance):
        riostestutils.report(TESTNAME, "Mismatched maxes: %s, %s" % 
            (repr(ratStats.Value.max), repr(trueMax)))
        allOK = False

    if allOK:
        riostestutils.report(TESTNAME, "Passed")
    
    return allOK


def equalTol(a, b, tol):
    """
    Compare two values to within a tolerance. If the difference
    between the two values is smaller than the tolerance, 
    then return True
    """
    diff = abs(a - b)
    return (diff < tol)


if __name__ == "__main__":
    run()

