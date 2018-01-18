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
Test the ratapplier functionality
"""
import os
import numpy

from rios import rat
from rios import calcstats
from rios import cuiprogress
from rios.riostests import riostestutils

from rios import ratapplier
from rios import applier

TESTNAME = 'TESTRATAPPLIER'

def run():
    """
    Run tests of the ratapplier
    """
    riostestutils.reportStart(TESTNAME)
    allOK = True

    imgfile = 'test.img'
    makeTestFile(imgfile)
    
    ok = testOutputSameFile(imgfile)
    if not ok: allOK = False
    
    imgfile2 = "test2.img"
    ok = testDifferentOutput(imgfile, imgfile2)
    if not ok: allOK = False
    
    imgfile3 = "test3.img"
    ok = testReduceRat(imgfile, imgfile3)
    if not ok: allOK = False
    
    imgfile4 = "test4.img"
    ok = testNewRat(imgfile4)
    if not ok: allOK = False
    
    for tmpfile in [imgfile, imgfile2, imgfile3, imgfile4]:
        os.remove(tmpfile)

    if allOK:
        riostestutils.report(TESTNAME, "Passed")

    return allOK
    
    

def testOutputSameFile(imgfile):
    # Now test the ratapplier
    inRats = ratapplier.RatAssociations()
    outRats = ratapplier.RatAssociations()
    controls = ratapplier.RatApplierControls()
    
    inRats.img = ratapplier.RatHandle(imgfile)
    outRats.img = inRats.img
    controls.setBlockLength(5)
    
    ratapplier.apply(myFunc, inRats, outRats, controls=controls)
    
    col = rat.readColumn(imgfile, 'Value')
    colSqrd = rat.readColumn(imgfile, 'sqrd')
    ok = True
    if (col**2 != colSqrd).any():
        riostestutils.report(TESTNAME, "sqrd incorrect, in sameFile output")
        ok = False
    return ok

def myFunc(info, inputs, outputs):
    outputs.img.sqrd = inputs.img.Value ** 2


def testDifferentOutput(imgfile, imgfile2):
    makeTestFile(imgfile2, withRat=False)

    inRats = ratapplier.RatAssociations()
    outRats = ratapplier.RatAssociations()
    controls = ratapplier.RatApplierControls()
    
    inRats.img = ratapplier.RatHandle(imgfile)
    outRats.outimg = ratapplier.RatHandle(imgfile2)
    controls.setBlockLength(3)
    
    ratapplier.apply(myFuncDiffFile, inRats, outRats, controls=controls)
    
    col = rat.readColumn(imgfile, 'Value')
    colSqrd = rat.readColumn(imgfile2, 'sqrd')
    ok = True
    if (col**2 != colSqrd).any():
        riostestutils.report(TESTNAME, "sqrd incorrect, in differentFile output")
        ok = False
    return ok

def myFuncDiffFile(info, inputs, outputs):
    outputs.outimg.sqrd = inputs.img.Value ** 2


def testReduceRat(imgfile, imgfile3):
    """
    This test creates a new output image, with all odd pixel values 
    replaced with the even number above it. The RAT must then be copied across
    with the same reduction performed. In this case, only the even numbered 
    rows are written
    """
    # First we copy the raster, with the reduction of pixel values
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    infiles.inimg = imgfile
    outfiles.outimg = imgfile3
    # Make sure we use a format which actually supports RAT's
    controls = applier.ApplierControls()
    controls.setOutputDriverName('HFA')
    applier.apply(rasterReduceFunc, infiles, outfiles, controls=controls)
    
    # Now use ratapplier to reduce the RAT
    inRats = ratapplier.RatAssociations()
    outRats = ratapplier.RatAssociations()
    controls = ratapplier.RatApplierControls()
    
    inRats.img = ratapplier.RatHandle(imgfile)
    outRats.outimg = ratapplier.RatHandle(imgfile3)
    controls.setBlockLength(3)
    
    ratapplier.apply(ratReduceFunc, inRats, outRats, controls=controls)
    
    col = rat.readColumn(imgfile, 'Value')
    colEven = col[::2]
    colReduced = rat.readColumn(imgfile3, 'Value')[:len(colEven)]
    ok = True
    if (colEven != colReduced).any():
        riostestutils.report(TESTNAME, "Reduced RAT incorrect: %s, %s"%(colEven, colReduced))
        ok = False
    return ok

def rasterReduceFunc(info, inputs, outputs):
    outputs.outimg = inputs.inimg // 2

def ratReduceFunc(info, inputs, outputs):
    evenMask = (((info.inputRowNumbers // 2) * 2) == info.inputRowNumbers)
    outputs.outimg.Value = inputs.img.Value[evenMask]

def testNewRat(imgfile4):
    makeTestFile(imgfile4, withRat=False)

    inRats = ratapplier.RatAssociations()
    outRats = ratapplier.RatAssociations()
    controls = ratapplier.RatApplierControls()
    controls.setRowCount(256)
    
    outRats.outimg = ratapplier.RatHandle(imgfile4)
    controls.setBlockLength(3)
    
    ratapplier.apply(myFuncNewRat, inRats, outRats, controls=controls)
    
    col = rat.readColumn(imgfile4, 'newCol')
    colIntended = numpy.arange(256, dtype=numpy.uint32)
    ok = (col == colIntended).all()
    if not ok:
        riostestutils.report(TESTNAME, "New RAT incorrect: %s, %s"%(col, colIntended))
    return ok

def myFuncNewRat(info, inputs, outputs):
    outputs.outimg.newCol = numpy.arange(info.blockLen, dtype=numpy.uint32) + info.startrow

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

    nullDN = 0
    calcstats.calcStats(ds, ignore=nullDN, progress=cuiprogress.SilentProgress())
    columnName = 'Value'
    # Note that the RAT has a row for lots of values which have no corresponding pixel
    ratValues = (numpy.mgrid[0:nRows]+10).astype(numpy.int32)
    ratValues[0] = 500
    if withRat:
        rat.writeColumnToBand(band, columnName, ratValues)
    band.SetMetadataItem('LAYER_TYPE', 'thematic')
    del ds

    
if __name__ == "__main__":
    run()
