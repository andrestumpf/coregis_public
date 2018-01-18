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
A simple set of tests of the access to coordinates from with RIOS. 
There are a number of routines on the ReaderInfo class which
allow the RIOS use function to access coordinate systems of the 
working grid, and the underlying rasters, and this tests whether 
these behave as expected, in a range of circumstances. 

"""
from __future__ import division

import os
import math

from rios import applier

from . import riostestutils

TESTNAME = "TESTCOORDS"

# The correct answers, hard-wired here. These will need to be updated if the coordinates
# are changed in riostestutils.genRampImageFile(). 
XSTART = float(riostestutils.DEFAULT_XLEFT)
YSTART = float(riostestutils.DEFAULT_YTOP)
PIX = riostestutils.DEFAULT_PIXSIZE
STEP = PIX * applier.DEFAULTWINDOWXSIZE
HALFPIX = 0.5 * PIX
OLAP = 2
MARGIN = OLAP * PIX
OFFSET_PIXELS = 100
OFFSET_2ND_IMAGE = OFFSET_PIXELS * PIX

TSTPT_X = XSTART + OFFSET_2ND_IMAGE + 10.5 * PIX
TSTPT_Y = YSTART - OFFSET_2ND_IMAGE - 10.5 * PIX

TSTPIX_1FILE_NOOVERLAP = [(110, 110)]
TSTPIX_1FILE_OVERLAP2 = [(112, 112)]
TSTPIX_2FILE_OVERLAP2 = [(12, 12)]
NUMBLOCKSX = int(math.ceil(riostestutils.DEFAULT_COLS /  applier.DEFAULTWINDOWXSIZE))
NUMBLOCKSY = int(math.ceil(riostestutils.DEFAULT_ROWS /  applier.DEFAULTWINDOWYSIZE))
CENTRES_1FILE_NOOVERLAP = [(XSTART + j * STEP + HALFPIX, YSTART - i * STEP - HALFPIX) 
    for i in range(NUMBLOCKSY) for j in range(NUMBLOCKSX)]
CENTRES_1FILE_OVERLAP2 = [(XSTART + j * STEP + HALFPIX - MARGIN, YSTART - i * STEP - HALFPIX + MARGIN) 
    for i in range(NUMBLOCKSY) for j in range(NUMBLOCKSX)]
NUMBLOCKSX_WITHOFFSET = int(math.ceil((riostestutils.DEFAULT_COLS - 2 * OFFSET_PIXELS)/  applier.DEFAULTWINDOWXSIZE))
NUMBLOCKSY_WITHOFFSET = int(math.ceil((riostestutils.DEFAULT_ROWS - 2 * OFFSET_PIXELS)/  applier.DEFAULTWINDOWYSIZE))
CENTRES_2FILE_OVERLAP2 = [(XSTART + OFFSET_2ND_IMAGE+ j * STEP + HALFPIX - MARGIN, 
    YSTART - OFFSET_2ND_IMAGE- i * STEP - HALFPIX + MARGIN) 
    for i in range(NUMBLOCKSY_WITHOFFSET) for j in range(NUMBLOCKSX_WITHOFFSET)]


def run():
    """
    Run the test
    """
    allOK = True
    
    riostestutils.reportStart(TESTNAME)

    ramp1 = 'ramp1.img'
    ramp2 = 'ramp2.img'
    riostestutils.genRampImageFile(ramp1)

    # Second file is same as first, but shifted 100 pixels right and down. 
    xLeft = riostestutils.DEFAULT_XLEFT + OFFSET_2ND_IMAGE
    yTop = riostestutils.DEFAULT_YTOP - OFFSET_2ND_IMAGE
    riostestutils.genRampImageFile(ramp2, xLeft=xLeft, yTop=yTop)
    
    # Set up some RIOS calls
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    controls = applier.ApplierControls()

    # Simplest check. Can we get back the coordinates from a single input file.     
    infiles.img = ramp1
    otherargs.tstPixList = []
    otherargs.centresList = []
    applier.apply(getCoords, infiles, outfiles, otherargs, controls=controls)
    ok = checkCoords('1 file, overlap=0', otherargs, TSTPIX_1FILE_NOOVERLAP, CENTRES_1FILE_NOOVERLAP)
    if not ok: allOK = False

    # Single input file, with an overlap set
    otherargs.tstPixList = []
    otherargs.centresList = []
    controls.setOverlap(OLAP)
    applier.apply(getCoords, infiles, outfiles, otherargs, controls=controls)
    ok = checkCoords('1 file, overlap=2', otherargs, TSTPIX_1FILE_OVERLAP2, CENTRES_1FILE_OVERLAP2)
    if not ok: allOK = False

    # Two input files, with an overlap set
    infiles.img = [ramp1, ramp2]
    otherargs.tstPixList = []
    otherargs.centresList = []
    controls.setOverlap(OLAP)
    applier.apply(getCoords, infiles, outfiles, otherargs, controls=controls)
    ok = checkCoords('2 files, overlap=2', otherargs, TSTPIX_2FILE_OVERLAP2, CENTRES_2FILE_OVERLAP2)
    if not ok: allOK = False
    
    # Clean up
    for filename in [ramp1, ramp2]:
        os.remove(filename)
    
    if allOK:
        riostestutils.report(TESTNAME, "Passed")

    
    return allOK


def getCoords(info, inputs, outputs, otherargs):
    """
    Accumulate some lists of the coordinates of the top-left corner pixel
    of each block, expressed in a variety of ways
    """
    (x, y) = info.getBlockCoordArrays()
    otherargs.centresList.append((x[0,0], y[0,0]))
    
    (r, c) = info.getPixRowColBlock(TSTPT_X, TSTPT_Y)
    if r is not None and c is not None:
        otherargs.tstPixList.append((r, c))


def checkCoords(testConditionStr, otherargs, tstPixList, centresList):
    """
    Check that the coordinate lists created in otherargs are the same as 
    the correct answers passed in to this routine. 
    """
    ok = checkCoordList(otherargs.tstPixList, tstPixList)
    if not ok:
        riostestutils.report(TESTNAME, '%s: Single pixel row/col mis-match. %s, %s' % (testConditionStr,
            otherargs.tstPixList, tstPixList))
    
    ok = checkCoordList(otherargs.centresList, centresList)
    if not ok:
        riostestutils.report(TESTNAME, '%s: Centre coordinates mis-match. %s, %s' % (testConditionStr,
            otherargs.centresList, centresList))
    
    return ok


def checkCoordList(list1, list2):
    """
    Check equality of given lists of (x,y) coords. Return True if they
    are the same, False otherwise. 
    """
    ok = True
    if len(list1) != len(list2):
        ok = False
    if ok:
        for i in range(len(list1)):
            if list1[i] != list2[i]:
                ok = False
    return ok
