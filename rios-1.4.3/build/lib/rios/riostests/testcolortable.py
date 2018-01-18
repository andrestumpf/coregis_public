"""
A simple test of color-table funcionality

Creates a raster image, and adds a color table to it. Then
reads it back, and checks that we get the same thing we
put there. 

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
import os
import numpy

from rios import rat

from . import riostestutils

TESTNAME = "TESTCOLORTABLE"

def run():
    """
    Run the color table test
    """
    riostestutils.reportStart(TESTNAME)
    
    filename = 'test.img'
    ds = riostestutils.genThematicFile(filename)
    del ds
    
    clrTbl = numpy.array([
        [0, 0, 0, 0, 0],
        [1, 255, 0, 0, 0],
        [2, 0, 255, 0, 0],
        [3, 0, 0, 255, 0],
    ], dtype=numpy.uint8)

    rat.setColorTable(filename, clrTbl)
    clrTbl2 = rat.getColorTable(filename)
    
    # Check the first N rows of the color table. This assumes
    # that we gave contiguous values in the input table (which we did)
    n = len(clrTbl)
    if (clrTbl == clrTbl2[:n]).all():
        riostestutils.report(TESTNAME, 'Passed')
        ok = True
    else:
        riostestutils.report(TESTNAME, 'Retrieved color table not equal to that written')
        ok = False

    os.remove(filename)

    return ok
