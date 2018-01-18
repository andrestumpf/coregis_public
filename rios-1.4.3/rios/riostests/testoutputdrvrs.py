"""
Test writing output with a range of different formats
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

from osgeo import gdal

from rios import applier
from rios.riostests import riostestutils

TESTNAME = 'TESTOUTPUTDRVRS'


def run():
    """
    Run tests of a number of more common output format drivers
    """
    riostestutils.reportStart(TESTNAME)
    
    usingExceptions = gdal.GetUseExceptions()
    gdal.UseExceptions()

    driverTestList = [
        ('HFA', ['COMPRESS=YES'], '.img'),
        ('GTiff', ['COMPRESS=LZW', 'TILED=YES', 'INTERLEAVE=BAND'], '.tif'),
        ('ENVI', ['INTERLEAVE=BSQ'], ''), 
        ('KEA', [], '.kea')
    ]
    # Remove any which current GDAL not suporting
    driverTestList = [(drvrName, options, suffix) for (drvrName, options, suffix) in driverTestList
        if gdal.GetDriverByName(drvrName) is not None]
    riostestutils.report(TESTNAME, 'Testing drivers {}'.format(str([d[0] for d in driverTestList])))


    filename = 'test.img'
    riostestutils.genRampImageFile(filename)
    
    ok = True
    outfileList = []
    errorList = []
    for (drvr, creationOptions, suffix) in driverTestList:
        infiles = applier.FilenameAssociations()
        outfiles = applier.FilenameAssociations()
        infiles.inimg = filename
        outfiles.outimg = "testout"+suffix
        outfileList.append(outfiles.outimg)
        
        controls = applier.ApplierControls()
        controls.setOutputDriverName(drvr)
        
        try:
            applier.apply(copyImg, infiles, outfiles, controls=controls)
        except Exception as e:
            ok = False
            errorList.append("{}:{}".format(drvr, str(e)))


    if ok:
        riostestutils.report(TESTNAME, "Passed")
    else:
        riostestutils.report(TESTNAME, "Resulted in these apparent errors:\n  {}".format('\n  '.join(errorList)))
    
    for fn in [filename] + outfileList:
        if os.path.exists(fn):
            drvr = gdal.Open(fn).GetDriver()
#            drvr.Delete(fn)
    
    if not usingExceptions:
        gdal.DontUseExceptions()
            
    return ok


def copyImg(info, inputs, outputs):
    """
    Copy the input data to the output file
    """
    outputs.outimg = inputs.inimg


if __name__ == "__main__":
    run()
