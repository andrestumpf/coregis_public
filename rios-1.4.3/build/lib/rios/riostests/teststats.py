"""
Test the calculation of statistics by rios.calcstats. 
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
from osgeo import gdal

from rios import calcstats, cuiprogress

from . import riostestutils

TESTNAME = 'TESTSTATS'

def run():
    """
    Run a test of statistics calculation
    """
    riostestutils.reportStart(TESTNAME)
    
    nullVal = 0
    
    # We repeat the basic test for a number of different GDAL datatypes, with different
    # ranges of data. Each element of the following list is a tuple of
    #    (gdalDataType, numpyDataType, scalefactor)
    # for which the test is run. The original data being scaled is in 
    # the range 25-100 (after clobbering half the array as nulls, to ensure that
    # the nulls are enough to make a difference). 
    dataTypesList = [
        (gdal.GDT_Byte, numpy.uint8, 1),
        (gdal.GDT_UInt16, numpy.uint16, 1),
        (gdal.GDT_Int16, numpy.int16, 300),
        (gdal.GDT_UInt16, numpy.uint16, 300),
        (gdal.GDT_Int32, numpy.int32, 30000),
        (gdal.GDT_UInt32, numpy.uint32, 30000),
        (gdal.GDT_Float32, numpy.float32, 1), 
        (gdal.GDT_Float32, numpy.float32, 100), 
        (gdal.GDT_Float32, numpy.float32, 0.01) 
    ]
    
    # We repeat these tests on a number of different drivers, if they are available,
    # as some stats-related things may work fine on some drivers but not on others. 
    driverTestList = [
        ('HFA', ['COMPRESS=YES']),
        ('GTiff', ['COMPRESS=LZW', 'TILED=YES', 'INTERLEAVE=BAND']),
        ('KEA', [])
    ]
    # Remove any which current GDAL not suporting
    driverTestList = [(drvrName, options) for (drvrName, options) in driverTestList
        if gdal.GetDriverByName(drvrName) is not None]
    
    # Loop over all drivers
    for (driverName, creationOptions) in driverTestList:
        # Loop over all datatype tuples in the list
        for (fileDtype, arrDtype, scalefactor) in dataTypesList:
            imgfile = 'test.img'
            ds = riostestutils.createTestFile(imgfile, dtype=fileDtype, driverName=driverName, 
                    creationOptions=creationOptions)
            rampArr = riostestutils.genRampArray().astype(arrDtype) * scalefactor
            (nRows, nCols) = rampArr.shape
            # Set half of it to null
            rampArr[:, :nCols//2] = nullVal
            band = ds.GetRasterBand(1)
            band.WriteArray(rampArr)
            del ds

            # Calculate  the stats on the file
            ds = gdal.Open(imgfile, gdal.GA_Update)
            calcstats.calcStats(ds, progress=cuiprogress.SilentProgress(), 
                ignore=nullVal)
            del ds

            # Read back the data as a numpy array
            ds = gdal.Open(imgfile)
            band = ds.GetRasterBand(1)
            rampArr = band.ReadAsArray()

            # Get stats from file, and from array, and compare
            stats1 = getStatsFromBand(band)
            stats2 = getStatsFromArray(rampArr, nullVal)
            iterationName = "%s %s scale=%s"%(driverName, gdal.GetDataTypeName(fileDtype), scalefactor)
            # This relative tolerance is used for comparing the median and mode, 
            # because those are approximate only, and the likely error depends on the 
            # size of the numbers in question (thus it depends on the scalefactor). 
            # Please do not make it any larger unless you have a really solid reason. 
            relativeTolerance = 0.1 * scalefactor
            ok = compareStats(stats1, stats2, iterationName, relativeTolerance)
            del ds

        if os.path.exists(imgfile):
            os.remove(imgfile)
    
    if ok:
        riostestutils.report(TESTNAME, "Passed")
    else:
        riostestutils.report(TESTNAME, 
            ("Note that the mode and median tests will fail in GDAL < 2.0, "+
                "unless the GDAL fixes suggested in tickets "+
                "http://trac.osgeo.org/gdal/ticket/4750 and " +
                "http://trac.osgeo.org/gdal/ticket/5289 are applied"))
    return ok


def getStatsFromBand(band):
    """
    Get statistics from given band object, return Stats instance
    """
    mean = float(band.GetMetadataItem('STATISTICS_MEAN'))
    stddev = float(band.GetMetadataItem('STATISTICS_STDDEV'))
    minval = float(band.GetMetadataItem('STATISTICS_MINIMUM'))
    maxval = float(band.GetMetadataItem('STATISTICS_MAXIMUM'))
    median = float(band.GetMetadataItem('STATISTICS_MEDIAN'))
    mode = float(band.GetMetadataItem('STATISTICS_MODE'))
    statsObj = Stats(mean, stddev, minval, maxval, median, mode)
    return statsObj


def getStatsFromArray(arr, nullVal):
    """
    Work out the statistics directly from the image array. 
    Return a Stats instance
    """
    nonNullMask = (arr != nullVal)
    nonNullArr = arr[nonNullMask].astype(numpy.float64)
    
    # Work out what the correct answers should be
    mean = nonNullArr.mean()
    stddev = nonNullArr.std()
    minval = nonNullArr.min()
    maxval = nonNullArr.max()
    median = numpy.median(nonNullArr)
    mode = calcMode(nonNullArr, axis=None)[0][0]
    return Stats(mean, stddev, minval, maxval, median, mode)


def equalTol(a, b, tol):
    """
    Compare two values to within a tolerance. If the difference
    between the two values is smaller than the tolerance, 
    then return True
    """
    diff = abs(a - b)
    return (diff < tol)


def calcMode(a, axis=0):
    """
    Copied directly from scipy.stats.mode(), so as not to have a dependency on scipy. 
    """
    def _chk_asarray(a, axis):
        "Also copied from scipy.stats, and inserted into this function. "
        if axis is None:
            a = numpy.ravel(a)
            outaxis = 0
        else:
            a = numpy.asarray(a)
            outaxis = axis

        if a.ndim == 0:
            a = numpy.atleast_1d(a)

        return a, outaxis
        
    a, axis = _chk_asarray(a, axis)
    if a.size == 0:
        return numpy.array([]), numpy.array([])

    scores = numpy.unique(numpy.ravel(a))       # get ALL unique values
    testshape = list(a.shape)
    testshape[axis] = 1
    oldmostfreq = numpy.zeros(testshape, dtype=a.dtype)
    oldcounts = numpy.zeros(testshape, dtype=int)
    for score in scores:
        template = (a == score)
        counts = numpy.expand_dims(numpy.sum(template, axis), axis)
        mostfrequent = numpy.where(counts > oldcounts, score, oldmostfreq)
        oldcounts = numpy.maximum(counts, oldcounts)
        oldmostfreq = mostfrequent

    return (mostfrequent, oldcounts)


class Stats(object):
    def __init__(self, mean, stddev, minval, maxval, median, mode):
        self.mean = mean
        self.stddev = stddev
        self.minval = minval
        self.maxval = maxval
        self.median = median
        self.mode = mode
    
    def __str__(self):
        return ' '.join(['%s:%s'%(n, repr(getattr(self, n)))
            for n in ['mean', 'stddev', 'minval', 'maxval', 'median', 'mode']])


def compareStats(stats1, stats2, dtypeName, relativeTolerance):
    """
    Compare two Stats instances, and report differences. Also
    return True if all OK. 
    """
    ok = True
    msgList = []
    absoluteTolerance = 0.000001
    for statsName in ['mean', 'stddev', 'minval', 'maxval', 'median', 'mode']:
        value1 = getattr(stats1, statsName)
        value2 = getattr(stats2, statsName)
        
        tol = absoluteTolerance
        if statsName in ['median', 'mode']:
            tol = relativeTolerance
        if not equalTol(value1, value2, tol):
            msgList.append("Error in %s: %s (from file) != %s (from array)" % 
                (statsName, repr(value1), repr(value2)))
    
    if len(msgList) > 0:
        ok = False
        riostestutils.report(TESTNAME, 
            'Datatype=%s\n%s'%(dtypeName, '\n'.join(msgList)))
    return ok


if __name__ == "__main__":
    run()
