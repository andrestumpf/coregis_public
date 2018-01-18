"""
This module creates pyramid layers and calculates statistics for image
files. Much of it was originally for ERDAS Imagine files but should work 
with any other format that supports pyramid layers and statistics

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
from . import cuiprogress
from .rioserrors import ProcessCancelledError

# Test whether we have access to the GDAL RFC40 facilities
haveRFC40 = False
if (os.getenv('RIOS_HISTOGRAM_IGNORE_RFC40') is None and 
        hasattr(gdal.RasterAttributeTable, 'ReadAsArray')):
    haveRFC40 = True


# When calculating overviews (i.e. pyramid layers), default behaviour
# is controlled by these
dfltOverviewLvls = os.getenv('RIOS_DFLT_OVERVIEWLEVELS')
if dfltOverviewLvls is None:
    DEFAULT_OVERVIEWLEVELS = [ 4, 8, 16, 32, 64, 128, 256, 512 ]
else:
    DEFAULT_OVERVIEWLEVELS = [int(i) for i in dfltOverviewLvls.split(',')]
DEFAULT_MINOVERVIEWDIM = int(os.getenv('RIOS_DFLT_MINOVERLEVELDIM', default=33))
DEFAULT_OVERVIEWAGGREGRATIONTYPE = os.getenv('RIOS_DFLT_OVERVIEWAGGTYPE', 
    default="AVERAGE")


def progressFunc(value,string,userdata):
    """
    Progress callback for BuildOverviews
    """
    percent = (userdata.curroffset + (value / userdata.nbands) * 100)
    userdata.progress.setProgress(percent)
    if value == 1.0:
        userdata.curroffset = userdata.curroffset + 100.0 / userdata.nbands
    return not userdata.progress.wasCancelled()
  
# make userdata object with progress and num bands
class ProgressUserData:
    pass

def addPyramid(ds, progress, 
        minoverviewdim=DEFAULT_MINOVERVIEWDIM, 
        levels=DEFAULT_OVERVIEWLEVELS,
        aggregationType=None):
    """
    Adds Pyramid layers to the dataset. Adds levels until
    the raster dimension of the overview layer is < minoverviewdim,
    up to a maximum level controlled by the levels parameter. 
    
    Uses gdal.Dataset.BuildOverviews() to do the work. 
    
    """
    progress.setLabelText("Computing Pyramid Layers...")
    progress.setProgress(0)

    # first we work out how many overviews to build based on the size
    if ds.RasterXSize < ds.RasterYSize:
        mindim = ds.RasterXSize
    else:
        mindim = ds.RasterYSize
    
    nOverviews = 0
    for i in levels:
        if (mindim // i ) > minoverviewdim:
            nOverviews = nOverviews + 1

    # Need to find out if we are thematic or continuous. 
    tmpmeta = ds.GetRasterBand(1).GetMetadata()
    if aggregationType is None:
        if 'LAYER_TYPE' in tmpmeta:
            if tmpmeta['LAYER_TYPE'] == 'athematic':
                aggregationType = "AVERAGE"
            else:
                aggregationType = "NEAREST"
        else:
            aggregationType = DEFAULT_OVERVIEWAGGREGRATIONTYPE
    
    userdata = ProgressUserData()
    userdata.progress = progress
    userdata.nbands = ds.RasterCount
    userdata.curroffset = 0
   
    ds.BuildOverviews(aggregationType, levels[:nOverviews], progressFunc, userdata )
  
    if progress.wasCancelled():
        raise ProcessCancelledError()

    # make sure it goes to 100%
    progress.setProgress(100)

def findOrCreateColumn(ratObj, usage, name, dtype):
    """
    Returns the index of an existing column matched
    on usage. Creates it if not already existing using 
    the supplied name and dtype
    Returns a tupe with index and a boolean specifying if 
    it is a new column or not
    """
    ncols = ratObj.GetColumnCount()
    for col in range(ncols):
        if ratObj.GetUsageOfCol(col) == usage:
            return col, False

    # got here so can't exist
    ratObj.CreateColumn(name, dtype, usage)
    # new one will be last col
    return ncols, True

gdalLargeIntTypes = set([gdal.GDT_Int16, gdal.GDT_UInt16, gdal.GDT_Int32, gdal.GDT_UInt32])
gdalFloatTypes = set([gdal.GDT_Float32, gdal.GDT_Float64])
def addStatistics(ds,progress,ignore=None):
    """
    Calculates statistics and adds them to the image
    
    Uses gdal.Band.ComputeStatistics() for mean, stddev, min and max,
    and gdal.Band.GetHistogram() to do histogram calculation. 
    The median and mode are estimated using the histogram, and so 
    for larger datatypes, they will be approximate only. 
    
    For thematic layers, the histogram is calculated with as many bins 
    as required, for athematic integer and float types, a maximum
    of 256 bins is used. 
    
    """
    progress.setLabelText("Computing Statistics...")
    progress.setProgress(0)
    percent = 0
    percentstep = 100.0 / (ds.RasterCount * 2) # 2 steps for each layer
  
    for bandnum in range(ds.RasterCount):
        band = ds.GetRasterBand(bandnum + 1)

        # fill in the metadata
        tmpmeta = band.GetMetadata()
    
        if ignore is not None:
            # tell QGIS that the ignore value was ignored
            band.SetNoDataValue(ignore)
            tmpmeta["STATISTICS_EXCLUDEDVALUES"] = repr(ignore) # doesn't seem to do anything
      
        # get GDAL to calculate statistics - force recalculation. Trap errors 
        useExceptions = gdal.GetUseExceptions()
        gdal.UseExceptions()
        try:
            (minval,maxval,meanval,stddevval) = band.ComputeStatistics(False)
        except RuntimeError as e:
            if str(e).endswith('Failed to compute statistics, no valid pixels found in sampling.'):
                minval = ignore
                maxval = ignore
                meanval = ignore
                stddevval = 0
        if not useExceptions:
            gdal.DontUseExceptions()

        percent = percent + percentstep
        progress.setProgress(percent)
    
        tmpmeta["STATISTICS_MINIMUM"] = repr(minval)
        tmpmeta["STATISTICS_MAXIMUM"] = repr(maxval)
        tmpmeta["STATISTICS_MEAN"]    = repr(meanval)
        tmpmeta["STATISTICS_STDDEV"]  = repr(stddevval)
        # because we did at full res - these are the default anyway
        tmpmeta["STATISTICS_SKIPFACTORX"] = "1"
        tmpmeta["STATISTICS_SKIPFACTORY"] = "1"

        # create a histogram so we can do the mode and median
        if band.DataType == gdal.GDT_Byte:
            # if byte data use 256 bins and the whole range
            histmin = 0
            histmax = 255
            histstep = 1.0
            histCalcMin = -0.5
            histCalcMax = 255.5
            histnbins = 256
            tmpmeta["STATISTICS_HISTOBINFUNCTION"] = 'direct'
        elif "LAYER_TYPE" in tmpmeta and tmpmeta["LAYER_TYPE"] == 'thematic':
            # all other thematic types a bin per value
            histmin = 0
            histmax = int(numpy.ceil(maxval))
            histstep = 1.0
            histCalcMin = -0.5
            histCalcMax = maxval + 0.5
            histnbins = histmax + 1
            tmpmeta["STATISTICS_HISTOBINFUNCTION"] = 'direct'
        elif band.DataType in gdalLargeIntTypes:
            histrange = int(numpy.ceil(maxval) - numpy.floor(minval)) + 1
            (histmin, histmax) = (minval, maxval)
            if histrange <= 256:
                histnbins = histrange
                histstep = 1.0
                tmpmeta["STATISTICS_HISTOBINFUNCTION"] = 'direct'
                histCalcMin = histmin - 0.5
                histCalcMax = histmax + 0.5
            else:
                histnbins = 256
                tmpmeta["STATISTICS_HISTOBINFUNCTION"] = 'linear'
                histCalcMin = histmin
                histCalcMax = histmax
                histstep = float(histCalcMax - histCalcMin) / histnbins
        elif band.DataType in gdalFloatTypes:
            histnbins = 256
            (histmin, histmax) = (minval, maxval)
            tmpmeta["STATISTICS_HISTOBINFUNCTION"] = 'linear'
            histCalcMin = minval
            histCalcMax = maxval
            histstep = float(histCalcMax - histCalcMin) / histnbins
        # Note that the complex number data types are not handled, as I am not sure
        # what a histogram or a median would mean for such types. 
      
        userdata = ProgressUserData()
        userdata.progress = progress
        userdata.nbands = ds.RasterCount * 2
        userdata.curroffset = percent
      
        # get histogram and force GDAL to recalculate it
        hist = band.GetHistogram(histCalcMin, histCalcMax, histnbins, False, 
                        False, progressFunc, userdata)
        
        # Check if GDAL's histogram code overflowed. This is not a fool-proof test,
        # as some overflows will not result in negative counts. 
        histogramOverflow = (min(hist) < 0)
        
        # we may use this ratObj reference for the colours below also
        # may be None if format does not support RATs
        ratObj = band.GetDefaultRAT()

        if not histogramOverflow:
            # comes back as a list for some reason
            hist = numpy.array(hist)

            # Note that we have explicitly set histstep in each datatype case 
            # above. In principle, this can be calculated, as it is done in the 
            # float case, but for some of the others we need it to be exactly
            # equal to 1, so we set it explicitly there, to avoid rounding
            # error problems. 

            # do the mode - bin with the highest count
            modebin = numpy.argmax(hist)
            modeval = modebin * histstep + histmin
            if band.DataType == gdal.GDT_Float32 or band.DataType == gdal.GDT_Float64:
                tmpmeta["STATISTICS_MODE"] = repr(modeval)
            else:
                tmpmeta["STATISTICS_MODE"] = repr(int(round(modeval)))

            tmpmeta["STATISTICS_HISTOMIN"] = repr(histmin)
            tmpmeta["STATISTICS_HISTOMAX"] = repr(histmax)
            tmpmeta["STATISTICS_HISTONUMBINS"] = repr(histnbins)

            if haveRFC40 and ratObj is not None:
                histIndx, histNew = findOrCreateColumn(ratObj, gdal.GFU_PixelCount, 
                                        "Histogram", gdal.GFT_Real)
                # write the hist in a single go
                ratObj.SetRowCount(histnbins)
                ratObj.WriteArray(hist, histIndx)

                # The HFA driver still honours the STATISTICS_HISTOBINVALUES
                # metadata item. If we are recalculating the histogram the old
                # values will be copied across with the metadata so clobber it
                if "STATISTICS_HISTOBINVALUES" in tmpmeta:
                    del tmpmeta["STATISTICS_HISTOBINVALUES"]
            else:
                # old method
                tmpmeta["STATISTICS_HISTOBINVALUES"] = '|'.join(map(repr,hist)) + '|'

            # estimate the median - bin with the middle number
            middlenum = hist.sum() / 2
            gtmiddle = hist.cumsum() >= middlenum
            medianbin = gtmiddle.nonzero()[0][0]
            medianval = medianbin * histstep + histmin
            if band.DataType == gdal.GDT_Float32 or band.DataType == gdal.GDT_Float64:
                tmpmeta["STATISTICS_MEDIAN"]  = repr(medianval)
            else:
                tmpmeta["STATISTICS_MEDIAN"]  = repr(int(round(medianval)))
    
        # set the data
        band.SetMetadata(tmpmeta)

        if haveRFC40 and ratObj is not None and not ratObj.ChangesAreWrittenToFile():
            # For drivers that require the in memory thing
            band.SetDefaultRAT(ratObj)

        percent = percent + percentstep
        progress.setProgress(percent)

        if progress.wasCancelled():
            raise ProcessCancelledError()
    
    progress.setProgress(100)
    
    
def calcStats(ds,progress=None,ignore=None,
        minoverviewdim=DEFAULT_MINOVERVIEWDIM, 
        levels=DEFAULT_OVERVIEWLEVELS,
        aggregationType=None):
    """
    Does both the stats and pyramid layers. Calls addStatistics()
    and addPyramid() functions. See their docstrings for details. 
    
    """
    if progress is None:
        progress = cuiprogress.CUIProgressBar()
        
    addStatistics(ds,progress,ignore)
    
    addPyramid(ds, progress, minoverviewdim=minoverviewdim, levels=levels, 
        aggregationType=aggregationType)

