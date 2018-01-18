
"""
This module contains the ReaderInfo class
which holds information about the area being
read and info on the current block

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

import math

import numpy

from . import imageio
from . import rat

class StatisticsCache(object):
    """
    Allows global statistics for all the files used to be cached
    Statistics are stored here associated with the filename
    so if there are multiple calls to global_stats for one
    file the statistics will only get calculated once. 
    """
    def __init__(self):
        self.stats = {}

    def getStats(self,fname,band,ignore):
        """
        Returns tuple with (min,max,mean,stddev) if 
        previously cached, None otherwise.
        """
        if ignore is None:
            key = '%s %d' % (fname,band)
        else:
            key = '%s %d %f' % (fname,band,ignore)
        if key in self.stats:
            return self.stats[key]
        else:
            return None

    def setStats(self,fname,band,ignore,stats):
        """
        Sets tuple with (min,max,mean,stddev) in cache
        """
        if ignore is None:
            key = '%s %d' % (fname,band)
        else:
            key = '%s %d %f' % (fname,band,ignore)
        self.stats[key] = stats

class AttributeTableCache(object):
    """
    Allows attribute columns to be cached. This means
    that an attribute column can be requested for each
    block through the image, but the actual column data
    will only be read the first time and cached for 
    subsequent calls.
    """
    def __init__(self):
        self.rats = {}

    def getColumn(self, fname, band, colName):
        """
        Returns the column if previously cached, otherwise
        None.
        """
        key = '%s %d %s' % (fname, band, colName)
        if key in self.rats:
            return self.rats[key]
        else:
            return None

    def setColumn(self, fname, band, colName, column):
        """
        Puts the column into the cache
        """
        key = '%s %d %s' % (fname, band, colName)
        self.rats[key] = column


class ReaderInfo(object):
    """
    ReaderInfo class. Holds information about the area being
    read and info on the current block
    
    """
    def __init__(self, workingGrid, statscache, ratcache, 
                    windowxsize, windowysize,
                    overlap, loggingstream):
                    
        self.loggingstream = loggingstream
        # grab the working grid
        self.workingGrid = workingGrid
        
        # save the window size and overlap
        self.windowxsize = windowxsize
        self.windowysize = windowysize
        self.overlap = overlap
        
        # work out the area being read
        self.xsize = int(round((self.workingGrid.xMax - 
                        self.workingGrid.xMin) / self.workingGrid.xRes))
        self.ysize = int(round((self.workingGrid.yMax - 
                        self.workingGrid.yMin) / self.workingGrid.yRes))
        
        # total number of blocks
        self.xtotalblocks = int(math.ceil(float(self.xsize) / self.windowxsize))
        self.ytotalblocks = int(math.ceil(float(self.ysize) / self.windowysize))
        
        # save the statistics cache. 
        self.statscache = statscache

        # save the RAT cache
        self.ratcache = ratcache
        
        # The feilds below apply to a particular block
        # and are filled in after this object is copied 
        # to make it specific fir each block
        self.blockwidth = None
        self.blockheight = None
        
        self.blocktl = None
        self.blockbr = None
        
        self.xblock = None
        self.yblock = None
        
        # dictionary keyed by id() of the number array
        # value is a tuple with the GDAL dataset object 
        # that corresponds to it, and the original filename
        self.blocklookup = {}
        
    def setBlockDataset(self,block,dataset,filename):
        """
        Saves a match between the numpy block read
        and it's GDAL dataset. So we can look up the
        dataset later given a block.
        
        This routine is for internal use by RIOS. Its use in any other
        context is not sensible. 
        
        """
        self.blocklookup[id(block)] = (dataset,filename)
        
    def getWindowSize(self):
        """
        Returns the size of the current window. Returns a 
        tuple (numCols, numRows)
        
        """
        return (self.windowxsize, self.windowysize)
        
    def getOverlapSize(self):
        """
        Returns the size of the pixel overlap between
        each window. This is the number of pixels added as 
        margin around each block
        """
        return self.overlap
        
    def getTotalSize(self):
        """
        Returns the total size (in pixels) of the dataset
        being processed
        """
        return (self.xsize,self.ysize)
        
    def getTransform(self):
        """
        Return the current transform between world
        and pixel coords. This is as defined by GDAL. 
        """
        return self.workingGrid.makeGeoTransform()
        
    def getProjection(self):
        """
        Return the WKT describing the current
        projection system
        """
        return self.workingGrid.projection

    def getTotalBlocks(self):
        """
        Returns the total number of blocks the dataset
        has been split up into for processing
        """
        return (self.xtotalblocks,self.ytotalblocks)

    def setBlockSize(self,blockwidth,blockheight):
        """
        Sets the size of the current block

        This routine is for internal use by RIOS. Its use in any other
        context is not sensible. 
        
        """
        self.blockwidth = blockwidth
        self.blockheight = blockheight

    def getBlockSize(self):
        """
        Get the size of the current block. Returns a tuple::

            (numCols, numRows)

        for the current block. Mostly the same as the window size, 
        except on the edge of the raster. 
        """
        return (self.blockwidth,self.blockheight)

    def setBlockBounds(self,blocktl,blockbr):
        """
        Sets the coordinate bounds of the current block

        This routine is for internal use by RIOS. Its use in any other
        context is not sensible. 
        
        """
        self.blocktl = blocktl
        self.blockbr = blockbr

    def getBlockCoordArrays(self):
        """
        Return a tuple of the world coordinates for every pixel
        in the current block. Each array has the same shape as the 
        current block. Return value is a tuple::

            (xBlock, yBlock)

        where the values in xBlock are the X coordinates of the centre
        of each pixel, and similarly for yBlock. 
        
        The coordinates returned are for the pixel centres. This is 
        slightly inconsistent with usual GDAL usage, but more likely to
        be what one wants. 
        
        """
        (tl, br) = (self.blocktl, self.blockbr)
        (nCols, nRows) = self.getBlockSize()
        nCols += 2*self.overlap
        nRows += 2*self.overlap
        (xRes, yRes) = self.getPixelSize()
        (rowNdx, colNdx) = numpy.mgrid[0:nRows, 0:nCols]
        xBlock = tl.x - self.overlap*xRes + xRes/2.0 + colNdx * xRes
        yBlock = tl.y + self.overlap*yRes - yRes/2.0 - rowNdx * yRes
        return (xBlock, yBlock)
        
    def setBlockCount(self,xblock,yblock):
        """
        Sets the count of the current block

        This routine is for internal use by RIOS. Its use in any other
        context is not sensible. 
        
        """
        self.xblock = xblock
        self.yblock = yblock

    def getBlockCount(self):
        """
        Gets the count of the current block
        """
        return (self.xblock,self.yblock)
    
    def getPixelSize(self):
        """
        Gets the current pixel size and returns it as a tuple (x and y)
        """
        return (self.workingGrid.xRes,self.workingGrid.yRes)
    
    def getPixRowColBlock(self, x, y):
        """
        Return the row/column numbers, within the current block,
        for the pixel which contains the given (x, y) coordinate.
        The coordinates of (x, y) are in the world coordinate
        system of the reference grid. The row/col numbers are 
        suitable to use as array indices in the array(s) for the 
        current block. If the nominated pixel is not contained
        within the current block, the row and column numbers are
        both None (hence this should be checked). 
        
        Return value is a tuple of 2 int values
            (row, col)
        
        """
        transform = self.workingGrid.makeGeoTransform()
        imgRowCol = imageio.wld2pix(transform, x, y)
        imgRow = imgRowCol.y
        imgCol = imgRowCol.x
        
        blockStartRow = self.yblock * self.windowysize - self.overlap
        blockStartCol = self.xblock * self.windowxsize - self.overlap
        
        blockRow = int(imgRow - blockStartRow)
        blockCol = int(imgCol - blockStartCol)
        
        if ((blockRow < 0 or blockRow > (self.windowysize + 2 * self.overlap)) or
           (blockCol < 0 or blockCol > (self.windowxsize + 2 * self.overlap))):
            blockRow = None
            blockCol = None
        
        return (blockRow, blockCol)
    
    def getPixColRow(self,x,y):
        """
        This function is for internal use only. The user should
        be looking at getBlockCoordArrays() or getPixRowColBlock()
        for dealing with blocks and coordinates.
        
        Get the (col, row) relative to the current image grid,
        for the nominated pixel within the current block. The
        given (x, y) are column/row numbers (starting at zero),
        and the return is a tuple::

            (column, row)

        where these are relative to the whole of the current
        working grid. If working with a single raster, this is the same
        as for that raster, but if working with multiple rasters, 
        the working grid is the intersection or union of them. 
        
        Note that this function will give incorrect/misleading results
        if used in conjunction with a block overlap. 
         
        """
        col = self.xblock * self.windowxsize + x
        row = self.yblock * self.windowysize + y
        return (col,row)

    def isFirstBlock(self):
        """
        Returns True if this is the first block to be processed
        """
        return self.xblock == 0 and self.yblock == 0

    def isLastBlock(self):
        """
        Returns True if this is the last block to be processed
        """
        xtotalblocksminus1 = self.xtotalblocks - 1
        ytotalblocksminus1 = self.ytotalblocks - 1
        return self.xblock == xtotalblocksminus1 and self.yblock == ytotalblocksminus1

    def getFilenameFor(self,block):
        """
        Get the input filename of a dataset
        """
        # can't use ds.GetDescription() as may have been resampled
        (ds,fname) = self.blocklookup[id(block)]
        return fname


    def getGDALDatasetFor(self,block):
        """
        Get the underlying GDAL handle of a dataset
        """
        (ds,fname) = self.blocklookup[id(block)]
        return ds

    def getGDALBandFor(self,block,band):
        """
        Get the underlying GDAL handle for a band of a dataset
        """
        ds = self.getGDALDatasetFor(block)
        return ds.GetRasterBand(band)

    def getNoDataValueFor(self,block,band=1):
        """
        Returns the 'no data' value for the dataset
        underlying the block. This should be the
        same as what was set for the stats ignore value
        when that dataset was created. 
        The value is cast to the same data type as the 
        dataset.
        """
        ds = self.getGDALDatasetFor(block)
        band = ds.GetRasterBand(band)
        novalue = band.GetNoDataValue()

        # if there is a valid novalue, cast it to the type
        # of the dataset. Note this creates a numpy 0-d array
        if novalue is not None:
            numpytype = imageio.GDALTypeToNumpyType(band.DataType)
            novalue = numpy.cast[numpytype](novalue)

        return novalue
        
    def getPercent(self):
        """
        Returns the percent complete. 
        """
        percent = int(float(self.yblock * self.xtotalblocks + self.xblock) / 
                    float(self.xtotalblocks * self.ytotalblocks) * 100)
        return percent

    def getAttributeColumn(self, block, colName, band=1):
        """
        Gets the attribute for the given block and column name
        Caches columns so only first call actually extracts data
        """
        (ds,fname) = self.blocklookup[id(block)]

        column = self.ratcache.getColumn(fname, band, colName)
        if column is None:
            column = rat.readColumn(ds, colName, band)
            self.ratcache.setColumn(fname, band, colName, column)
        
        return column
        
    def global_stats(self,block,band=1,ignore=None):
        """
        Returns the (min,max,mean,stddev) for the whole band
        """
        fname = self.getFilenameFor(block)
        
        # see if we have the stats in our cache
        values = self.statscache.getStats(fname,band,ignore)
        
        if values is None:
            # no, get the gdal band
            bandhandle = self.getGDALBandFor(block,band)
            
            # set ignore value if specified so that 
            # GDAL ignores it when calculating stats
            if ignore is not None:
                bandhandle.SetNoDataValue(ignore)
                
            self.loggingstream.write("Calculating global statistics...\n")
                
            # get GDAL to calc the stats
            values = bandhandle.GetStatistics(False,True)
            
            # set it back in our cache for next time
            self.statscache.setStats(fname,band,ignore,values)
            
        return values

    def global_min(self,block,band=1,ignore=None):
        """
        Returns the min for the whole band
        """
        return self.global_stats(block,band,ignore)[0]

    def global_max(self,block,band=1,ignore=None):
        """
        Returns the max for the whole band
        """
        return self.global_stats(block,band,ignore)[1]

    def global_mean(self,block,band=1,ignore=None):
        """
        Returns the mean for the whole band
        """
        return self.global_stats(block,band,ignore)[2]

    def global_stddev(self,block,band=1,ignore=None):
        """
        Returns the stddev for the whole band
        """
        return self.global_stats(block,band,ignore)[3]
