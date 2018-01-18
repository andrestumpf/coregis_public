"""
Utility classes for accessing information from files, ouside of the
main RIOS applier structure. Typically these are used to access information
required to set up the call to :func:`rios.applier.apply`, passing some of the 
information in via the otherargs parameter. 

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
from osgeo import ogr
from osgeo import osr

from . import rioserrors
from . import rat

# List of datatype names corresponding to GDAL datatype numbers. 
# The index of this list corresponds to the gdal datatype number. Not sure if this 
# is a bit obscure and cryptic.....
GDALdatatypeNames = ['Unknown', 'UnsignedByte', 'UnsignedInt16', 'SignedInt16', 
    'UnsignedInt32', 'SignedInt32', 'Float32', 'Float64', 'ComplexInt16', 'ComplexInt32', 
    'ComplexFloat32', 'ComplexFloat64']

class ImageInfo(object):
    """
    An object with the bounds and other info for the given image, 
    in GDAL conventions. 
    
    Object contains the following fields
        * **xMin**            Map X coord of left edge of left-most pixel
        * **xMax**            Map X coord of right edge of right-most pixel
        * **yMin**            Map Y coord of bottom edge of bottom pixel
        * **yMax**            Map Y coord of top edge of top-most pixel
        * **xRes**            Map coord size of each pixel, in X direction
        * **yRes**            Map coord size of each pixel, in Y direction
        * **nrows**           Number of rows in image
        * **ncols**           Number of columns in image
        * **transform**       Transformation params to map between pixel and map coords, in GDAL form
        * **projection**      WKT string of projection
        * **rasterCount**     Number of rasters in file
        * **lnames**          Names of the layers as a list.
        * **layerType**       "thematic" or "athematic", if it is set
        * **dataType**        Data type for the first band (as a GDAL integer constant)
        * **dataTypeName**    Data type for the first band (as a human-readable string)
        * **nodataval**       Value used as the no-data indicator (per band)
    
    The omitPerBand argument on the constructor is provided in order to speed up the 
    access of very large VRT stacks. The information which is normally extracted 
    from each band will, in that case, trigger a gdal.Open() for each band, which 
    can be quite slow. So, if none of that information is actually required, then 
    setting omitPerBand=True will omit that information, but will return as quickly 
    as for a normal single file. 

    """
    def __init__(self, filename, omitPerBand=False):
        is_HDF_EOS_subdataset = (filename.startswith("HDF4_EOS:") or 
                filename.startswith("HDF5_EOS:"))
        if not is_HDF_EOS_subdataset and not os.path.exists(filename):
            raise rioserrors.FileOpenError("Unable to open file %s"%filename)
            
        ds = gdal.Open(str(filename), gdal.GA_ReadOnly)
        if ds is None:
            raise rioserrors.FileOpenError("Unable to open file %s"%filename)

        geotrans = ds.GetGeoTransform()
        (ncols, nrows) = (ds.RasterXSize, ds.RasterYSize)
        self.rasterCount = ds.RasterCount    
        
        self.xMin = geotrans[0]
        self.xRes = geotrans[1]
        self.yMax = geotrans[3]
        self.yRes = abs(geotrans[5])
        self.xMax = self.xMin + ncols * self.xRes
        self.yMin = self.yMax - nrows * self.yRes
        self.ncols = ncols
        self.nrows = nrows
        
        # Projection, etc. 
        self.transform = geotrans
        self.projection = ds.GetProjection()
        
        # Per-band stuff, including layer names and no data values, and stats
        self.lnames = []
        self.nodataval = []
        if not omitPerBand:
            for band in range(ds.RasterCount):
                bandObj = ds.GetRasterBand(band + 1)
                self.lnames.append(bandObj.GetDescription())
                self.nodataval.append(bandObj.GetNoDataValue())
        
        gdalMeta = ds.GetRasterBand(1).GetMetadata()
        if 'LAYER_TYPE' in gdalMeta:
            self.layerType = gdalMeta['LAYER_TYPE']
        else:
            self.layerType = None
        
        # Pixel datatype, stored as a GDAL enum value. 
        self.dataType = ds.GetRasterBand(1).DataType
        self.dataTypeName = GDALdatatypeNames[self.dataType]
        
        del ds
    
    
    def __str__(self):
        """
        Print a readable version of the object
        """
        lines = []
        for attribute in ['nrows', 'ncols', 'rasterCount', 'xMin', 'xMax', 'yMin', 'yMax', 
                'xRes', 'yRes', 'lnames', 'layerType', 'dataType', 'dataTypeName', 
                'nodataval', 'transform', 'projection']:
            value = self.__dict__[attribute]
            lines.append("%-20s%s" % (attribute, value))
        result = '\n'.join(lines)
        return result
    
    
    def layerNumberFromName(self, layerName):
        """
        Return the layer number corresponding to the given layer name.
        Valid layer numbers are as per GDAL conventions, i.e. starting at 1. 
        If the given layer name is not found in this file info, then zero is returned. 
        
        """
        try:
            ndx = self.lnames.index(layerName)
        except ValueError:
            ndx = -1
        layerNumber = ndx + 1
        return layerNumber
    
    
    def layerNameFromNumber(self, layerNumber):
        """
        Return the layer name corresponding to the given layer number. 
        Valid layer numbers are as per GDAL conventions, i.e. starting at 1.
        If the given layer number is not valid for this file info, an exception is 
        raised. 
        
        """
        ndx = layerNumber - 1
        if ndx >= 0 and ndx < self.rasterCount:
            layerName = self.lnames[ndx]
        else:
            raise rioserrors.GDALLayerNumberError("Layer number %s is outside range 1-%d"%
                    (layerNumber, self.rasterCount))
        return layerName
    
    def getCorners(self, outWKT=None, outEPSG=None, outPROJ4=None):
        """
        Return the coordinates of the image corners, possibly reprojected. 
        
        This is the same information as in the xMin, xMax, yMin, yMax fields, 
        but with the option to reproject them into a given output projection. 
        Because the output coordinate system will not in general align with the 
        image coordinate system, there are separate values for all four corners. 
        These are returned as::

            (ul_x, ul_y, ur_x, ur_y, lr_x, lr_y, ll_x, ll_y)
            
        The output projection can be given as either a WKT string, an 
        EPSG number, or a PROJ4 string. If none of those is given, then 
        bounds are not reprojected, but will be in the same coordinate 
        system as the image corners. 
        
        """
        if outWKT is not None:
            outSR = osr.SpatialReference(wkt=outWKT)
        elif outEPSG is not None:
            outSR = osr.SpatialReference()
            outSR.ImportFromEPSG(int(outEPSG))
        elif outPROJ4 is not None:
            outSR = osr.SpatialReference()
            outSR.ImportFromProj4(outPROJ4)
        else:
            outSR = None
        
        if outSR is not None:
            inSR = osr.SpatialReference(wkt=self.projection)
            t = osr.CoordinateTransformation(inSR, outSR)
            (ul_x, ul_y, z) = t.TransformPoint(self.xMin, self.yMax)
            (ll_x, ll_y, z) = t.TransformPoint(self.xMin, self.yMin)
            (ur_x, ur_y, z) = t.TransformPoint(self.xMax, self.yMax)
            (lr_x, lr_y, z) = t.TransformPoint(self.xMax, self.yMin)
        else:
            (ul_x, ul_y) = (self.xMin, self.yMax)
            (ll_x, ll_y) = (self.xMin, self.yMin)
            (ur_x, ur_y) = (self.xMax, self.yMax)
            (lr_x, lr_y) = (self.xMax, self.yMin)
        
        return (ul_x, ul_y, ur_x, ur_y, lr_x, lr_y, ll_x, ll_y)
        

class ImageLayerStats(object):
    """
    Hold the stats for a single image layer. These are as retrieved
    from the given image file, and are not calculated again. If they
    are not present in the file, they will be None. 
    Typically this class is not used separately, but only instantiated 
    as a part of the ImageFileStats class. 
    
    The object contains the following fields
        * **mean**        Mean value over all non-null pixels
        * **min**         Minimum value over all non-null pixels
        * **max**         Maximum value over all non-null pixels
        * **stddev**      Standard deviation over all non-null pixels
        * **median**      Median value over all non-null pixels
        * **mode**        Mode over all non-null pixels
        * **skipfactorx** The statistics skip factor in the X direction
        * **skipfactory** The statistics skip factor in the Y direction
        
    There are many ways to report a histogram. 
    The following attributes report it the way GDAL does. 
    See GDAL doco for precise details. 

        * **histoCounts**     Histogram counts (numpy array)
        * **histoMin**        Minimum edge of smallest bin
        * **histoMax**        Maximum edge of largest bin
        * **histoNumBins**    Number of histogram bins
        
    """
    def __init__(self, bandObj):
        metadata = bandObj.GetMetadata()
        self.mean = self.__getMetadataItem(metadata, 'STATISTICS_MEAN')
        self.stddev = self.__getMetadataItem(metadata, 'STATISTICS_STDDEV')
        self.max = self.__getMetadataItem(metadata, 'STATISTICS_MAXIMUM')
        self.min = self.__getMetadataItem(metadata, 'STATISTICS_MINIMUM')
        self.median = self.__getMetadataItem(metadata, 'STATISTICS_MEDIAN')
        self.mode = self.__getMetadataItem(metadata, 'STATISTICS_MODE')
        self.skipfactorx = self.__getMetadataItem(metadata, 'STATISTICS_SKIPFACTORX')
        self.skipfactory = self.__getMetadataItem(metadata, 'STATISTICS_SKIPFACTORY')
        
        self.histoMin = self.__getMetadataItem(metadata, 'STATISTICS_HISTOMIN')
        self.histoMax = self.__getMetadataItem(metadata, 'STATISTICS_HISTOMAX')
        self.histoNumBins = self.__getMetadataItem(metadata, 'STATISTICS_HISTONUMBINS')

        self.histoCounts = None
        # check the histo info - from RAT if available
        # and we are using RFC40
        rat = bandObj.GetDefaultRAT()
        if rat is not None and hasattr(rat, "ReadAsArray"):
            for col in range(rat.GetColumnCount()):
                if rat.GetUsageOfCol(col) == gdal.GFU_PixelCount:
                    self.histoCounts = rat.ReadAsArray(col)
                    break

        if self.histoCounts is None and 'STATISTICS_HISTOBINVALUES' in metadata:
            histoString = metadata['STATISTICS_HISTOBINVALUES']
            histoStringList = [c for c in histoString.split('|') if len(c) > 0]
            counts = [eval(c) for c in histoStringList]
            self.histoCounts = numpy.array(counts)
    
    @staticmethod
    def __getMetadataItem(metadata, key):
        "Return eval(item) by key, or None if not present"
        item = None
        if key in metadata:
            item = eval(metadata[key])
        return item
    
    def __str__(self):
        "Readable string representation of stats"
        fmt = "Mean: %s, Stddev: %s, Min: %s, Max: %s, Median: %s, Mode: %s"
        return (fmt % (self.mean, self.stddev, self.min, self.max, self.median, self.mode))


class ImageFileStats(object):
    """
    Hold the stats for all layers in an image file. This object can be indexed 
    with the layer index, and each element is an instance of :class:`rios.fileinfo.ImageLayerStats`. 
    
    """
    def __init__(self, filename):
        ds = gdal.Open(filename)
        self.statsList = []
        for i in range(ds.RasterCount):
            bandObj = ds.GetRasterBand(i+1)
            self.statsList.append(ImageLayerStats(bandObj))
        del ds
    
    def __getitem__(self, i):
        return self.statsList[i]
    
    def __len__(self):
        return len(self.statsList)

    def __str__(self):
        return '\n'.join([str(s) for s in self.statsList])


class VectorFileInfo(object):
    """
    Hold useful general information about a vector file. This object
    can be indexed with the layer index, and each element is
    an instance of :class:`rios.fileinfo.VectorLayerInfo`. 
    
    """
    def __init__(self, filename):
        ds = ogr.Open(filename)
        if ds is None:
            raise rioserrors.VectorLayerError("Unable to open vector dataset '%s'"%filename)
        layerCount = ds.GetLayerCount()
        self.layerInfo = [VectorLayerInfo(ds, i) for i in range(layerCount)]
    
    def __getitem__(self, i):
        return self.layerInfo[i]
    
    def __str__(self):
        return '\n'.join(['Layer:%s\n%s'%(i, str(self.layerInfo[i])) 
                for i in range(len(self.layerInfo))])


geometryTypeStringDict = {
    1:'Point',
    2:'Line',
    3:'Polygon'
}
class VectorLayerInfo(object):
    """
    Hold useful general information about a single vector layer. 
    
    Object contains the following fields
        * **featureCount**        Number of features in the layer
        * **xMin**                Minimum X coordinate
        * **xMax**                Maximum X coordinate
        * **yMin**                Minimum Y coordinate
        * **yMax**                Maximum Y coordinate
        * **geomType**            OGR geometry type code (integer)
        * **geomTypeStr**         Human-readable geometry type name (string)
        * **fieldCount**          Number of fields (i.e. columns) in attribute table
        * **fieldNames**          List of names of attribute table fields
        * **fieldTypes**          List of the type code numbers of each attribute table field
        * **fieldTypeNames**      List of the string names of the field types
        * **spatialRef**          osr.SpatialReference object of layer projection
    
    """
    def __init__(self, ds, i):
        lyr = ds.GetLayer(i)
        if lyr is None:
            raise rioserrors.VectorLayerError("Unable to open layer %s in dataset '%s'"%(i, ds.GetName()))
        
        self.featureCount = lyr.GetFeatureCount()
        extent = lyr.GetExtent()
        self.xMin = extent[0]
        self.xMax = extent[1]
        self.yMin = extent[2]
        self.yMax = extent[3]
        
        self.geomType = lyr.GetGeomType()
        if self.geomType in geometryTypeStringDict:
            self.geomTypeStr = geometryTypeStringDict[self.geomType]
        
        lyrDefn = lyr.GetLayerDefn()
        self.fieldCount = lyrDefn.GetFieldCount()
        fieldDefnList = [lyrDefn.GetFieldDefn(i) for i in range(self.fieldCount)]
        self.fieldNames = [fd.GetName() for fd in fieldDefnList]
        self.fieldTypes = [fd.GetType() for fd in fieldDefnList]
        self.fieldTypeNames = [fd.GetTypeName() for fd in fieldDefnList]
        
        self.spatialRef = lyr.GetSpatialRef()

    def __str__(self):
        valueList = []
        for valName in ['featureCount', 'xMin', 'xMax', 'yMin', 'yMax', 'geomType', 
                'geomTypeStr', 'fieldCount', 'fieldNames', 'fieldTypes', 'fieldTypeNames',
                'spatialRef']:
            valueList.append("  %s: %s" % (valName, getattr(self, valName)))
        return '\n'.join(valueList)


class ColumnStats(object):
    """
    Summary statistics for a single column of a raster attribute table
    
    Attributes are:
        * **count**
        * **mean**
        * **stddev**
        * **min**
        * **max**
        * **mode**            Not yet implemented
        * **median**          Not yet implemented
    
    Unless otherwise requested, the statistics are weighted by the Histogram column, so
    that they represent spatial statistics, e.g. the mean would correspond to the mean
    over all pixels, rather than the mean over all rows, and the count is effectively
    the number of pixels, rather than the number of rows. This can be disabled by
    setting histogramweighted=False when calling the constructor. 
    
    Unless otherwise requested, rows which correspond to the image null value are not
    included in stats calculation. This can be reversed by setting includeImageNull=True
    on the constructor, although normally the histogram does not include counts for the 
    null value either, so this only really makes sense when using histogramweighted=False. 
    (N.B. that last point assumes the presence of the GDAL patches for 
    tickets #4750 and #5289. )
        
    """
    def __init__(self, band, columnName, includeImageNull=False, histogramweighted=True):
        gdalRat = band.GetDefaultRAT()
        columnNdx = self.__findColumnNdx(gdalRat, columnName)
        if columnNdx == -1:
            raise rioserrors.AttributeTableColumnError("Cannot find RAT column '%s'"%columnName)
        
        if gdalRat.GetTypeOfCol(columnNdx) == gdal.GFT_String:
            # The values in the column are strings, so do nothing. This essentially
            # creates an empty object, hence the RatStats class checks this and
            # does not make use of it. 
            return

        histoColumnNdx = gdalRat.GetColOfUsage(gdal.GFU_PixelCount)
        # Kludge, if GDAL is missing HFA patch for ticket #5359
        if histoColumnNdx == -1:
            histoColumnNdx = self.__findColumnNdx(gdalRat, "Histogram")

        # Check if we have the features available in GDAL RFC40  
        # See http://trac.osgeo.org/gdal/wiki/rfc40_enhanced_rat_support
        haveRFC40 = hasattr(gdalRat, 'ReadAsArray')
        numRows = gdalRat.GetRowCount()
        if not haveRFC40:
            blocklen = numRows
        else:
            blocklen = 10000
        numBlocks = int(numpy.ceil(float(numRows) / blocklen))
        
        # Initialise sums and counters
        sumX = ssqX = count = 0
        self.min = self.max = None
        
        # Loop over all blocks of data
        for i in range(numBlocks):
            startrow = i * blocklen
            endrow = min(startrow + blocklen - 1, numRows-1)
            if haveRFC40:
                datablock = gdalRat.ReadAsArray(columnNdx, start=startrow, length=(endrow-startrow+1))
                histoblock = gdalRat.ReadAsArray(histoColumnNdx, start=startrow, length=(endrow-startrow+1))
            else:
                # Without RFC40, we have no choice but to read the whole column
                datablock = rat.readColumnFromBand(band, columnName)
                histoColumnName = gdalRat.GetNameOfCol(histoColumnNdx)
                histoblock = rat.readColumnFromBand(band, histoColumnName)
        
            imgNullVal = band.GetNoDataValue()
            if not includeImageNull and imgNullVal is not None:
                pixelvals = numpy.arange(startrow, endrow+1, dtype=numpy.uint32)
                nonNullMask = (pixelvals != imgNullVal)
                datablock = datablock[nonNullMask]
                histoblock = histoblock[nonNullMask]
                del nonNullMask, pixelvals

            # Do sum calculations on float64 array, to avoid integer wrap-arounds and overflows. 
            # Note that this still implies a loss of precision, but this is a known feature
            # of doing incremental mean/stddev calculations. 
            dataAsFloat = datablock.astype(numpy.float64)

            # If requested, weight the values by their histogram counts
            weights = numpy.ones(len(datablock), dtype=numpy.uint8)
            if histogramweighted:
                weights = histoblock

            sumX += (dataAsFloat * weights).sum()
            ssqX += ((dataAsFloat**2) * weights).sum()
            count += weights.sum()

            # Min and max values
            blockMin = datablock[weights>0].min()
            blockMax = datablock[weights>0].max()
            if self.min is None or blockMin < self.min:
                self.min = blockMin
            if self.max is None or blockMax > self.max:
                self.max = blockMax
        
        self.count = count
        
        # Finish mean and stddev
        self.mean = None
        if count > 0:
            self.mean = sumX / count
        self.stddev = None
        if count > 0:
            # Use the naive formula, even though we know that it is prone to loss of precision. 
            # Should we instead be using a 2-pass method? 
            self.stddev = numpy.sqrt(ssqX / count - self.mean**2)
            
        self.median = None
        self.mode = None
    
    @staticmethod
    def __findColumnNdx(gdalRat, columnName):
        """
        Utility routine to find the column index for the given column name, in 
        the given gdalRat
        """
        columnNdx = -1
        for i in range(gdalRat.GetColumnCount()):
            if gdalRat.GetNameOfCol(i) == columnName:
                columnNdx = i
        return columnNdx
    
    
    def __str__(self):
        "Readable string representation of stats"
        fmt = "Count: %s, Mean: %s, Stddev: %s, Min: %s, Max: %s, Median: %s, Mode: %s"
        return (fmt % (self.count, self.mean, self.stddev, self.min, self.max, self.median, self.mode))
            

class RatStats(object):
    """
    Calculate statistics on columns in a Raster Attribute Table

    Normal usage is via the RatStats class, e.g.::

        columnsOfInterest = ['col1', 'col4']
        ratStatsObj = ratstats.RatStats('file.img', columnlist=columnsOfInterest)
    
        print ratStatsObj.col1.mean, ratStatsObj.col4.mean

    Contains an attribute for each named column. Each attribute is 
    a ColumnStats object, containing all relevant global stats 
    for that column. See docstring for that class for details. 

    """
    def __init__(self, filename, columnlist=None, layernum=1, includeImageNull=False, histogramweighted=True):
        """
        Default columnlist is all columns in the table. 
    
        """
        ds = gdal.Open(filename)
        band = ds.GetRasterBand(layernum)
        
        if columnlist is None:
            columnlist = rat.getColumnNamesFromBand(band)
        
        for columnName in columnlist:
            stats = ColumnStats(band, columnName, includeImageNull, histogramweighted)
            if hasattr(stats, 'mean'):
                setattr(self, columnName, stats)
