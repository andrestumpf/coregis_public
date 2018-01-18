
"""
This module contains routines for reading and writing Raster
Attribute Tables (RATs). These are designed to be able to 
be called from outside of RIOS.

Within RIOS, these are called from the ReaderInfo and ImageWriter
classes.

It is recommended that the ratapplier module be used instead of this
interface where possible. 
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

import sys
import os
from osgeo import gdal
import numpy
from . import rioserrors

# use turborat if available
try:
    from turbogdal import turborat
    HAVE_TURBORAT = True
except ImportError:
    HAVE_TURBORAT = False

if sys.version_info[0] > 2:
    # hack for Python 3 which uses str instead of basestring
    # we just use basestring
    basestring = str
    
# For supporting the automatic color table generation thing which Sam loves. 
DEFAULT_AUTOCOLORTABLETYPE = os.getenv('RIOS_DFLT_AUTOCOLORTABLETYPE', default=None)

def isColorColFromUsage(usage):
    "Tells if usage is one of the color column types"
    colorCol = (usage == gdal.GFU_Red or usage == gdal.GFU_Green or
                usage == gdal.GFU_Blue or usage == gdal.GFU_Alpha)
    return colorCol

def readColumnFromBand(gdalBand, colName):
    """
    Given a GDAL Band, extract the Raster Attribute with the
    given name. Returns an array of ints or floats for numeric
    data types, or a list of strings.
    
    """
    # get the RAT for this band
    rat = gdalBand.GetDefaultRAT()

    # get the size of the RAT  
    numCols = rat.GetColumnCount()
    numRows = rat.GetRowCount()
  
    # if this is still None at the end
    # we didn't find the column
    colArray = None

    # loop thru the columns looking for the right one
    for col in range(numCols):
        if rat.GetNameOfCol(col) == colName:
            # found it - create the output array
            # and fill in the values

            # use RFC40 function if available
            if hasattr(rat, "ReadAsArray"):
                colArray = rat.ReadAsArray(col)

            elif HAVE_TURBORAT:
                # if turborat is available use that
                colArray = turborat.readColumn(rat, col)
            else:
                # do it the slow way
                dtype = rat.GetTypeOfCol(col)
                if dtype == gdal.GFT_Integer:
                    colArray = numpy.zeros(numRows,int)
                elif dtype == gdal.GFT_Real:
                    colArray = numpy.zeros(numRows,float)
                elif dtype == gdal.GFT_String:
                    # for string attributes, create a list - convert later
                    colArray = []
                else:
                    msg = "Can't interpret data type of attribute"
                    raise rioserrors.AttributeTableTypeError(msg)
            
        
                # do it checking the type outside the loop for maximum speed
                if dtype == gdal.GFT_Integer:
                    for row in range(numRows):
                        val = rat.GetValueAsInt(row,col)
                        colArray[row] = val
                elif dtype == gdal.GFT_Real:
                    for row in range(numRows):
                        val = rat.GetValueAsDouble(row,col)
                        colArray[row] = val
                else:
                    for row in range(numRows):
                        val = rat.GetValueAsString(row,col)
                        colArray.append(val)

                if isinstance(colArray, list):
                    # convert to array - numpy can handle this now it can work out the lengths
                    colArray = numpy.array(colArray)

            # one last little hack - if is a colour column, but type
            # was float, multiply by 255. This is so that HFA etc that stores values
            # between 0 and 1 is consistant with its color table and other formats
            usage = rat.GetUsageOfCol(col)
            if isColorColFromUsage(usage) and rat.GetTypeOfCol(col) == gdal.GFT_Real:
                colArray = colArray * 255
                colArray = colArray.astype(int)
                
            # exit loop
            break

    # couldn't find named column - raise exception
    if colArray is None:
        msg = "Unable to find column named '%s'" % colName
        raise rioserrors.AttributeTableColumnError(msg)
    
    # return the lookup array to the caller
    return colArray

def readColumn(imgFile, colName, bandNumber=1):
    """
    Given either an open gdal dataset, or a filename,
    extract the Raster Attribute with the
    given name. Returns an array of ints or floats for numeric
    data types, or a list of strings.
    """
    if isinstance(imgFile, basestring):
        ds = gdal.Open(str(imgFile))
    elif isinstance(imgFile, gdal.Dataset):
        ds = imgFile

    gdalBand = ds.GetRasterBand(bandNumber) 

    return readColumnFromBand(gdalBand, colName)

def getColumnNamesFromBand(gdalBand):
    """
    Return the names of the columns in the attribute table
    associated with the gdalBand as a list.
    """
    # get the RAT for this band
    rat = gdalBand.GetDefaultRAT()

    colNames = []

    numCols = rat.GetColumnCount()
    for col in range(numCols):
        name = rat.GetNameOfCol(col)
        colNames.append(name)

    return colNames

def getColumnNames(imgFile, bandNumber=1):
    """
    Given either an open gdal dataset, or a filename,
    Return the names of the columns in the attribute table
    associated with the file as a list.
    """
    if isinstance(imgFile, basestring):
        ds = gdal.Open(str(imgFile))
    elif isinstance(imgFile, gdal.Dataset):
        ds = imgFile

    gdalBand = ds.GetRasterBand(bandNumber) 

    return getColumnNamesFromBand(gdalBand)


def inferColumnType(sequence):
    """
    Infer from the type of the first element in the sequence
    """
    if isinstance(sequence[0],int) or isinstance(sequence[0],numpy.integer):
        colType = gdal.GFT_Integer
    elif isinstance(sequence[0],float) or isinstance(sequence[0],numpy.floating):
        colType = gdal.GFT_Real
    elif isinstance(sequence[0],basestring) or isinstance(sequence[0],bytes):
        colType = gdal.GFT_String
    else:
        colType = None
    return colType


def writeColumnToBand(gdalBand, colName, sequence, colType=None, 
                    colUsage=gdal.GFU_Generic):
    """
    Given a GDAL band, Writes the data specified in sequence 
    (can be list, tuple or array etc)
    to the named column in the attribute table assocated with the
    gdalBand. colType must be one of gdal.GFT_Integer,gdal.GFT_Real,gdal.GFT_String.
    can specify one of the gdal.GFU_* constants for colUsage - default is 'generic'
    GDAL dataset must have been created, or opened with GA_Update
    """

    if colType is None:
        colType = inferColumnType(sequence)
    if colType is None:
        msg = "Can't infer type of column for sequence of %s" % type(sequence[0])
        raise rioserrors.AttributeTableTypeError(msg)

    # check it is acually a valid type
    elif colType not in (gdal.GFT_Integer,gdal.GFT_Real,gdal.GFT_String):
        msg = "coltype must be a valid gdal column type"
        raise rioserrors.AttributeTableTypeError(msg)

    # things get a bit weird here as we need different
    # behaviour depending on whether we have an RFC40
    # RAT or not.
    if hasattr(gdal.RasterAttributeTable, "WriteArray"):
        # new behaviour
        attrTbl = gdalBand.GetDefaultRAT()
        if attrTbl is None:
            # some formats eg ENVI return None
            # here so we need to be able to cope
            attrTbl = gdal.RasterAttributeTable()
            isFileRAT = False
        else:

            isFileRAT = True

            # but if it doesn't support dynamic writing
            # we still ahve to call SetDefaultRAT
            if not attrTbl.ChangesAreWrittenToFile():
                isFileRAT = False

    else:
        # old behaviour
        attrTbl = gdal.RasterAttributeTable()
        isFileRAT = False

    # thanks to RFC40 we need to ensure colname doesn't already exist
    colExists = False
    for n in range(attrTbl.GetColumnCount()):
        if attrTbl.GetNameOfCol(n) == colName:
            colExists = True
            colNum = n
            break
    if not colExists:
        # preserve usage
        attrTbl.CreateColumn(colName, colType, colUsage)
        colNum = attrTbl.GetColumnCount() - 1

    rowsToAdd = len(sequence)
    # Imagine has trouble if not 256 items for byte
    if gdalBand.DataType == gdal.GDT_Byte:
        rowsToAdd = 256

    # another hack to hide float (0-1) and int (0-255)
    # color table handling.
    # we assume that the column has already been created
    # of the right type appropriate for the format (maybe by calcstats)
    # Note: this only works post RFC40 when we have an actual reference
    # to the RAT rather than a new one so we can ask GetTypeOfCol
    usage = attrTbl.GetUsageOfCol(colNum)
    if (isColorColFromUsage(usage) and 
            attrTbl.GetTypeOfCol(colNum) == gdal.GFT_Real and
            colType == gdal.GFT_Integer):
        sequence = numpy.array(sequence, dtype=numpy.float)
        sequence = sequence / 255.0

    if hasattr(attrTbl, "WriteArray"):
        # if GDAL > 1.10 has these functions
        # thanks to RFC40
        attrTbl.SetRowCount(rowsToAdd)
        attrTbl.WriteArray(sequence, colNum)

    elif HAVE_TURBORAT:
        # use turborat to write values to RAT if available
        if not isinstance(sequence, numpy.ndarray):
            # turborat.writeColumn needs an array
            sequence = numpy.array(sequence)
            
        # If the dtype of the array is some unicode type, then convert to simple string type,
        # as turborat does not cope with the unicode variant. 
        if 'U' in str(sequence.dtype):
            sequence = sequence.astype(numpy.character)
            
        turborat.writeColumn(attrTbl, colNum, sequence, rowsToAdd)
    else:
        defaultValues = {gdal.GFT_Integer:0, gdal.GFT_Real:0.0, gdal.GFT_String:''}

        # go thru and set each value into the RAT
        for rowNum in range(rowsToAdd):
            if rowNum >= len(sequence):
                # they haven't given us enough values - fill in with default
                val = defaultValues[colType]
            else:
                val = sequence[rowNum]

            if colType == gdal.GFT_Integer:
                # appears that swig cannot convert numpy.int64
                # to the int type required by SetValueAsInt
                # so we need to cast. 
                # This is a problem as readColumn returns numpy.int64 
                # for integer columns. 
                # Seems fine converting numpy.float64 to 
                # float however for SetValueAsDouble.
                attrTbl.SetValueAsInt(rowNum, colNum, int(val))
            elif colType == gdal.GFT_Real:
                attrTbl.SetValueAsDouble(rowNum, colNum, float(val))
            else:
                attrTbl.SetValueAsString(rowNum, colNum, val)

    if not isFileRAT:
        # assume existing bands re-written
        # Use GDAL's exceptions to trap the error message which arises when 
        # writing to a format which does not support it
        usingExceptions = gdal.GetUseExceptions()
        gdal.UseExceptions()
        try:
            gdalBand.SetDefaultRAT(attrTbl)
        except Exception:
            pass
        if not usingExceptions:
            gdal.DontUseExceptions()


def writeColumn(imgFile, colName, sequence, colType=None, bandNumber=1, 
        colUsage=gdal.GFU_Generic):
    """
    Given either an open gdal dataset, or a filename,
    writes the data specified in sequence (can be list, tuple or array etc)
    to the named column in the attribute table assocated with the
    file. colType must be one of gdal.GFT_Integer,gdal.GFT_Real,gdal.GFT_String.
    can specify one of the gdal.GFU_* constants for colUsage - default is 'generic'
    """
    if isinstance(imgFile, basestring):
        ds = gdal.Open(str(imgFile), gdal.GA_Update)
    elif isinstance(imgFile, gdal.Dataset):
        ds = imgFile

    gdalBand = ds.GetRasterBand(bandNumber) 

    writeColumnToBand(gdalBand, colName, sequence, colType, colUsage) 

def getUsageOfColumnFromBand(gdalBand, colName):
    """
    Given a gdalBand returns the usage of the named column
    """
    # get the RAT for this band
    rat = gdalBand.GetDefaultRAT()

    # get the size of the RAT  
    numCols = rat.GetColumnCount()
  
    # loop thru the columns looking for the right one
    for col in range(numCols):
        if rat.GetNameOfCol(col) == colName:
            return rat.GetUsageOfCol(col)
    
    msg = "Unable to find column named '%s'" % colName
    raise rioserrors.AttributeTableColumnError(msg)

def getUsageOfColumn(imgFile, colName, bandNumber=1):
    """
    Given either an open gdal dataset, or a filename,
    returns the 'usage' of the column which can then be passed
    to writeColumn to preserve usage when copying
    """
    if isinstance(imgFile, basestring):
        ds = gdal.Open(str(imgFile), gdal.GA_Update)
    elif isinstance(imgFile, gdal.Dataset):
        ds = imgFile

    gdalBand = ds.GetRasterBand(bandNumber) 

    return getUsageOfColumnFromBand(gdalBand, colName)


def getColorTable(imgFile, bandNumber=1):
    """
    Given either an open gdal dataset, or a filename,
    reads the color table as an array that can be passed
    to ImageWriter.setColorTable() or rat.setColorTable()
    
    The returned colour table is a numpy array, described in detail
    in the docstring for rat.setColorTable(). 
    
    """
    if isinstance(imgFile, basestring):
        ds = gdal.Open(str(imgFile))
    elif isinstance(imgFile, gdal.Dataset):
        ds = imgFile

    gdalBand = ds.GetRasterBand(bandNumber)
    colorTable = gdalBand.GetColorTable()
    if colorTable is None:
        raise rioserrors.AttributeTableColumnError("Image has no color table")

    count = colorTable.GetCount()
    # count could be any size so we have to go with int
    colorArray = numpy.zeros((count, 5), dtype=numpy.int)
    for index in range(count):
        colorEntry = colorTable.GetColorEntry(index)
        arrayEntry = [index] + list(colorEntry)
        colorArray[index] = numpy.array(arrayEntry)

    return colorArray


def setColorTable(imgfile, colorTblArray, layernum=1):
    """
    Set the color table for the specified band. You can specify either 
    the imgfile as either a filename string or a gdal.Dataset object. The
    layer number defaults to 1, i.e. the first layer in the file. 
    
    The color table is given as a numpy array of 5 columns. There is one row 
    (i.e. first array index) for every value to be set, and the columns
    are:

        * pixelValue
        * Red
        * Green
        * Blue
        * Opacity

    The Red/Green/Blue values are on the range 0-255, with 255 meaning full 
    color, and the opacity is in the range 0-255, with 255 meaning fully 
    opaque. 
    
    The pixels values in the first column must be in ascending order, but do 
    not need to be a complete set (i.e. you don't need to supply a color for 
    every possible pixel value - any not given will default to transparent black).
    It does not even need to be contiguous. 
    
    For reasons of backwards compatability, a 4-column array will also be accepted, 
    and will be treated as though the row index corresponds to the pixelValue (i.e. 
    starting at zero). 
    
    """
    arrayShape = colorTblArray.shape
    if len(arrayShape) != 2:
        raise rioserrors.ArrayShapeError("ColorTableArray must be 2D. Found shape %s instead"%arrayShape)
        
    (numRows, numCols) = arrayShape
    # Handle the backwards-compatible case of a 4-column array
    if numCols == 4:
        pixVals = numpy.arange(numRows)
        colorTbl4cols = colorTblArray
    elif numCols == 5:
        pixVals = colorTblArray[:, 0]
        colorTbl4cols = colorTblArray[:, 1:]
    else:
        raise rioserrors.ArrayShapeError("Color table array has %d columns, expecting 4 or 5"%numCols)
    
    # Open the image file and get the band object
    if isinstance(imgfile, gdal.Dataset):
        ds = imgfile
    elif isinstance(imgfile, basestring):
        ds = gdal.Open(imgfile, gdal.GA_Update)
    
    bandobj = ds.GetRasterBand(layernum)
    
    clrTbl = gdal.ColorTable()
    maxPixVal = pixVals.max()
    i = 0
    # This loop sets an entry for every pixel value up to the largest given. Imagine
    # bitches if we don't do this. 
    tblMaxVal = maxPixVal
    if bandobj.DataType == gdal.GDT_Byte:
        # For Byte files, we always add rows for entries up to 255. Imagine gets 
        # confused if we don't
        tblMaxVal = 255
        
    for pixVal in range(tblMaxVal+1):
        while  i < numRows and pixVals[i] < pixVal:
            i += 1
        if i < numRows:
            tblPixVal = pixVals[i]
            if tblPixVal == pixVal:
                colEntry = tuple(colorTbl4cols[i])
            else:
                colEntry = (0, 0, 0, 0)
        else:
            colEntry = (0, 0, 0, 0)
        clrTbl.SetColorEntry(pixVal, colEntry)
    
    bandobj.SetRasterColorTable(clrTbl)


def genColorTable(numEntries, colortype):
    """
    Generate a colour table array. The type of colour table generated
    is controlled by the colortype string. Possible values are:
    
        * "rainbow"
        * "gray"
        * "random"
    
    See corresponding genColorTable_<colortype> function for details of each. 
    
    """
    if colortype == "rainbow":
        clrTbl = genColorTable_rainbow(numEntries)
    elif colortype == "gray":
        clrTbl = genColorTable_gray(numEntries)
    elif colortype == "random":
        clrTbl = genColorTable_random(numEntries)
    else:
        raise rioserrors.ColorTableGenerationError("Unknown colortype '{}'".format(colortype))
    return clrTbl


def genColorTable_random(numEntries):
    """
    Generate a color table array with the given number of entries by assigning 
    random red/green/blue values. No attempt is made to always generate
    unique colours, i.e. it is randomly possibly for different entries to
    have the same colour. 

    Returns the array as a 4 column array, suitable for use with the 
    :func:`rios.rat.setColorTable` function. 
    
    """
    colorTable = numpy.zeros((numEntries, 4), dtype=numpy.uint8)
    # RGB entries, random numbers. 
    colorTable[:, :3] = numpy.random.random_integers(0, 255, size=(numEntries, 3))
    # Opacity
    colorTable[:, 3] = 255
    
    return colorTable


def genColorTable_rainbow(numEntries):
    """
    Generate a color table array with the given number of entries, with colors notionally
    describing a rainbow (i.e. red-orange-yellow-green-blue-indigo-violet). Probably
    not what a painter would call a rainbow, but it will do. 
    
    Returns the array as a 4 column array, suitable for use with the 
    :func:`rios.rat.setColorTable` function. 

    """
    colorTable = numpy.zeros((numEntries, 4), dtype=numpy.uint32)
    # RGB entries. Start at red, blend through to green, then through to blue, in linear steps
    midNdx = numEntries // 2
    # Red to green
    colorTable[0:midNdx, 0] = numpy.mgrid[255:0:midNdx*1j].astype(numpy.uint8)
    colorTable[0:midNdx, 1] = numpy.mgrid[0:255:midNdx*1j].astype(numpy.uint8)
    # Green to blue
    colorTable[midNdx:, 1] = numpy.mgrid[255:0:(numEntries-midNdx)*1j].astype(numpy.uint8)
    colorTable[midNdx:, 2] = numpy.mgrid[0:255:(numEntries-midNdx)*1j].astype(numpy.uint8)
    # Opacity
    colorTable[:, 3] = 255
    
    return colorTable


def genColorTable_gray(numEntries):
    """
    Generate a color table array with the given number of entries, with all 
    colors being shades of grey. First entry is black, last entry is white. 
    
    Returns the array as a 4 column array, suitable for use with the 
    :func:`rios.rat.setColorTable` function. 

    """
    colorTable = numpy.zeros((numEntries, 4), dtype=numpy.uint8)
    # RGB entries. Start at black, blend through to white, in linear steps
    colorTable[:, 0] = numpy.mgrid[0:255:numEntries*1j].astype(numpy.uint8)
    colorTable[:, 1] = colorTable[:, 0]
    colorTable[:, 2] = colorTable[:, 0]
    # Opacity
    colorTable[:, 3] = 255
    
    return colorTable
