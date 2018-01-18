
"""
This module contains the Vector and VectorReader
class that perform on-the-fly rasterization of
vectors into raster blocks to fit in with
ImageReader and ImageWriter classes.

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
import sys
import tempfile
import subprocess
from .imagewriter import DEFAULTDRIVERNAME
from .imagewriter import DEFAULTCREATIONOPTIONS
from . import rioserrors
from . import cuiprogress
from .imagereader import ImageReader
from .imageio import NumpyTypeToGDALType
from osgeo import ogr
from osgeo import gdal
from osgeo import osr
import numpy

DEFAULTBURNVALUE = 1

class Vector(object):
    """
    Class that holds information about a vector dataset and how it
    should be rasterized. Used for passing to VectorReader.
    """
    def __init__(self, filename, inputlayer=0, burnvalue=DEFAULTBURNVALUE,
                    attribute=None, filter=None, alltouched=False, datatype=numpy.uint8, 
                    tempdir='.', driver=DEFAULTDRIVERNAME, 
                    driveroptions=DEFAULTCREATIONOPTIONS,
                    nullval=0):
        """
        Constructs a Vector object. filename should be a path to an OGR readable
        dataset. 
        inputlayer should be the OGR layer number (or name) in the dataset to rasterize.
        burnvalue is the value that gets written into the raster inside a polygon.
        Alternatively, attribute may be the name of an attribute that is looked up
        for the value to write into the raster inside a polygon.
        If you want to filter the attributes in the vector, pass filter which
        is in the format of an SQL WHERE clause.
        By default, only pixels whose centre lies within the a polygon get
        rasterised. To have all pixels touched by a polygon rasterised, set
        alltouched=True. 
        datatype is the numpy type to rasterize to - byte by default.
        tempdir is the directory to create the temporary raster file.
        driver and driveroptions set the GDAL raster driver name and options
        for the temporary rasterised file.

        """
        # open the file and get the requested layer
        self.filename = filename
        self.layerid = inputlayer
        self.ds = ogr.Open(filename)
        if self.ds is None:
            raise rioserrors.ImageOpenError("Unable to open OGR dataset: %s" % filename)
        self.layer = self.ds.GetLayer(inputlayer)
        if self.layer is None:
            raise rioserrors.VectorLayerError("Unable to find layer: %s" % inputlayer)
        layerdefn = self.layer.GetLayerDefn()

        # check the attribute exists
        if attribute is not None:
            fieldidx = layerdefn.GetFieldIndex(attribute)
            if fieldidx == -1:
                raise rioserrors.VectorAttributeError("Attribute does not exist in file: %s" % attribute)

        # check they have passed a polygon type
        validtypes = [ogr.wkbMultiPolygon,ogr.wkbMultiPolygon25D,ogr.wkbPolygon,ogr.wkbPolygon25D]
        if layerdefn.GetGeomType() not in validtypes:
            gdalVersion = None
            if hasattr(gdal, '__version__'):
                gdalVersion = gdal.__version__
            # This seems to be the only way to reliably deal with 
            # GDAL 1.10.0 < 1.9.0 comparisons...
            from distutils.version import LooseVersion
            if gdalVersion is None or LooseVersion(gdalVersion) < LooseVersion('1.9.0'):
                raise rioserrors.VectorGeometryTypeError("Can only rasterize polygon types "+
                    "with this version of gdal. Need gdal version >= 1.9.0")

        # apply the attribute filter if passed
        if filter is not None:
            self.layer.SetAttributeFilter(filter)     
        # store it in case of reprojection
        self.filter = filter

        # create a temporary file name based on raster
        # driver extension
        # save GDAL driver object for dataset creation later
        self.driver = gdal.GetDriverByName(driver)
        drivermeta = self.driver.GetMetadata()
        ext = ''
        if gdal.DMD_EXTENSION in drivermeta:
            ext = '.' + drivermeta[gdal.DMD_EXTENSION]
        # save the driver options
        self.driveroptions = driveroptions

        self.tempdir = tempdir
        (fileh,self.temp_image) = tempfile.mkstemp(ext,dir=tempdir)
        # close the file so we can get GDAL to clobber it
        # probably a security hole - not sure
        os.close(fileh)

        # create the options string
        self.options = []
        if attribute is not None:
            self.options.append('ATTRIBUTE=%s' % attribute)
        if alltouched:
            self.options.append('ALL_TOUCHED=TRUE')

        # store the data type
        self.datatype = datatype
        # burnvalue
        self.burnvalue = burnvalue
        # Value used for area not burned
        self.nullval = nullval

    def cleanup(self):
        """
        Remove temporary file(s) and close dataset

        """
        if os.path.exists(self.temp_image):
            self.rasterDS = None
            drvr = gdal.IdentifyDriver(self.temp_image)
            drvr.Delete(self.temp_image)

        del self.layer
        self.layer = None
        del self.ds
        self.ds = None
        if hasattr(self, 'reprojectedFile'):
            drvr = self.reprojectedDS.GetDriver()
            delattr(self, 'reprojectedDS')
            drvr.DeleteDataSource(self.reprojectedFile)
            delattr(self, 'reprojectedFile')

    def __del__(self):
        # destructor - call cleanup
        self.cleanup()
    
    def matchingProj(self, proj):
        """
        Returns True if the current vector has the same projection
        as the given projection. The projectino can be given as either
        a WKT string, or an osr.SpatialReference instance. 
        
        """
        if isinstance(proj, osr.SpatialReference):
            sr = proj
        else:
            sr = osr.SpatialReference(wkt=proj)
        selfSr = self.layer.GetSpatialRef()
        return selfSr.IsSame(sr)
    
    def reproject(self, proj):
        """
        Reproject the current vector to the given projection. Places the
        result in a temporary shapefile, and opens it. Returns the ogr.Layer
        object resulting. 
         
        Assumes that shapefile format will always be sufficient - I can't think 
        of any reason why not. 
        
        This reprojection becomes a part of the "state" of the current object, i.e.
        there can only be one "reprojected" copy of the current vector. This should
        be fine. 
        
        """
        if isinstance(proj, osr.SpatialReference):
            projWKT = proj.ExportAsWkt()
        else:
            projWKT = proj
        
        (fd, tmpVectorfile) = tempfile.mkstemp(prefix='tmp', suffix='.shp',dir=self.tempdir)
        os.close(fd)
        # This is naughty, but otherwise ogr2ogr won't work
        os.remove(tmpVectorfile)
        cmdList = ["ogr2ogr", '-f', "ESRI Shapefile", '-t_srs', projWKT,
            tmpVectorfile, self.filename]
        proc = subprocess.Popen(cmdList, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
        (stdoutStr, stderrStr) = proc.communicate()
        # sometimes warnings etc printed to stderr so we can't
        # rely on that for testing success. Use returncode instead.
        if proc.returncode != 0:
            msg = "Trouble reprojecting vector\n\n"+stdoutStr+'\n'+stderrStr
            raise rioserrors.VectorProjectionError(msg)
        
        self.reprojectedFile = tmpVectorfile
        self.reprojectedDS = ogr.Open(tmpVectorfile)
        layer = self.reprojectedDS.GetLayer(self.layerid)

        # apply the attribute filter if in original
        if self.filter is not None:
            layer.SetAttributeFilter(self.filter)

        return layer


def rasterizeProgressFunc(value, string, progressObj):
    """
    called by gdal.RasterizeLayer
    """
    percent = int(value * 100)
    progressObj.setProgress(percent)

class VectorReader(object):
    """
    Class that performs rasterization of Vector objects.

    """
    def __init__(self, vectorContainer, progress=None):
        """
        vectorContainer is a single Vector object, or a 
        list or dictionary that contains
        the Vector objects of the files to be read.
        If a Vector object is passed, a single block is returned
        from rasterize(), if a list is passed, 
        a list of blocks is returned, if a dictionary a dictionary is
        returned for each call to rasterize() with the same keys.
        progress is an instance of a Progress class, if none 
        an instance of cuiprogress.CUIProgress is created an used
        """
        self.vectorContainer = vectorContainer
        if progress is None:
            self.progress = cuiprogress.CUIProgressBar()
        else:
            self.progress = progress

    @staticmethod
    def rasterizeSingle(info, vector, progress):
        """
        Static method to rasterize a single Vector for the extents
        specified in the info object. 
        
        For efficiency, it rasterizes the whole working grid, and caches 
        this, and then reads the relevant section of the gridd for the 
        current block. 
        
        Will reproject the vector into the working projection, if required. 
        
        A single numpy array is returned of rasterized data.
        
        """
        try:
            if info.isFirstBlock():
                projection = info.getProjection()
                veclayer = vector.layer
                if not vector.matchingProj(projection):
                    veclayer = vector.reproject(projection)
                    
                # Haven't yet rasterized, so do this for the whole workingGrid
                (nrows, ncols) = info.workingGrid.getDimensions()
                numLayers = 1
                gdaldatatype = NumpyTypeToGDALType(vector.datatype)
                outds = vector.driver.Create(vector.temp_image, ncols, nrows, numLayers, 
                    gdaldatatype, vector.driveroptions)
                if outds is None:
                    raise rioserrors.ImageOpenError("Unable to create temporary file %s" % vector.temp_image)
                outds.SetGeoTransform(info.getTransform())
                outds.SetProjection(projection)
                # Fill raster with vector null value
                for i in range(numLayers):
                    band = outds.GetRasterBand(i+1)
                    band.Fill(vector.nullval)

                progress.setLabelText("Rasterizing...")
                err = gdal.RasterizeLayer(outds, [1], veclayer, burn_values=[vector.burnvalue], 
                                        options=vector.options, callback=rasterizeProgressFunc,
                                        callback_data=progress)
                progress.reset()

                if err != gdal.CE_None:
                    raise rioserrors.VectorRasterizationError("Rasterization failed")
                
                vector.rasterDS = outds
        except Exception:
            # if there has been an exception
            # ensure all the files are cleaned up
            vector.cleanup()
            # and the exception raised again
            raise
        
        xoff, yoff = info.getPixColRow(0, 0)
        blockcols, blockrows = info.getBlockSize()
        margin = info.getOverlapSize()
        block = ImageReader.readBlockWithMargin(vector.rasterDS, xoff, yoff, blockcols, blockrows, 
            vector.datatype, margin, [vector.nullval])

        return block
    
    def rasterize(self, info):
        """
        Rasterize the container of Vector objects passed to the 
        constuctor. Returns blocks in the same form as the 
        container passed to the constructor.

        """
        if isinstance(self.vectorContainer, dict):
            blockContainer = {}
            for key in self.vectorContainer:
                vector = self.vectorContainer[key]
                if isinstance(vector, list):
                    block = [self.rasterizeSingle(info, v, self.progress) for v in vector]
                else:
                    block = self.rasterizeSingle(info, vector, self.progress)
                blockContainer[key] = block

        elif isinstance(self.vectorContainer, Vector):
            blockContainer = self.rasterizeSingle(info, self.vectorContainer, self.progress)
    
        else:
            blockContainer = []
            for vector in self.vectorContainer:
                block = self.rasterizeSingle(info, vector, self.progress)
                blockContainer.append(block)

        return blockContainer
        
    def close(self):
        """
        Closes all datasets and removes temporary files.

        """
        if isinstance(self.vectorContainer, dict):
            for key in self.vectorContainer:
                vector = self.vectorContainer[key]
                vector.cleanup()

        elif isinstance(self.vectorContainer, Vector):
            self.vectorContainer.cleanup()

        else:
            for vector in self.vectorContainer:
                vector.cleanup()
