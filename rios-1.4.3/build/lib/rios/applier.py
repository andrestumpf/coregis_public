#!/usr/bin/env python
"""
Basic tools for setting up a function to be applied over 
a raster processing chain. The :func:`rios.applier.apply` function is the main
point of entry in this module. 

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

import numpy
from osgeo import gdal
from osgeo import ogr

from . import imagereader
from . import imagewriter
from . import imageio
from . import rioserrors
from . import vectorreader
from . import cuiprogress
from . import calcstats
from . import rat
from .parallel import jobmanager

# All default values, etc., copied in from their appropriate rios modules. 
DEFAULT_RESAMPLEMETHOD = "near"
DEFAULTFOOTPRINT = imagereader.DEFAULTFOOTPRINT
DEFAULTWINDOWXSIZE = imagereader.DEFAULTWINDOWXSIZE
DEFAULTWINDOWYSIZE = imagereader.DEFAULTWINDOWYSIZE
DEFAULTOVERLAP = imagereader.DEFAULTOVERLAP
DEFAULTLOGGINGSTREAM = imagereader.DEFAULTLOGGINGSTREAM
DEFAULTDRIVERNAME = imagewriter.DEFAULTDRIVERNAME
# Use of DEFAULTCREATIONOPTIONS is now deprecated, in favour of dfltDriverOptions[driver]
DEFAULTCREATIONOPTIONS = imagewriter.DEFAULTCREATIONOPTIONS
dfltDriverOptions = imagewriter.dfltDriverOptions
INTERSECTION = imageio.INTERSECTION
"Use the spatial intersection of inputs"
UNION = imageio.UNION
"Use the spatial union of inputs"
BOUNDS_FROM_REFERENCE = imageio.BOUNDS_FROM_REFERENCE
"Use the spatial extent of the reference file"

# From the calcstats module
DEFAULT_OVERVIEWLEVELS = calcstats.DEFAULT_OVERVIEWLEVELS
"Default overview levels on output images"
DEFAULT_MINOVERVIEWDIM = calcstats.DEFAULT_MINOVERVIEWDIM
"Default minimum dimension of highest overview level calculated"
DEFAULT_OVERVIEWAGGREGRATIONTYPE = calcstats.DEFAULT_OVERVIEWAGGREGRATIONTYPE
"Default agregation type when using formats without LAYER_TYPE"

# For supporting that automatic color table thing which Sam loves
DEFAULT_AUTOCOLORTABLETYPE = rat.DEFAULT_AUTOCOLORTABLETYPE


if sys.version_info[0] > 2:
    # hack for Python 3 which uses str instead of basestring
    # we just use basestring
    basestring = str

class FilenameAssociations(object): 
    """
    Class for associating external image filenames with internal
    names, which are then the same names used inside a function given
    to the :func:`rios.applier.apply` function. 
    
    Each attribute created on this object should be a filename, or a 
    list of filenames. The corresponding attribute names will appear 
    on the 'inputs' or 'outputs' objects inside the applied function. 
    Each such attribute will be an image data block or a list of image 
    data blocks, accordingly. 
    
    """
    def __len__(self):
        "Number of names defined on this instance (a list counts as only one name)"
        return len(self.__dict__.keys())

class BlockAssociations(object): 
    """
    Generic object to store the image blocks used within
    the applied function. The attributes are named the same way 
    as in the corresponding FilenameAssociations object, but are
    blocks of image data, instead of filenames. Where lists of 
    filenames were used, lists of image blocks are used here. 
    """
    pass
    
class OtherInputs(object): 
    """
    Generic object to store any extra inputs and outputs used 
    inside the function being applied. This class was originally
    named for inputs, but in fact works just as well for outputs, 
    too. Any items stored on this will be persistent between 
    iterations of the block loop. 
    """
    pass

class ApplierControls(object):
    """
    Controls for the operation of rios, for use with 
    the :func:`rios.applier.apply` function. 
    
    This object starts with default values for all controls, and 
    has methods for setting each of them to something else. 
    
    Attributes are:
        * **windowxsize**     X size of rios block (pixels)
        * **windowysize**     Y size of rios block (pixels)
        * **overlap**         Number of pixels in margin for block overlaps
        * **footprint**       :data:`rios.applier.INTERSECTION` or :data:`rios.applier.UNION` or :data:`rios.applier.BOUNDS_FROM_REFERENCE`
        * **drivername**      GDAL driver short name for output
        * **creationoptions** GDAL creation options for output
        * **thematic**        True/False for thematic outputs
        * **layernames**      List of layer names for outputs
        * **referenceImage**  Image for reference projection and grid that inputs will be resampled to.
        * **referencePixgrid** pixelGrid for reference projection and grid
        * **loggingstream**   file-like for logging of messages
        * **progress**        progress object
        * **statsIgnore**     global stats ignore value for output (i.e. null value)
        * **statscache**      stats cache if pre-calculated
        * **calcStats**       True/False to signal calculate statistics and pyramids
        * **omitPyramids**    True/False to omit pyramids when doing stats
        * **overviewLevels**  List of level factors used when calculating output image overviews
        * **overviewMinDim**  Minimum dimension of highest overview level
        * **overviewAggType** Aggregation type for calculating overviews
        * **tempdir**         Name of directory for temp files (resampling, etc.)
        * **resampleMethod**  String for resample method, when required (as per GDAL)
        * **numThreads**      Number of parallel threads used for processing each image block
        * **jobManagerType**  Which :class:`rios.parallel.jobmanager.JobManager` sub-class to use for parallel processing (by name)
        * **autoColorTableType** Type of color table to be automatically added to thematic output rasters
    
    Options relating to vector input files
        * **burnvalue**       Value to burn into raster from vector
        * **filtersql**       SQL where clause used to filter vector features
        * **alltouched**      Boolean. If True, all pixels touched are included in vector. 
        * **burnattribute**   Name of vector attribute used to supply burnvalue
        * **vectorlayer**     Number (or name) of vector layer
        * **vectordatatype**  Numpy datatype to use for raster created from vector
        * **vectornull**      Rasterised vector is initialised to this value, before burning
        
    
    Default values are provided for all attributes, and can then be over-ridden
    with the 'set' methods given. 
    
    Some 'set' methods take an optional imagename argument. If given, this should be 
    the same internal name used for a given image as in the :class:`rios.applier.FilenameAssociations`
    objects. This is the internal name for that image, and the method will set 
    the parameter in question for that specific image, which will over-ride the
    global value set when no imagename is given. 
    
    """
    def __init__(self):
        self.loggingstream = sys.stdout
        self.drivername = DEFAULTDRIVERNAME
        self.overlap = DEFAULTOVERLAP
        self.windowxsize = DEFAULTWINDOWXSIZE
        self.windowysize = DEFAULTWINDOWYSIZE
        self.footprint = DEFAULTFOOTPRINT
        self.referenceImage = None
        self.referencePixgrid = None
        self.progress = cuiprogress.SilentProgress()
        self.creationoptions = None
        self.statscache = None
        self.statsIgnore = 0
        self.calcStats = True
        self.omitPyramids = False
        self.overviewLevels = DEFAULT_OVERVIEWLEVELS
        self.overviewMinDim = DEFAULT_MINOVERVIEWDIM
        self.overviewAggType = None
        self.thematic = False
        self.layernames = None
        self.tempdir = '.'
        self.resampleMethod = DEFAULT_RESAMPLEMETHOD
        self.numThreads = 1
        self.jobManagerType = os.getenv('RIOS_DFLT_JOBMGRTYPE', default=None)
        self.autoColorTableType = DEFAULT_AUTOCOLORTABLETYPE

        # Vector fields
        self.burnvalue = 1
        self.vectornull = 0
        self.burnattribute = None
        self.filtersql = None
        self.alltouched = False
        self.vectordatatype = numpy.uint8
        self.vectorlayer = 0
        self.layerselection = None
        
        
        # Options specific to a named image. This was added on later, and is 
        # only valid for some of the attributes, so it looks a bit out-of-place.
        # Instead of the options being attributes of self, they are keys in a
        # dictionary. This dictionary is managed by the two methods
        # setOptionForImagename() and getOptionForImagename(). 
        self.optionsByImage = {}
    
    def setOptionForImagename(self, option, imagename, value):
        """
        Set the given option specifically for the given imagename. This 
        method is for internal use only. If you wish to set a particular 
        attribute, use the corresponding 'set' method. 
        """
        if imagename is None:
            setattr(self, option, value)
        else:
            if not option in self.optionsByImage:
                self.optionsByImage[option] = {}
            self.optionsByImage[option][imagename] = value
            
    def getOptionForImagename(self, option, imagename):
        """
        Returns the value of a particular option for the 
        given imagename. If only the global option has been set,
        then that is returned, but if a specific value has been set for 
        the given imagename, then use that. 
        
        The imagename is the same internal name as used for the image
        in the :class:`rios.applier.FilenameAssociations` objects. 
        
        """
        value = getattr(self, option)
        if option in self.optionsByImage:
            if imagename in self.optionsByImage[option]:
                value = self.optionsByImage[option][imagename]
        return value
        
    def setLoggingStream(self, loggingstream):
        """
        Set the rios logging stream to the given file-like object. 
        """
        self.loggingstream = loggingstream
        
    def setOverlap(self, overlap):
        """
        Set the overlap to the given value.
        Overlap is a number of pixels, and is somewhat mis-named. It refers 
        to the amount of margin added to each block of input, so that the blocks 
        will overlap, hence the actual amount of overlap is really more like
        double this value (allowing for odd and even numbers, etc). 
        """
        self.overlap = overlap
        
    def setOutputDriverName(self, drivername, imagename=None):
        """
        Set the output driver name to the given GDAL shortname.
        
        Note that the GDAL creation options have defaults suitable only 
        for the default driver, so if one sets the output driver, then 
        the creation options should be reviewed too. 
        
        In more recent versions of RIOS, the addition of driver-specific
        default creation options ($RIOS_DFLT_CREOPT_<driver>) allows for
        multiple default creation optionsto be set up. 
        
        """
        self.setOptionForImagename('drivername', imagename, drivername)
        
    def setWindowXsize(self, windowxsize):
        """
        Set the X size of the blocks used. Images are processed in 
        blocks (windows) of 'windowxsize' columns, and 'windowysize' rows. 
        """
        self.windowxsize = windowxsize
        
    def setWindowYsize(self, windowysize):
        """
        Set the Y size of the blocks used. Images are processed in 
        blocks (windows) of 'windowxsize' columns, and 'windowysize' rows. 
        """
        self.windowysize = windowysize
        
    def setFootprintType(self, footprint):
        """
        Set type of footprint, one of INTERSECTION, UNION or 
        BOUNDS_FROM_REFERENCE from this module
        """
        self.footprint = footprint
        
    def setReferenceImage(self, referenceImage):
        """
        Set the name of the image to use for the reference pixel grid and 
        projection. If neither referenceImage nor referencePixgrid are set, 
        then no resampling will be allowed. Only set one of referenceImage or
        referencePixgrid. 
        
        Note that this is the external filename, not the internal name (which 
        unfortunately is a bit inconsistent with everything else). 
        
        """
        self.referenceImage = referenceImage
        
    def setReferencePixgrid(self, referencePixgrid):
        """
        Set the reference pixel grid. If neither referenceImage nor 
        referencePixgrid are set, then no resampling will be allowed. 
        Only set one of referenceImage or referencePixgrid. The referencePixgrid
        argument is of type :class:`rios.pixelgrid.PixelGridDefn`. 

        """
        self.referencePixgrid = referencePixgrid

    def setProgress(self, progress):
        """
        Set the progress display object. Default is no progress
        object. 
        """
        self.progress = progress
        
    def setCreationOptions(self, creationoptions, imagename=None):
        """
        Set a list of GDAL creation options (should match with the driver). 
        Each list element is a string of the form "NAME=VALUE". 
        
        Defaults are suitable for the default driver, and need to be changed
        if that is changed. However, if an appropriate driver-specific default 
        environment variable ($RIOS_DFLT_CREOPT_<driver>) is given, this 
        will be used. 
        
        """
        self.setOptionForImagename('creationoptions', imagename, creationoptions)
        
    def setStatsCache(self, statscache, imagename=None):
        "Set the stats cache, if statistics are known from some other source."
        self.setOptionForImagename('statscache', imagename, statscache)
        
    def setStatsIgnore(self, statsIgnore, imagename=None):
        """
        Set the global default value to use as the 
        null value when calculating stats.
        """
        self.setOptionForImagename('statsIgnore', imagename, statsIgnore)
        
    def setCalcStats(self, calcStats, imagename=None):
        """
        Set True to calc stats, False otherwise. If True, then statistics and 
        pyramid layers are calculated (if supported by the driver
        """
        self.setOptionForImagename('calcStats', imagename, calcStats)
        
    def setOmitPyramids(self, omitPyramids, imagename=None):
        """
        Set True to omit pyramid layers, False otherwise. If True, then when
        statistics are being calculated, pyramid layers will be omitted, 
        otherwise they will be created at the same time. 
        Usual default is False. 
        """
        self.setOptionForImagename('omitPyramids', imagename, omitPyramids)
    
    def setOverviewLevels(self, overviewLevels, imagename=None):
        """
        Set the overview levels to be used on output images (i.e. pyramid layers). 
        Levels are specified as a list of integer factors, with the same meanings 
        as given to the gdaladdo command. 
        
        """
        self.setOptionForImagename('overviewLevels', imagename, overviewLevels)
    
    def setOverviewMinDim(self, overviewMinDim, imagename=None):
        """
        Set minimum dimension allowed for output overview. Overview levels (i.e. pyramid
        layers) will be calculated as per the overviewLevels list of factors, but 
        only until the minimum dimension falls below the value of overviewMinDim
        
        """
        self.setOptionForImagename('overviewMinDim', imagename, overviewMinDim)
    
    def setOverviewAggregationType(self, overviewAggType, imagename=None):
        """
        Set the type of aggregation used when computing overview images (i.e. pyramid 
        layers). Normally a thematic image should be aggregated using "NEAREST", while a 
        continuous image should be aggregated using "AVERAGE". When the setting is 
        given as None, then a default is used. If using an output format which 
        supports LAYER_TYPE, the default is based on this, but if not, it comes from 
        the value of the environment variable $RIOS_DEFAULT_OVERVIEWAGGREGATIONTYPE.
        
        This method should usually be used to set when writing an output to a format
        which does not support LAYER_TYPE, and which is not appropriate for the
        setting given by the environment default. 
        
        """
        self.setOptionForImagename('overviewAggType', imagename, overviewAggType)
        
    def setThematic(self, thematicFlag, imagename=None):
        "Set boolean value of thematic flag (may not be supported by the GDAL driver)"
        self.setOptionForImagename('thematic', imagename, thematicFlag)

    def setLayerNames(self, layerNames, imagename=None):
        """
        Set list of layernames to be given to the output file(s)
        """
        self.setOptionForImagename('layernames', imagename, layerNames)
        
    def setTempdir(self, tempdir):
        "Set directory to use for temporary files for resampling, etc. "
        self.tempdir = tempdir
        
    def setResampleMethod(self, resampleMethod, imagename=None):
        """
        Set resample method to be used for all resampling. Possible 
        options are those defined by gdalwarp, i.e. 'near', 'bilinear', 
        'cubic', 'cubicspline', 'lanczos'. 
        """
        self.setOptionForImagename('resampleMethod', imagename, resampleMethod)
    
    def setBurnValue(self, burnvalue, vectorname=None):
        """
        Set the burn value to be used when rasterizing the input vector(s).
        If vectorname given, set only for that vector. Default is 1. 
        """
        self.setOptionForImagename('burnvalue', vectorname, burnvalue)
    
    def setBurnAttribute(self, burnattribute, vectorname=None):
        """
        Set the vector attribute name from which to get the burn value
        for each vector feature. If vectorname is given, set only for that
        vector input. Default is to use burnvalue instead of burnattribute. 
        """
        self.setOptionForImagename('burnattribute', vectorname, burnattribute)
    
    def setVectorNull(self, vectornull, vectorname=None):
        """
        Set the vector null value. This is used to initialise the
        rasterised vector, before burning in the burn value. This is of most
        importance when burning values from a vector attribute column, as 
        this should be a distinct value from any of the values in the column. 
        If this is not so, then polygons can end up blending with the background,
        resulting in incorrect answers. 
        """
        self.setOptionForImagename('vectornull', vectorname, vectornull)
    
    def setFilterSQL(self, filtersql, vectorname=None):
        """
        Set an SQL WHERE clause which will be used to filter vector features.
        If vectorname is given, then set only for that vector
        """
        self.setOptionForImagename('filtersql', vectorname, filtersql)
    
    def setAlltouched(self, alltouched, vectorname=None):
        """
        Set boolean value of alltouched attribute. If alltouched is True, then
        pixels will count as "inside" a vector polygon if they touch the polygon,
        rather than only if their centre is inside. 
        If vectornmame given, then set only for that vector. 
        """
        self.setOptionForImagename('alltouched', vectorname, alltouched)
    
    def setVectorDatatype(self, vectordatatype, vectorname=None):
        """
        Set numpy datatype to use for rasterized vectors
        If vectorname given, set only for that vector
        """
        self.setOptionForImagename('vectordatatype', vectorname, vectordatatype)
    
    def setVectorlayer(self, vectorlayer, vectorname=None):
        """
        Set number/name of vector layer, for vector formats which have 
        multiple layers. Not required for plain shapefiles. 
        Can be either a layer number (start at zero) or 
        a layer name. If vectorname given, set only for that vector.
        """
        self.setOptionForImagename('vectorlayer', vectorname, vectorlayer)
    
    def makeResampleDict(self, imageDict):
        """
        Make a dictionary of resample methods, one for every image
        name in the given dictionary. This method is for internal use only. 
        """
        d = {}
        imagenamelist = imageDict.keys()
        for name in imagenamelist:
            method = self.getOptionForImagename('resampleMethod', name)
            if isinstance(imageDict[name], list):
                # We have a list of images for this name, so make a list of 
                # resample methods
                d[name] = [method] * len(imageDict[name])
            else:
                # We have just one image, so the corresponding entry is just one 
                # resample method
                d[name] = method
        return d
    
    def selectInputImageLayers(self, layerselection, imagename=None):
        """
        Set which layers are to be read from the input image(s). Default
        will read all layers. If imagename is given, selection will be for 
        that image only. The layerselection parameter should be a list
        of layer numbers. Layer numbers follow GDAL conventions, i.e. 
        a layer number of 1 refers to the first layer in the file. 
        Can  be much more efficient when only using a small subset of 
        layers from the inputs. 
        """
        self.setOptionForImagename('layerselection', imagename, layerselection)
    
    def setNumThreads(self, numThreads):
        """
        Set the number of 'threads' to be used when processing each block 
        of imagery. Note that these are not threads in the technical sense, 
        but are handled by the JobManager class, and are some form of 
        cooperating parallel processes, depending on the type of job 
        manager sub-class selected. See :mod:`rios.parallel.jobmanager` 
        for full details. Note that this is only worth using on very 
        computationally-intensive tasks. Default is 1, i.e. no parallel 
        processing. 
        
        """
        self.numThreads = numThreads
    
    def setJobManagerType(self, jobMgrType):
        """
        Set which type of JobManager is to be used for parallel processing.
        See :mod:`rios.parallel.jobmanager` for details. Default is taken from
        $RIOS_DFLT_JOBMGRTYPE. 
        
        """
        self.jobManagerType = jobMgrType
    
    def setAutoColorTableType(self, autoColorTableType, imagename=None):
        """
        If this option is set, then thematic raster outputs will have a
        color table automatically generated and attached to them. The type is
        passed to :func:`rios.rat.genColorTable` to determine what type of automatic
        color table is generated. 
        
        The default type will be taken from $RIOS_DFLT_AUTOCOLORTABLETYPE if it
        is set. If that is not set, then the default is not to automatically attached
        any color table to thematic output rasters.
        
        In practise, it is probably simpler to explicitly set the color table using 
        the :func:`rios.rat.setColorTable` function, after creating the file, but this
        is an alternative. 
        
        Note that the imagename parameter can be given, in which case the autoColorTableType 
        will only be applied to that raster. 
        
        None of this has any impact on athematic outputs. 
        
        """
        self.setOptionForImagename('autoColorTableType', imagename, autoColorTableType)


def apply(userFunction, infiles, outfiles, otherArgs=None, controls=None):
        """
        Apply the given 'userFunction' to the given
        input and output files. 
        
        infiles and outfiles are :class:`rios.applier.FilenameAssociations` objects to 
        define associations between internal variable names and
        external filenames, for the raster file inputs and outputs. 
        
        otherArgs is an object of extra arguments to be passed to the 
        userFunction, each with a sensible name on the object. These 
        can be either input or output arguments, entirely at the discretion
        of userFunction(). otherArgs should be in instance of :class:`rios.applier.OtherInputs`
        
        The userFunction has the following call sequence::

            userFunction(info, inputs, outputs)

        or::

            userFunction(info, inputs, outputs, otherArgs)

        if otherArgs is not None. 
        inputs and outputs are objects in which there are named attributes 
        with the same names as those given in the infiles and outfiles 
        objects. In the inputs and outputs objects, available inside 
        userFunction, these attributes contain numpy arrays of data read 
        from/written to the corresponding image file. 
        
        If the attributes given in the infiles or outfiles objects are 
        lists of filenames, the the corresponding attributes of the 
        inputs and outputs objects inside the applied function will be 
        lists of image data blocks instead of single blocks. 
        
        The numpy arrays are always 3-d arrays, with shape::

            (numBands, numRows, numCols)

        The datatype of the output image(s) is determined directly
        from the datatype of the numpy arrays in the outputs object. 
        
        The info object contains many useful details about the processing, 
        and will always be passed to the userFunction. It can, of course, 
        be ignored. It is an instance of the :class:`rios.readerinfo.ReaderInfo` class. 
        
        The controls argument, if given, is an instance of the 
        :class:`rios.applier.ApplierControls` class, which allows control of various 
        aspects of the reading and writing of images. See the class 
        documentation for further details. 

        There is a page dedicated to :doc:`applierexamples`.

        """
        # Get default controls object if none given. 
        if controls is None:
            controls = ApplierControls()
        
        (imagefiles, vectorfiles) = separateVectors(infiles)
        inputImageLayerSelection = makeInputImageLayerSelection(imagefiles, controls)
        reader = imagereader.ImageReader(imagefiles.__dict__, 
            controls.footprint, controls.windowxsize, controls.windowysize, 
            controls.overlap, controls.statscache, loggingstream=controls.loggingstream,
            layerselection=inputImageLayerSelection)

        vecreader = None
        if len(vectorfiles) > 0:
            vectordict = makeVectorObjects(vectorfiles, controls)
            vecreader = vectorreader.VectorReader(vectordict, progress=controls.progress)
        
        handleInputResampling(imagefiles, controls, reader)

        writerdict = {}
        
        if controls.progress is not None:
            controls.progress.setTotalSteps(100)
            controls.progress.setProgress(0)
        lastpercent = 0
        
        # Set up for parallel processing, if requested. 
        jobmgr = None
        if controls.numThreads > 1:
            jobmgr = jobmanager.getJobMgrObject(controls)

        done = False
        iterator = reader.__iter__()
        while not done:
            # list of RIOSJobInfo
            jobInputs = []

            for n in range(controls.numThreads):
                try:
                    info, blockdict = iterator.__next__()
                except StopIteration:
                    done = True
                    break

                inputBlocks = BlockAssociations()
                inputBlocks.__dict__.update(blockdict)
                if vecreader is not None:
                    vecblocks = vecreader.rasterize(info)
                    inputBlocks.__dict__.update(vecblocks)

                # build a RIOSJobInfo with the params
                jobInfo = RIOSJobInfo(info, inputBlocks, otherArgs)
                jobInputs.append(jobInfo)

            if len(jobInputs) == 0:
                break
            
            # Now call the function with those args
            if jobmgr is None:
                # single threaded - just call it
                params = jobInfo.getFunctionParams()
                userFunction(*params)
                outputBlocks = jobInfo.getFunctionResult(params)
                writeOutputBlocks(writerdict, outfiles, outputBlocks, 
                                controls, info)
            else:
                # multi threaded - get the job manager to run jobs
                outBlocksList = jobmgr.runSubJobs(userFunction, jobInputs)
                for outputBlocks in outBlocksList:
                    writeOutputBlocks(writerdict, outfiles, outputBlocks, 
                                controls, info)
            
            lastpercent = updateProgress(controls, info, lastpercent)
                
        if controls.progress is not None:
            controls.progress.setProgress(100)

        closeOutputImages(writerdict, outfiles, controls)


def closeOutputImages(writerdict, outfiles, controls):
    """
    Called by :func:`rios.applier.apply` to close all output image files. 
    """
    for name in outfiles.__dict__.keys():
        writer = writerdict[name]
        if isinstance(writer, list):
            writerList = writer
        else:
            writerList = [writer]

        for singleWriter in writerList:
            singleWriter.close(calcStats=controls.getOptionForImagename('calcStats', name), 
                statsIgnore=controls.getOptionForImagename('statsIgnore', name), 
                progress=controls.progress,
                omitPyramids=controls.getOptionForImagename('omitPyramids', name),
                overviewLevels=controls.getOptionForImagename('overviewLevels', name),
                overviewMinDim=controls.getOptionForImagename('overviewMinDim', name), 
                overviewAggType=controls.getOptionForImagename('overviewAggType', name),
                autoColorTableType=controls.getOptionForImagename('autoColorTableType', name))


def updateProgress(controls, info, lastpercent):
    """
    Called by :func:`rios.applier.apply` to update progress
    """
    if controls.progress is not None:
        percent = info.getPercent()
        if percent != lastpercent:
            controls.progress.setProgress(percent)
            lastpercent = percent
    return lastpercent


def handleInputResampling(infiles, controls, reader):
    """
    Called by :func:`rios.applier.apply` to handle automatic resampling of input rasters.
    Most of the work is done by the read.allowResample() method. 
    """
    if controls.referenceImage is not None:
        resampleDict = controls.makeResampleDict(infiles.__dict__)
        reader.allowResample(refpath=controls.referenceImage, tempdir=controls.tempdir,
            resamplemethod=resampleDict, useVRT=True)
    elif controls.referencePixgrid is not None:
        resampleDict = controls.makeResampleDict(infiles.__dict__)
        reader.allowResample(refPixgrid=controls.referencePixgrid, 
            tempdir=controls.tempdir, 
            resamplemethod=resampleDict, useVRT=True)


def writeOutputBlocks(writerdict, outfiles, outputBlocks, controls, info):
    """
    Called by :func:`rios.applier.apply`, to write the output blocks, after
    they have been created by the user function. 
    For internal use only. 
    
    For all names given in outfiles object, look for a data block 
    of the same name in the outputBlocks object. If the given name
    is a list, then the corresponding name should be a list of blocks. 
    
    """
    for name in outfiles.__dict__.keys():
        if name not in outputBlocks.__dict__:
            msg = 'Output key %s not found in output blocks' % name
            raise rioserrors.KeysMismatch(msg)

        outblock = outputBlocks.__dict__[name]
        outfileName = getattr(outfiles, name)
        if name not in writerdict:
            # We have not yet created the output writers
            if isinstance(outfileName, list):
                # We have a list of filenames under this name in the dictionary,
                # and so we must create a list of writers. The outblock will also be 
                # a list of blocks
                writerdict[name] = []
                numFiles = len(outfileName)
                if len(outblock) != numFiles:
                    raise rioserrors.MismatchedListLengthsError(("Output '%s' writes %d files, "+
                        "but only %d blocks given")%(name, numFiles, len(outblock)))
                for i in range(numFiles):
                    filename = outfileName[i]
                    writer = imagewriter.ImageWriter(filename, info=info, 
                        firstblock=outblock[i], 
                        drivername=controls.getOptionForImagename('drivername', name), 
                        creationoptions=controls.getOptionForImagename('creationoptions', name))
                    writerdict[name].append(writer)
                    if controls.getOptionForImagename('thematic', name):
                        writer.setThematic()

                    layernames = controls.getOptionForImagename('layernames', name)
                    if layernames is not None:
                        writer.setLayerNames(layernames)
            else:
                # This name in the dictionary is just a single filename
                writer = imagewriter.ImageWriter(outfileName, info=info, firstblock=outblock,
                    drivername=controls.getOptionForImagename('drivername', name), 
                    creationoptions=controls.getOptionForImagename('creationoptions', name))
                writerdict[name] = writer
                if controls.getOptionForImagename('thematic', name):
                    writer.setThematic()

                layernames = controls.getOptionForImagename('layernames', name)
                if layernames is not None:
                    writer.setLayerNames(layernames)
        else:
            # The output writers exist, so select the correct one and write the block
            if isinstance(outfileName, list):
                # We have a list of files for this name, and a list of blocks to write
                numFiles = len(outfileName)
                if len(outblock) != numFiles:
                    raise rioserrors.MismatchedListLengthsError(("Output '%s' writes %d files, "+
                        "but only %d blocks given")%(name, numFiles, len(outblock)))
                for i in range(numFiles):
                    writerdict[name][i].write(outblock[i])
            else:
                # This name is just a single file, and we write a single block
                writerdict[name].write(outblock)


def separateVectors(infiles):
    """
    Given a :class:`rios.applier.FilenameAssociations` object, separate out the files which 
    are raster, and the files which are vectors. Returns two :class:`rios.applier.FilenameAssociations`
    objects, carrying the same attribute names, but each has only the raster
    or the vectors. 
    
    """
    imagefiles = FilenameAssociations()
    vectorfiles = FilenameAssociations()
    
    nameList = sorted(infiles.__dict__.keys())
    for name in nameList:
        fileValue = getattr(infiles, name)
        if isinstance(fileValue, basestring):
            testFilename = fileValue
        elif isinstance(fileValue, list):
            # We only check the first filename in a list. If the user
            # mixed rasters and vectors in one list, things would go horribly wrong
            testFilename = fileValue[0]
        else:
            testFilename = None

        if opensAsRaster(testFilename):
            setattr(imagefiles, name, fileValue)
        elif opensAsVector(testFilename):
            setattr(vectorfiles, name, fileValue)
        else:
            raise rioserrors.FileOpenError("Failed to open file '%s' as either raster or vector"%testFilename)
        
    return (imagefiles, vectorfiles)


def opensAsRaster(filename):
    """
    Return True if filename opens as a GDAL raster, False otherwise
    """
    usingExceptions = False
    if hasattr(gdal, 'GetUseExceptions'):
        usingExceptions = gdal.GetUseExceptions()
    gdal.UseExceptions()
    try:
        ds = gdal.Open(filename)
    except Exception:
        ds = None
    opensOK = (ds is not None)
    
    if not usingExceptions:
        gdal.DontUseExceptions()
    return opensOK


def opensAsVector(filename):
    """
    Return True if filename opens as an OGR vector, False otherwise
    """
    usingExceptions = False
    if hasattr(ogr, 'GetUseExceptions'):
        usingExceptions = ogr.GetUseExceptions()
    ogr.UseExceptions()
    try:
        ds = ogr.Open(filename)
    except Exception:
        ds = None
    opensOK = (ds is not None)
    
    if not usingExceptions:
        ogr.DontUseExceptions()
    return opensOK


def makeVectorObjects(vectorfiles, controls):
    """
    Returns a dictionary of :class:`rios.vectorreader.Vector` objects,
    with the keys being the attribute names used
    on the vectorfiles object. This is then ready to
    go into the :class:`rios.vectorreader.VectorReader` constructor. 
    
    """
    vectordict = {}
    namelist = sorted(vectorfiles.__dict__.keys())
    for name in namelist:
        burnvalue = controls.getOptionForImagename('burnvalue', name)
        vectordatatype = controls.getOptionForImagename('vectordatatype', name)
        alltouched = controls.getOptionForImagename('alltouched', name)
        vectorlayer = controls.getOptionForImagename('vectorlayer', name)
        burnattribute = controls.getOptionForImagename('burnattribute', name)
        filtersql = controls.getOptionForImagename('filtersql', name)
        tempdir = controls.tempdir
        vectornull = controls.getOptionForImagename('vectornull', name)
        
        fileValue = getattr(vectorfiles, name)
        if isinstance(fileValue, list):
            veclist = []
            for filename in fileValue:
                vec = vectorreader.Vector(filename, burnvalue=burnvalue, datatype=vectordatatype,
                    attribute=burnattribute, filter=filtersql, inputlayer=vectorlayer,
                    alltouched=alltouched, tempdir=tempdir, nullval=vectornull)
                veclist.append(vec)
            vectordict[name] = veclist
        elif isinstance(fileValue, basestring):
            vectordict[name] = vectorreader.Vector(fileValue, burnvalue=burnvalue, 
                datatype=vectordatatype, attribute=burnattribute, filter=filtersql, 
                inputlayer=vectorlayer, alltouched=alltouched, tempdir=tempdir, 
                nullval=vectornull)

    return vectordict


def makeInputImageLayerSelection(imagefiles, controls):
    """
    Make a dictionary with the same image name keys as imagefiles, but with
    layerselection lists for each entry, as per the controls object. If only
    some images have a layerselection set, then the remaining entries are None. 
    
    """
    layerselection = {}
    for name in imagefiles.__dict__.keys():
        layerselection[name] = controls.getOptionForImagename('layerselection', name)
    return layerselection

class RIOSJobInfo(jobmanager.JobInfo):
    """
    Class that contains information for parameters to a RIOS
    function
    """
    def __init__(self, info, inputs, otherargs=None):
        self.info = info
        self.inputs = inputs
        self.otherargs = otherargs
        # we don't bother pickling the outputs - start again 
        # with a fresh BlockAssociations

    def prepareForPickling(self):
        """
        GDAL datasets cannot be pickled
        and neither stderr (in info.logginstream)
        so we clean up the info object a bit

        """
        self.info.blocklookup = {}
        self.info.loggingstream = None
        return self

    def getFunctionParams(self):
        """
        Return the parameters as a tuple.

        When this is run in the subprocess, we can 
        reset the info.loggingstream to something
        other than None

        """
        if self.info.loggingstream is None:
            self.info.loggingstream = imagereader.DEFAULTLOGGINGSTREAM
        outputs = BlockAssociations()
        params = (self.info, self.inputs, outputs)
        if self.otherargs is not None:
            params += (self.otherargs,)

        return params

    def getFunctionResult(self, params):
        """
        Return the ouputs parameter

        """
        return params[2]

