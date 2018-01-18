
"""
Contains the ImageReader class

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
import copy
import numpy
from osgeo import gdal
from . import imageio
from . import inputcollection
from . import readerinfo
from . import rioserrors

if sys.version_info[0] > 2:
    # hack for Python 3 which uses str instead of basestring
    # we just use basestring
    basestring = str

DEFAULTFOOTPRINT = int(os.getenv('RIOS_DFLT_FOOTPRINT', 
                            default=imageio.INTERSECTION))
DEFAULTWINDOWXSIZE = int(os.getenv('RIOS_DFLT_BLOCKXSIZE', default=256))
DEFAULTWINDOWYSIZE = int(os.getenv('RIOS_DFLT_BLOCKYSIZE', default=256))
DEFAULTOVERLAP = int(os.getenv('RIOS_DFLT_OVERLAP', default=0))
DEFAULTLOGGINGSTREAM = sys.stdout

class ImageIterator(object):
    """
    Class to allow iteration across an ImageReader instance.
    Do not instantiate this class directly - it is created
    by ImageReader.__iter__().
    
    See http://docs.python.org/library/stdtypes.html#typeiter
    for a description of how this works. There is another way,
    see: http://docs.python.org/reference/expressions.html#yieldexpr
    but it seemed too much like Windows 3.1 programming which
    scared me!
    
    Returns a tuple containing an ReaderInfo class, plus a numpy
    array for each iteration
    
    """
    def __init__(self,reader):
        # reader = an ImageReader instance
        self.reader = reader
        self.nblock = 0 # start at first block
        
    def __iter__(self):
        # For iteration support - just return self.
        return self

    def next(self):
        # for Python 2.x
        return self.__next__()
        
    def __next__(self):
        # for iteration support. Raises a StopIteration
        # if we have read beyond the end of the image
        try:
            # get ImageReader.readBlock() to do the work
            # this raises a OutsideImageBounds exception,
            # but the iteration protocol expects a 
            # StopIteration exception.
            returnTuple = self.reader.readBlock(self.nblock)
        except rioserrors.OutsideImageBoundsError:
            raise StopIteration()
            
        # look at the next block next time
        self.nblock += 1
        
        return returnTuple

class ImageReader(object):
    """
    Class that reads a single file, a list or dictionary of files and 
    iterates through them block by block

    **Example**
    
    ::
    
        import sys
        from rios.imagereader import ImageReader

        reader = ImageReader(sys.argv[1]) 
        for (info, block) in reader:     
            block2 = block * 2

    """
    def __init__(self, imageContainer,
				footprint=DEFAULTFOOTPRINT,
				windowxsize=DEFAULTWINDOWXSIZE, windowysize=DEFAULTWINDOWYSIZE,
				overlap=DEFAULTOVERLAP, statscache=None,
                loggingstream=sys.stdout, layerselection=None):
        """
        imageContainer is a filename or list or dictionary that contains
        the filenames of the images to be read.
        If a list is passed, a list of blocks is returned at 
        each iteration, if a dictionary a dictionary is
        returned at each iteration with the same keys.
        
        footprint can be either INTERSECTION, UNION or BOUNDS_FROM_REFERENCE
        
        windowxsize and windowysize specify the size
        of the block to be read at each iteration
        
        overlap specifies the number of pixels to overlap
        between each block
        
        statscache if specified, should be an instance of 
        readerinfo.StatisticsCache. If None, cache is
        created per instance of this class. If doing
        multiple reads on same datasets, consider having 
        a single instance of statscache between all instances
        of this class.
        
        Set loggingstream to a file like object if you wish
        logging about resampling to be sent somewhere else
        rather than stdout.
        
        layerselection, if given, should be of the same type as imageContainer, 
        that is, if imageContainer is a dictionary, then layerselection 
        should be a dictionary with the same keys, and if imageContainer 
        is a list, then layerselection should be a list of the same length. 
        The elements in layerselection should always be lists of layer numbers, 
        used to select only particular layers to read from the corresponding 
        input image. Layer numbers use GDAL conventions, i.e. start at 1. 
        Default reads all layers. 
        
        """

        # grab the imageContainer so we can always know what 
        # type of container they passed in
        self.imageContainer = imageContainer
      
        if isinstance(imageContainer,dict):
            # Convert the given imageContainer into a list suitable for
            # the standard InputCollection. 
            imageList = []
            self.layerselectionList = []
            for name in imageContainer.keys():
                filename = imageContainer[name]
                if isinstance(filename, list):
                    # We have actually been given a list of filenames, so tack then all on to the imageList
                    imageList.extend(filename)
                elif isinstance(filename, basestring):
                    # We just have a single filename
                    imageList.append(filename)
                else:
                    msg = "Dictionary must contain either lists or strings. Got '%s' instead" % type(filename)
                    raise rioserrors.ParameterError(msg)

                # Layer selection for this filename. 
                thisLayerSelection = None
                if layerselection is not None and name in layerselection:
                    thisLayerSelection = layerselection[name]

                if isinstance(filename, list):
                    self.layerselectionList.extend([thisLayerSelection for fn in filename])
                else:
                    self.layerselectionList.append(thisLayerSelection)

        
        elif isinstance(imageContainer,basestring):
            # they passed a string, just make a list out of it
            imageList = [imageContainer]
            if layerselection is not None:
                self.layerselectionList = [layerselection]
            else:
                self.layerselectionList = [None]
        else:
            # we hope they passed a tuple or list. Don't need to do much
            imageList = imageContainer
            if layerselection is not None:
                self.layerselectionList = layerselection
            else:
                self.layerselectionList = [None for fn in imageList]
        
        # create an InputCollection with our inputs
        self.inputs = inputcollection.InputCollection(imageList,loggingstream=loggingstream)
        
        # save the other vars
        self.footprint = footprint
        self.windowxsize = windowxsize
        self.windowysize = windowysize
        self.overlap = overlap
        self.statscache = statscache
        # just create a new instance of the AttributeTableCache
        # possibly this should be passed in like statscache so the
        # attributes can be cached between instances of ImageReader
        # but considering retrieving an attribute nowhere as expensive
        # as getting global statistics, probably overkill.
        self.ratcache = readerinfo.AttributeTableCache()
        self.loggingstream = loggingstream
        
        # these are None until prepare() is called
        self.workingGrid = None
        self.info = None

    def __len__(self):
        # see http://docs.python.org/reference/datamodel.html#emulating-container-types
        
        # need self.info to be created so run prepare()
        if self.info is None:
            self.prepare()

        # get the total number of blocks for image            
        (xtotalblocks,ytotalblocks) = self.info.getTotalBlocks()
        
        # return the total number of blocks as our len()
        return xtotalblocks * ytotalblocks
        
    def __getitem__(self,key):
        # see http://docs.python.org/reference/datamodel.html#emulating-container-types
        # for indexing, returns tuple from readBlock()

        # need self.info to be created so run prepare()
        if self.info is None:
            self.prepare()

        # if they have passed a negative block - count
        # back from the end           
        if key < 0:
            # get total number of blocks
            (xtotalblocks,ytotalblocks) = self.info.getTotalBlocks()
            # add the key (remember, its negative)
            key = (xtotalblocks * ytotalblocks) + key
            if key < 0:
                # still negative - not enough blocks
                raise KeyError()
        
        try:
            # get readBlock() to do the work
            # this raises a OutsideImageBounds exception,
            # but the container protocol expects a 
            # KeyError exception.
            returnTuple = self.readBlock(key)
        except rioserrors.OutsideImageBoundsError:
            raise KeyError()
            
        return returnTuple
            

    def __iter__(self):
        # see http://docs.python.org/reference/datamodel.html#emulating-container-types

        # need self.info to be created so run prepare()
        if self.info is None:
            self.prepare()

        # return in ImageIterator instance
        # with a reference to this object            
        return ImageIterator(self)

    def allowResample(self, resamplemethod="near", refpath=None, refgeotrans=None, 
            refproj=None, refNCols=None, refNRows=None, refPixgrid=None, 
            tempdir='.', useVRT=False):
        """
        By default, resampling is disabled (all datasets must
        match). Calling this enables it. 
        Either refgeotrans, refproj, refNCols and refNRows must be passed, 
        or refpath passed and the info read from that file.

        tempdir is the temporary directory where the resampling happens. By 
        default the current directory.

        resamplemethod is the method used - must be supported by gdalwarp.
        This can be a single string if all files are to be resampled by the
        same method, or a list or dictionary (to match what passed to the 
        constructor) contain the methods for each file.
        
        If resampling is needed it will happen before the call returns.
        
        """
        # set the reference in our InputCollection
        self.inputs.setReference(refpath, refgeotrans, refproj,
                refNCols, refNRows, refPixgrid)
             
        if isinstance(resamplemethod, basestring):
            # turn it into a list with the same method repeated
            resamplemethodlist = [resamplemethod] * len(self.inputs)
        elif isinstance(resamplemethod, dict):
            # dictionary - check they passed a dictionary to the constructor
            # and the keys match
            if not isinstance(self.imageContainer, dict):
                msg = 'Can only pass a dictionary if a dictionary passed to the constructor'
                raise rioserrors.ParameterError(msg)
            elif sorted(self.imageContainer.keys()) != sorted(resamplemethod.keys()):
                msg = ('Dictionary keys must match those passed to the constructor, '+
                    'constructor keys = %s, resample keys = %s') % (self.imageContainer.keys(),
                    resamplemethod.keys())
                raise rioserrors.ParameterError(msg)
            else:
                # create a list out of the dictionary in the same way as the constructor does
                resamplemethodlist = []
                for name in resamplemethod.keys():
                    method = resamplemethod[name]
                    if isinstance(method, list):
                        # We have actually been given a list of method, so tack then all on to the resamplemethodlist
                        resamplemethodlist.extend(method)
                    elif isinstance(method, basestring):
                        # We just have a single method
                        resamplemethodlist.append(method)
                    else:
                        msg = "Dictionary must contain either lists or strings. Got '%s' instead" % type(method)
                        raise rioserrors.ParameterError(msg)

        else:
            # we assume they have passed a list/tuple
            if len(resamplemethod) != len(self.inputs):
                msg = 'must pass correct number of resample methods'
                raise rioserrors.ParameterError(msg)
            resamplemethodlist = resamplemethod

        try:   
            # resample all in collection to reference
            self.inputs.resampleAllToReference(self.footprint, resamplemethodlist, tempdir, useVRT)
        finally:
            # if the user interrupted, then ensure all temp
            # files removed.
            self.inputs.cleanup()
        
    def prepare(self, workingGrid=None):
        """
        Prepare to read from images. These steps are not
        done in the constructor, but are done just before
        reading in case allowResample() is called which
        will resample the inputs.
        
        The pixelGrid instance to use as the working grid can
        be passed in case it is not to be derived from the images
        to be read or is different from that passed to allowResample
        """
    
        # if resampled has happened then they should all match
        if not self.inputs.checkAllMatch():
            msg = 'Inputs do not match - must enable resampling'
            raise rioserrors.ResampleNeededError(msg)
        
        if workingGrid is None:
            # set the working grid based on the footprint
            self.workingGrid = self.inputs.findWorkingRegion(self.footprint)
        else:
            # user supplied
            self.workingGrid = workingGrid
        
        # create a statscache if not passed to constructor.
        # Created once per dataset so stats
        # only have to be calculated once per image - it
        # returns cached value for subsequent calls.
        if self.statscache is None:
            self.statscache = readerinfo.StatisticsCache()
        
        # create a ReaderInfo class with the info it needs
        # a copy of this class is passed with each iteration
        self.info = readerinfo.ReaderInfo(self.workingGrid, self.statscache, self.ratcache,
                        self.windowxsize, self.windowysize, self.overlap, self.loggingstream)
        
    def readBlock(self,nblock):
        """
        Read a block. This is normally called from the
        __getitem__ method when this class is indexed, 
        or from the ImageIterator when this class is 
        being iterated through.
        
        A block is read from each image and returned
        in a tuple along with a ReaderInfo instance.
        
        nblock is a single index, and will be converted
        to row/column.
        
        """
        
        # need self.info to be created so run prepare()
        if self.info is None:
            self.prepare()
           
        # do a shallow copy of the ReaderInfo.
        # this copy will have the fields filled in
        # that relate to the whole image.
        # We will then fill in the fields that relate
        # to this block. 
        # This means that calls to read other blocks
        # wont clobber the per block info, and user 
        # writing back into the object wont stuff up
        # the system 
        # because it is a shallow copy, statscache should
        # still be pointing to a single object
        info = copy.copy(self.info)
        
        # get the size of the are we are to read
        (xsize,ysize) = info.getTotalSize()
        
        # get the number of blocks are to read
        (xtotalblocks,ytotalblocks) = info.getTotalBlocks()
        
        # check they asked for block is valid
        if nblock >= (xtotalblocks * ytotalblocks):
            raise rioserrors.OutsideImageBoundsError()
        
        # convert the block to row/column
        yblock = nblock // xtotalblocks
        xblock = nblock % xtotalblocks
        
        # set this back to our copy of the info object
        info.setBlockCount(xblock,yblock)
    
        # calculate the coords of this block in pixels
        xcoord = xblock * self.windowxsize
        ycoord = yblock * self.windowysize
        
        # convert this to world coords
        blocktl = imageio.pix2wld( info.getTransform(), xcoord, ycoord )

        # work out the bottom right coord for this block
        nBlockBottomX = (( xblock + 1 ) * self.windowxsize)
        nBlockBottomY = (( yblock + 1 ) * self.windowysize)
        
        # make adjuctment if we are at the edge of the image
        # and there are smaller blocks
        if nBlockBottomX > xsize:
          nBlockBottomX = xsize
        if nBlockBottomY > ysize:
          nBlockBottomY = ysize

        # work out the world coords for the bottom right
        blockbr = imageio.pix2wld( info.getTransform(), nBlockBottomX, nBlockBottomY )
        
        # set this back to our copy of the info object
        info.setBlockBounds(blocktl,blockbr)

        # work out number of pixels of this block
        blockwidth = nBlockBottomX - xcoord
        blockheight = nBlockBottomY - ycoord
        
        # set this back to our copy of the info object
        info.setBlockSize(blockwidth,blockheight)
        
        # start creating our tuple. Start with an empty list
        # and append the blocks.
        blockList = []
        
        try:
            i = 0
            
            # read all the files using our iterable InputCollection
            for (image,ds,pixgrid,nullValList,datatype) in self.inputs:
            
                # get the pixel coords for this block for this file
                tl = imageio.wld2pix(pixgrid.makeGeoTransform(),blocktl.x,blocktl.y)
            
                # just read in the dataset (will return how many layers it has)
                # will just use the datatype of the image
                block = self.readBlockWithMargin(ds,int(round(tl.x)),int(round(tl.y)),blockwidth,blockheight,
                             datatype, margin=self.overlap, nullValList=nullValList,
                             layerselection=self.layerselectionList[i])

                # add this block to our list
                blockList.append(block)
            
                # set the relationship between numpy array
                # and dataset in case the user needs the dataset object
                # and/or the original filename
                info.setBlockDataset(block, ds, image)
                
                i += 1
                
        finally:
            # if there is any exception thrown here, make
            # sure temporary resampled files are deleted.
            # doesn't seem the destructor is called in this case.
            self.inputs.cleanup()
        
        
        if isinstance(self.imageContainer,dict):
            # we need to use the original keys passed in
            # to the constructor and return a dictionary
            blockDict = {}
            i = 0
            for name in self.imageContainer.keys():
                filename = self.imageContainer[name]
                if isinstance(filename, list):
                    listLen = len(filename)
                    blockDict[name] = []
                    for j in range(listLen):
                        blockDict[name].append(blockList[i])
                        i += 1
                elif isinstance(filename, basestring):
                    blockDict[name] = blockList[i]
                    i += 1
                                    
            # blockContainer is a dictionary
            blockContainer = blockDict
         
        elif isinstance(self.imageContainer,basestring):
            # blockContainer is just a single block
            blockContainer = blockList[0]

        else:   
            # blockContainer is a tuple
            blockContainer = tuple(blockList)
            
        # return a tuple with the info object and
        # our blockContainer
        return (info, blockContainer)
        
        
    @staticmethod
    def readBlockWithMargin(ds, xoff, yoff, xsize, ysize, datatype, margin=0, nullValList=None,
            layerselection=None):
        """
        A 'drop-in' look-alike for the ReadAsArray function in GDAL,
        but with the option of specifying a margin width, such that
        the block actually read and returned will be larger by that many pixels. 
        The returned array will ALWAYS contain these extra rows/cols, and 
        if they do not exist in the file (e.g. because the margin would push off 
        the edge of the file) then they will be filled with the given nullVal. 
        Otherwise they will be read from the file along with the rest of the block. 
        
        Variables within this function which have _margin as suffix are intended to 
        designate variables which include the margin, as opposed to those without. 
        
        This routine will cope with any specified region, even if it is entirely outside
        the given raster. The returned block would, in that case, be filled
        entirely with the null value. 
        
        """
        if layerselection is None:
            layerselection = [i+1 for i in range(ds.RasterCount)]
        nLayers = len(layerselection)
        
        # Create the final array, with margin, but filled with the null value. 
        xSize_margin = xsize + 2 * margin
        ySize_margin = ysize + 2 * margin
        outBlockShape = (nLayers, ySize_margin, xSize_margin)
        
        # Create the empty output array, filled with the appropriate null value. 
        block_margin = numpy.zeros(outBlockShape, dtype=datatype)
        if nullValList is not None and len(nullValList) > 0:
            # We really need something as a fill value, so if any of the 
            # null values in the list is None, then replace it with 0. 
            fillValList = [nullVal for nullVal in nullValList]
            for i in range(len(fillValList)):
                if fillValList[i] is None:
                    fillValList[i] = 0
            # Now use the appropriate null value for each layer as the 
            # initial value in the output array for the block. 
            if len(outBlockShape) == 2:
                block_margin.fill(fillValList[0])
            else:
                for i in range(nLayers):
                    block_margin[i].fill(fillValList[layerselection[i]-1])
#                for (i, fillVal) in enumerate(fillValList):
#                    block_margin[i].fill(fillVal)
        
        
        # Calculate the bounds of the block which we will actually read from the file,
        # based on what we have been asked for, what margin size, and how close we
        # are to the edge of the file. 
        
        # The bounds of the whole image in the file
        imgLeftBound = 0
        imgTopBound = 0
        imgRightBound = ds.RasterXSize
        imgBottomBound = ds.RasterYSize
        
        # The region we will, in principle, read from the file. Note that xSize_margin 
        # and ySize_margin are already calculated above
        xoff_margin = xoff - margin
        yoff_margin = yoff - margin
        
        # Restrict this to what is available in the file
        xoff_margin_file = max(xoff_margin, imgLeftBound)
        xoff_margin_file = min(xoff_margin_file, imgRightBound)
        xright_margin_file = xoff_margin + xSize_margin
        xright_margin_file = min(xright_margin_file, imgRightBound)
        xSize_margin_file = xright_margin_file - xoff_margin_file

        yoff_margin_file = max(yoff_margin, imgTopBound)
        yoff_margin_file = min(yoff_margin_file, imgBottomBound)
        ySize_margin_file = min(ySize_margin, imgBottomBound - yoff_margin_file)
        ybottom_margin_file = yoff_margin + ySize_margin
        ybottom_margin_file = min(ybottom_margin_file, imgBottomBound)
        ySize_margin_file = ybottom_margin_file - yoff_margin_file
        
        # How many pixels on each edge of the block we end up NOT reading from 
        # the file, and thus have to leave as null in the array
        notRead_left = xoff_margin_file - xoff_margin
        notRead_right = xSize_margin - (notRead_left + xSize_margin_file)
        notRead_top = yoff_margin_file - yoff_margin
        notRead_bottom = ySize_margin - (notRead_top + ySize_margin_file)
        
        # The upper bounds on the slices specified to receive the data
        slice_right = xSize_margin - notRead_right
        slice_bottom = ySize_margin - notRead_bottom
        
        if xSize_margin_file > 0 and ySize_margin_file > 0:
            # Now read in the part of the array which we can actually read from the file.
            # Read each layer separately, to honour the layerselection
            
            # The part of the final array we are filling
            imageSlice = (slice(notRead_top, slice_bottom), slice(notRead_left, slice_right))
            
            for i in range(nLayers):
                band = ds.GetRasterBand(layerselection[i])
                block_margin[i][imageSlice] = band.ReadAsArray(xoff_margin_file, yoff_margin_file, 
                    xSize_margin_file, ySize_margin_file)

        return block_margin
        
    def close(self):
        """
        Closes all open datasets
        """
        self.inputs.close()

