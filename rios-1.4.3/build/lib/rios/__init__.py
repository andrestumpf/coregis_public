"""
Raster Input Output Simplification

A Python package to simplify raster input/output
using GDAL, allowing a programmer to focus on the
processing of the data rather than the mechanics of
raster I/O. 

Rios is built on top of GDAL, and handles, among other things: 
    - Opening raster files and reading/writing raster data
    - Determining the area of intersection (or union) of 
      input rasters
    - Reprojecting inputs to be on the same projection and
      pixel alignment
    - Stepping through the rasters in small blocks, to avoid 
      large memory usage

Most common entry point is the apply function in the applier 
module. For more subtle and complex work the imagereader and 
imagewriter modules provide the ImageReader and ImageWriter 
classes, respectively. 

"""

from distutils.version import LooseVersion

RIOS_VERSION = '1.4.3'
RIOS_VERSION_OBJ = LooseVersion(RIOS_VERSION)
__version__ = RIOS_VERSION
