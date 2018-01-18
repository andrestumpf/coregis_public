#!/usr/bin/env python
"""
Testing deleteIfExisting()

"""

from __future__ import print_function, division

from osgeo import gdal

from rios import applier

def main():
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    
    infiles.img = 'crap.img'
    #outfiles.outimg = 'crapz.img'
    outfiles.outimg = "/apollo/imagery/rsc/landsat/landsat57tm/wrs2/090_079/2015/201512/l7tmpa_p090r079_20151228_da0m6.img"
    
    applier.apply(doit, infiles, outfiles)

def doit(info, inputs, outputs):
    outputs.outimg = inputs.img

if __name__ == "__main__":
    main()
