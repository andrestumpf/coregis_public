#!/usr/bin/env python
"""
Test the input layer selection feature. 

Generates a multi-band image, then applies a RIOS function which only reads 
two of the layers. Checks that it gets the correct answer for adding them together. 

"""
import os
import numpy

from rios.riostests import riostestutils
from rios import applier

TESTNAME = "TESTLAYERSELECT"

def run():
    """
    Run the test
    """
    riostestutils.reportStart(TESTNAME)

    # Create a multi-band file with some data in it. 
    tstfile = 'multilayer.img'
    numBands = 5
    ds = riostestutils.createTestFile(tstfile, numBands=numBands)
    
    onelayerArr = riostestutils.genRampArray()
    lyrList = []
    for i in range(numBands):
        lyr = (onelayerArr + 1)
        band = ds.GetRasterBand(i+1)
        band.WriteArray(lyr)
        lyrList.append(lyr)
    del ds
        
    stack = numpy.array(lyrList)

    # Sum of all pixels in bands 2 & 4
    layerList = [2, 4]
    correctSum = sum([(stack[i-1].astype(numpy.float64)).sum() for i in layerList])
    
    # Now do it using RIOS
    infiles = applier.FilenameAssociations()
    outfiles = applier.FilenameAssociations()
    otherargs = applier.OtherInputs()
    controls = applier.ApplierControls()
    
    infiles.img = tstfile
    controls.selectInputImageLayers(layerList)
    otherargs.total = 0
    # We will use this to check the number of layers being read
    otherargs.numLayers = len(layerList)
    otherargs.numLayersIsOK = True
    
    applier.apply(doSum, infiles, outfiles, otherargs, controls=controls)
    
    if correctSum != otherargs.total:
        riostestutils.report(TESTNAME, "Totals do not match: %s != %s"%(correctSum, otherargs.total))
        ok = False
    else:
        riostestutils.report(TESTNAME, "Passed")
        ok = True
    
    os.remove(tstfile)
    
    return ok


def doSum(info, inputs, outputs, otherargs):
    """
    Should be given only the two layers we want, therefore can calculate the sum 
    by adding all layers
    """
    # First check that the number of layers is right
    if len(inputs.img) != otherargs.numLayers:
        otherargs.numLayersIsOK = False
    
    # Now accumulate the total
    otherargs.total += inputs.img.astype(numpy.float64).sum()


if __name__ == "__main__":
    run()
