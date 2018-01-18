#!/usr/bin/env python
"""
Use rios.fileinfo to print the statistics for the given image(s). 

"""
from __future__ import print_function

import argparse

from rios import fileinfo

def getCmdargs():
    """
    Get commandline arguments
    """
    p = argparse.ArgumentParser()
    p.add_argument("imgfile", nargs='*', help="Name of input image file")
    p.add_argument("--printfilename", default=False, action="store_true",
        help=("Print each filename at the start of each line of output "+
            "(only useful with multiple input files, to distinguish "+
            "which line belongs with which file)"))
    cmdargs = p.parse_args()
    return cmdargs
    

def main():
    """
    Main routine
    """
    cmdargs = getCmdargs()
    
    for filename in cmdargs.imgfile:
        stats = fileinfo.ImageFileStats(filename)
        for layerStats in stats:
            outStr = str(layerStats)
            if cmdargs.printfilename:
                outline = "File:%s, %s" % (filename, outStr)
            else:
                outline = outStr
            
            print(outline)


if __name__ == "__main__":
    main()
