"""
Progress bars using CUI interface. 
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
from osgeo.gdal import TermProgress_nocb

class CUIProgressBar(object):
    """
    Simple progress as a percentage pronted to the terminal
    """
    def __init__(self):
        self.totalsteps = 100

    def setTotalSteps(self,steps):
        self.totalsteps = steps

    def setProgress(self,progress):
        if sys.stdout.isatty(): # don't write if a log file etc
            progress = int(float(progress) / self.totalsteps * 100)
            sys.stdout.write('%d%%\r' % progress)

    def reset(self):
        sys.stdout.write('\n')

    def setLabelText(self,text):
        sys.stdout.write('\n%s\n' % text)

    def wasCancelled(self):
        return False

    def displayException(self,trace):
        sys.stdout.write(trace)

    def displayWarning(self,text):
        sys.stdout.write("Warning: %s\n" % text)

    def displayError(self,text):
        sys.stdout.write("Error: %s\n" % text)

    def displayInfo(self,text):
        sys.stdout.write("Info: %s\n" % text)

class GDALProgressBar(object):
    """
    Uses GDAL's TermProgress to print a progress bar to the terminal
    """
    def __init__(self):
        self.totalsteps = 100

    def setTotalSteps(self,steps):
        self.totalsteps = steps

    def setProgress(self,progress):
        if sys.stdout.isatty(): # don't write if a log file etc
            TermProgress_nocb(float(progress) / self.totalsteps)

    def reset(self):
        if sys.stdout.isatty(): # don't write if a log file etc
            TermProgress_nocb(100.0)
        sys.stdout.write('\n')

    def setLabelText(self,text):
        sys.stdout.write('\n%s\n' % text)

    def wasCancelled(self):
        return False

    def displayException(self,trace):
        sys.stdout.write(trace)

    def displayWarning(self,text):
        sys.stdout.write("Warning: %s\n" % text)

    def displayError(self,text):
        sys.stdout.write("Error: %s\n" % text)

    def displayInfo(self,text):
        sys.stdout.write("Info: %s\n" % text)

class SilentProgress(object):
    """
    A progress object which is completely silent. 
    """
    def __init__(self):
        pass
    def setTotalSteps(self, steps):
        pass
    def setProgress(self,progress):
        pass
    def reset(self):
        pass
    def setLabelText(self,text):
        pass
    def wasCancelled(self):
        return False
    def displayException(self,trace):
        pass
    def displayWarning(self,text):
        pass
    def displayError(self,text):
        pass
    def displayInfo(self,text):
        pass

