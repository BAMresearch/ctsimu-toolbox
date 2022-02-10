import os    # File and path handling
import numpy
from ..image import *
from ..helpers import *
from ..geometry import Geometry   # for analytical flat field

from .pipeline import Pipeline
from .step import Step

class Step_FlatFieldCorrection(Step):
    """ Flat field correction for the processing pipeline. """

    def __init__(self, flatFileStack=None, darkFileStack=None):
        Step.__init__(self, "Flatfield Correction")

        # Parameters:
        self.flatFileStack        = None
        self.darkFileStack        = None
        self.flatFieldRescaleFactor = 1
        self.offsetAfterRescale     = 0
        self.doOffsetCorrection   = False
        self.doGainCorrection     = False

        self.analyticalCorrection = False
        self.jsonScene = None
        self.analyticalOffset = 0

        # Stored image data:
        self.darkImage = None
        self.gainImage = None

        self.setFlatFileStack(flatFileStack)
        self.setDarkFileStack(darkFileStack)

    def setFlatFileStack(self, flatFileStack):
        self.flatFileStack = createImageStack(flatFileStack)
        self.setPrepared(False)
        if(flatFileStack != None):
            self.setGainCorrection(True)
        else:
            self.setGainCorrection(False)

    def setDarkFileStack(self, darkFileStack):
        self.darkFileStack = createImageStack(darkFileStack)
        self.setPrepared(False)
        if(darkFileStack != None):
            self.setOffsetCorrection(True)
        else:
            self.setOffsetCorrection(False)

    def setGainCorrection(self, doGainCorrection):
        self.doGainCorrection = doGainCorrection
        self.setPrepared(False)

    def setOffsetCorrection(self, doOffsetCorrection):
        self.doOffsetCorrection = doOffsetCorrection
        self.setPrepared(False)

    def setFlatFieldRescaleFactor(self, factor, offsetAfterRescale=0):
        self.flatFieldRescaleFactor = factor
        self.offsetAfterRescale     = offsetAfterRescale
        self.setPrepared(False)

    def setAnalyticalCorrection(self, doAnalytical, jsonScenarioFile=None, offsetValue=0, flatFieldRescaleFactor=1):
        """ Configure the analytical flat field correction. """
        self.analyticalCorrection = doAnalytical
        self.jsonScene = jsonScenarioFile
        self.analyticalOffset = offsetValue
        self.setFlatFieldRescaleFactor(flatFieldRescaleFactor)
        self.setPrepared(False)
        if doAnalytical:
            self.setGainCorrection(True)
            if self.analyticalOffset != 0:
                self.setOffsetCorrection(True)
            else:
                self.setOffsetCorrection(False)

    def prepare(self):
        """ Prepare by reading/generating dark and flat field images. """
        if not isinstance(self.pipe, Pipeline):
            self.prepared = False
            raise Exception("Step must be part of a processing pipeline before it can prepare. Current: {}".format(self.pipe))

        if not self.prepared:
            outputFiles = None
            if isinstance(self.pipe.outputFileStack, ImageStack):
                outputFiles = self.pipe.outputFileStack.files

            # Read dark field image if offset correction is active:
            if(self.doOffsetCorrection):
                if(not self.analyticalCorrection):
                    self.darkFileStack.buildStack()
                    if self.darkFileStack.nSlices > 0:
                        self.darkImage = self.darkFileStack.getMeanImage(outputFiles)
                    else:
                        log("No offset file(s) found that match the given name pattern.")

            # Read flat field image if gain correction is active:
            if(self.doGainCorrection):
                if(not self.analyticalCorrection):
                    self.flatFileStack.buildStack()
                    if self.flatFileStack.nSlices > 0:
                        self.gainImage = self.flatFileStack.getMeanImage(outputFiles)
                    else:
                        log("No flat field file(s) found that match the given name pattern.")
                   
                    # Do an offset correction of the flat field image:
                    if(self.doOffsetCorrection):
                        self.gainImage.applyDark(self.darkImage)

                elif(self.analyticalCorrection): # Analytical flat field
                    if(self.jsonScene != None):
                        ctsimuGeometry = Geometry(self.jsonScene)
                        self.gainImage = ctsimuGeometry.createDetectorFlatField_analytical()
                        if self.flatFieldRescaleFactor != 1:
                            print("WARNING: For analytical correction, a rescale factor of 1 is recommended. Your current choice: {}".format(self.flatFieldRescaleFactor))
                    else:
                        raise Exception("A JSON scenario description must be specified when using analytical flat field correction.")

            self.prepared = True

    def run(self, image):
        """ Run offset and flat field correction. """
        self.prepare()
        if self.doOffsetCorrection:
            if self.analyticalCorrection:
                if self.analyticalOffset != 0:
                    image.subtract(self.analyticalOffset)
            else:
                image.applyDark(self.darkImage)

        if self.doGainCorrection:
            image.applyFlatfield(self.gainImage, self.flatFieldRescaleFactor)

        if self.offsetAfterRescale != 0:
            image.add(self.offsetAfterRescale)

        return image