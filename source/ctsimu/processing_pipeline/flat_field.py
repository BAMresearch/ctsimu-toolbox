import os    # File and path handling
import numpy
from ..image import *
from ..general import *
from ..geometry import Geometry   # for analytical flat field

from .pipeline import ProcessingPipeline
from .processing_step import ProcessingStep

class ProcessingStep_FlatFieldCorrection(ProcessingStep):
    """ Flat field correction for the processing pipeline. """

    def __init__(self, flatFileStack=None, darkFileStack=None):
        ProcessingStep.__init__(self, "Flatfield Correction")

        # Parameters:
        self._flatFileStack        = None
        self._darkFileStack        = None
        self._flatFieldRescaleFactor = 1
        self._offsetAfterRescale     = 0
        self._doOffsetCorrection   = False
        self._doGainCorrection     = False

        self._analyticalCorrection = False
        self._jsonScene = None
        self._analyticalOffset = 0

        # Stored image data:
        self._darkImage = None
        self._gainImage = None

        self.setFlatFileStack(flatFileStack)
        self.setDarkFileStack(darkFileStack)

    def setFlatFileStack(self, flatFileStack):
        self._flatFileStack = createImageStack(flatFileStack)
        self.setPrepared(False)
        if(flatFileStack != None):
            self.setGainCorrection(True)
        else:
            self.setGainCorrection(False)

    def setDarkFileStack(self, darkFileStack):
        self._darkFileStack = createImageStack(darkFileStack)
        self.setPrepared(False)
        if(darkFileStack != None):
            self.setOffsetCorrection(True)
        else:
            self.setOffsetCorrection(False)

    def setGainCorrection(self, doGainCorrection):
        self._doGainCorrection = doGainCorrection
        self.setPrepared(False)

    def setOffsetCorrection(self, doOffsetCorrection):
        self._doOffsetCorrection = doOffsetCorrection
        self.setPrepared(False)

    def setFlatFieldRescaleFactor(self, factor, offsetAfterRescale=0):
        self._flatFieldRescaleFactor = factor
        self._offsetAfterRescale     = offsetAfterRescale
        self.setPrepared(False)

    def setAnalyticalCorrection(self, doAnalytical, jsonScenarioFile=None, offsetValue=0, flatFieldRescaleFactor=1):
        """ Configure the analytical flat field correction. """
        self._analyticalCorrection = doAnalytical
        self._jsonScene = jsonScenarioFile
        self._analyticalOffset = offsetValue
        self.setFlatFieldRescaleFactor(flatFieldRescaleFactor)
        self.setPrepared(False)
        if doAnalytical:
            self.setGainCorrection(True)
            if self._analyticalOffset != 0:
                self.setOffsetCorrection(True)
            else:
                self.setOffsetCorrection(False)

    def prepare(self):
        """ Prepare by reading/generating dark and flat field images. """
        if not isinstance(self._pipe, ProcessingPipeline):
            self._prepared = False
            raise Exception("Step must be part of a processing pipeline before it can prepare. Current: {}".format(self._pipe))

        if not self._prepared:
            outputFiles = None
            if isinstance(self._pipe._outputFileStack, ImageStack):
                outputFiles = self._pipe._outputFileStack._files

            # Read dark field image if offset correction is active:
            if(self._doOffsetCorrection):
                if(not self._analyticalCorrection):
                    self._darkFileStack.buildStack()
                    if self._darkFileStack.nSlices() > 0:
                        self._darkImage = self._darkFileStack.getMeanImage(outputFiles)
                    else:
                        log("No offset file(s) found that match the given name pattern.")

            # Read flat field image if gain correction is active:
            if(self._doGainCorrection):
                if(not self._analyticalCorrection):
                    self._flatFileStack.buildStack()
                    if self._flatFileStack.nSlices() > 0:
                        self._gainImage = self._flatFileStack.getMeanImage(outputFiles)
                    else:
                        log("No flat field file(s) found that match the given name pattern.")
                   
                    # Do an offset correction of the flat field image:
                    if(self._doOffsetCorrection):
                        self._gainImage.applyDark(self._darkImage)

                elif(self._analyticalCorrection): # Analytical flat field
                    if(self._jsonScene != None):
                        ctsimuGeometry = Geometry(self._jsonScene)
                        self._gainImage = ctsimuGeometry.createDetectorFlatField_analytical()
                        if self._flatFieldRescaleFactor != 1:
                            print("WARNING: For analytical correction, a rescale factor of 1 is recommended. Your current choice: {}".format(self._flatFieldRescaleFactor))
                    else:
                        raise Exception("A JSON scenario description must be specified when using analytical flat field correction.")

            self._prepared = True

    def run(self, image):
        """ Run offset and flat field correction. """
        self.prepare()
        if self._doOffsetCorrection:
            if self._analyticalCorrection:
                if self._analyticalOffset != 0:
                    image.subtract(self._analyticalOffset)
            else:
                image.applyDark(self._darkImage)

        if self._doGainCorrection:
            image.applyFlatfield(self._gainImage, self._flatFieldRescaleFactor)

        if self._offsetAfterRescale != 0:
            image.add(self._offsetAfterRescale)

        return image