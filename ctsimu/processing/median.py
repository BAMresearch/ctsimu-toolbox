import os    # File and path handling
import numpy
from ..image import *
from ..helpers import *

from .pipeline import ProcessingPipeline
from .step import ProcessingStep

class ProcessingStep_Median(ProcessingStep):
    """ Binning operation for the processing pipeline. """

    def __init__(self, medianSize=3):
        ProcessingStep.__init__(self, "Median")
        self.size = medianSize

    def setSize(self, medianSize=3):
        self.setPrepared(False)
        self.size = medianSize

    def getSize(self):
        return self.size

    def prepare(self):
        """ Nothing to prepare for the binning module. """
        if isinstance(self.pipe, ProcessingPipeline):
            self.prepared = True
            return

        self.prepared = False
        raise Exception("Step must be part of a processing pipeline before it can prepare.")

    def run(self, image):
        """ Bin given image. """
        self.prepare()

        image.applyMedian(kernelSize=self.size)

        return image