import os    # File and path handling
import numpy
from ..image import *
from ..general import *

from .pipeline import ProcessingPipeline
from .processing_step import ProcessingStep

class ProcessingStep_Median(ProcessingStep):
    """ Binning operation for the processing pipeline. """

    def __init__(self, medianSize=3):
        ProcessingStep.__init__(self, "Median")
        self._size = medianSize

    def setSize(self, medianSize=3):
        self.setPrepared(False)
        self._size = medianSize

    def getSize(self):
        return self._size

    def prepare(self):
        """ Nothing to prepare for the binning module. """
        if isinstance(self._pipe, ProcessingPipeline):
            self._prepared = True
            return

        self._prepared = False
        raise Exception("Step must be part of a processing pipeline before it can prepare.")

    def run(self, image):
        """ Bin given image. """
        self.prepare()

        image.applyMedian(kernelSize=self._size)

        return image