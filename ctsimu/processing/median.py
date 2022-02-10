import os    # File and path handling
import numpy
from ..image import *
from ..helpers import *

from .pipeline import Pipeline
from .step import Step

class Step_Median(Step):
    """ Binning operation for the processing pipeline. """

    def __init__(self, medianSize=3):
        Step.__init__(self, "Median")
        self.size = medianSize

    def setSize(self, medianSize=3):
        self.setPrepared(False)
        self.size = medianSize

    def getSize(self):
        return self.size

    def prepare(self):
        """ Nothing to prepare for the binning module. """
        if isinstance(self.pipe, Pipeline):
            self.prepared = True
            return

        self.prepared = False
        raise Exception("Step must be part of a processing pipeline before it can prepare.")

    def run(self, image):
        """ Bin given image. """
        self.prepare()

        image.applyMedian(kernelSize=self.size)

        return image