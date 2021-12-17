import os    # File and path handling
import numpy
import copy
from ..image import *
from ..general import *

from .pipeline import ProcessingPipeline
from .processing_step import ProcessingStep

class ProcessingStep_Smoothing(ProcessingStep):
    """ Smooth an image. """

    def __init__(self, mode="gaussian", sigma=0):
        ProcessingStep.__init__(self, "Smoothing")
        self._sigma = None
        self._mode  = None

        self.setSigma(sigma)
        self.setMode(mode)

    def setSigma(self, sigma=None):
        if sigma is None:
            self._sigma = 0
        else:
            self._sigma = sigma

    def setMode(self, mode):
        self._mode = mode

    def prepare(self):
        """ Nothing to prepare for this module. """
        if isinstance(self._pipe, ProcessingPipeline):
            self._prepared = True
            return

        self._prepared = False
        raise Exception("Step must be part of a processing pipeline before it can prepare.")

    def run(self, image):
        """ Transform given image. """
        if self._mode is not None:
            if self._mode == "gaussian":
                image.smooth_gaussian(sigma=self._sigma)

        return image