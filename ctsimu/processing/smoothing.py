import os    # File and path handling
import numpy
import copy
from ..image import *
from ..helpers import *

from .pipeline import Pipeline
from .step import Step

class Step_Smoothing(Step):
    """ Smooth an image. """

    def __init__(self, mode="gaussian", sigma=0):
        Step.__init__(self, "Smoothing")
        self.sigma = None
        self.mode  = None

        self.setSigma(sigma)
        self.setMode(mode)

    def setSigma(self, sigma=None):
        if sigma is None:
            self.sigma = 0
        else:
            self.sigma = sigma

    def setMode(self, mode):
        self.mode = mode

    def prepare(self):
        """ Nothing to prepare for this module. """
        if isinstance(self.pipe, Pipeline):
            self.prepared = True
            return

        self.prepared = False
        raise Exception("Step must be part of a processing pipeline before it can prepare.")

    def run(self, image):
        """ Transform given image. """
        if self.mode is not None:
            if self.mode == "gaussian":
                image.smooth_gaussian(sigma=self.sigma)

        return image