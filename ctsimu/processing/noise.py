import os    # File and path handling
import numpy
import copy
from ..image import *
from ..helpers import *

from .pipeline import Pipeline
from .step import Step

class Step_Noise(Step):
    """ Add noise to image according to SNR characteristics.

        The noise characteristics must be specified by two lists:
            1. grey values
            2. assigned SNRs

        The characteristics must be sorted by grey values in ascending order.
        Linear interpolation will take place for missing grey values.
    """

    def __init__(self, sigma=None, greyValues=None, SNR=None):
        Step.__init__(self, "Noise")
        self.sigma = None
        self.greyValues = None
        self.SNR = None

        self.setSigma(sigma)
        self.setNoiseCharacteristics(greyValues, SNR)

    def setSigma(self, sigma=None):
        if sigma is None:
            self.sigma = 0
        else:
            self.sigma = sigma

    def setNoiseCharacteristics(self, greyValues=None, SNR=None):
        self.greyValues = greyValues
        self.SNR = SNR

    def prepare(self):
        """ Nothing to prepare for this module. """
        if isinstance(self.pipe, Pipeline):
            self.prepared = True
            return

        self.prepared = False
        raise Exception("Step must be part of a processing pipeline before it can prepare.")

    def run(self, image):
        """ Transform given image. """
        sigma = copy.deepcopy(image)

        if self.greyValues is None or self.SNR is None:
            # Assign constant sigma to each pixel:
            sigma.erase(self.sigma)
        else:
            # Map grey values to SNR:
            sigma.map(gv_from=self.greyValues, gv_to=self.SNR, bins=1000)

            # Calculate sigma from sigma = I / SNR where SNR>0:
            sigma.px = numpy.where(sigma.px > 0, image.px / sigma.px, 0)
        
        image.noise(sigma.px)

        return image