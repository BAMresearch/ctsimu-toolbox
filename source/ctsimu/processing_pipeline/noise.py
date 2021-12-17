import os    # File and path handling
import numpy
import copy
from ..image import *
from ..general import *

from .pipeline import ProcessingPipeline
from .processing_step import ProcessingStep

class ProcessingStep_Noise(ProcessingStep):
    """ Add noise to image according to SNR characteristics.

        The noise characteristics must be specified by two lists:
            1. grey values
            2. assigned SNRs

        The characteristics must be sorted by grey values in ascending order.
        Linear interpolation will take place for missing grey values.
    """

    def __init__(self, sigma=None, greyValues=None, SNR=None):
        ProcessingStep.__init__(self, "Noise")
        self._sigma = None
        self._greyValues = None
        self._SNR = None

        self.setSigma(sigma)
        self.setNoiseCharacteristics(greyValues, SNR)

    def setSigma(self, sigma=None):
        if sigma is None:
            self._sigma = 0
        else:
            self._sigma = sigma

    def setNoiseCharacteristics(self, greyValues=None, SNR=None):
        self._greyValues = greyValues
        self._SNR = SNR

    def prepare(self):
        """ Nothing to prepare for this module. """
        if isinstance(self._pipe, ProcessingPipeline):
            self._prepared = True
            return

        self._prepared = False
        raise Exception("Step must be part of a processing pipeline before it can prepare.")

    def run(self, image):
        """ Transform given image. """
        sigma = copy.deepcopy(image)

        if self._greyValues is None or self._SNR is None:
            # Assign constant sigma to each pixel:
            sigma.erase(self._sigma)
        else:
            # Map grey values to SNR:
            sigma.map(gv_from=self._greyValues, gv_to=self._SNR, bins=1000)

            # Calculate sigma from sigma = I / SNR where SNR>0:
            sigma._px = numpy.where(sigma._px > 0, image._px / sigma._px, 0)
        
        image.noise(sigma._px)

        return image