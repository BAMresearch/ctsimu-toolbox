import os    # File and path handling
import numpy
from ..image import *
from ..general import *

from .pipeline import ProcessingPipeline
from .processing_step import ProcessingStep

class ProcessingStep_Binning(ProcessingStep):
    """ Binning operation for the processing pipeline. """

    def __init__(self, binSizeX=1, binSizeY=1, binningOperation="mean"):
        ProcessingStep.__init__(self, "Binning")
        self._binSizeX = 1
        self._binSizeY = 1
        self._binningOperation = "mean"

        self._validBinningOperations = ["mean", "min", "max", "sum", None]

        self.setBinning(binSizeX, binSizeY, binningOperation)

    def setBinning(self, x, y, operation='mean'):
        self.setPrepared(False)

        if x == None:
            x = 1

        if y == None:
            y = 1

        if (x >= 0) and (y >= 0):
            if x == 0:
                x = 1
            if y == 0:
                y = 1

            self._binSizeX = int(x)
            self._binSizeY = int(y)
        else:
            raise Exception("The bin size must be >= 1.")

        if operation in self._validBinningOperations:
            self._binningOperation = operation
        else:
            raise Exception("'{}' is not a valid binning operation. Options are: {}".format(operation, self._validBinningOperations))

    def getBinSizeX(self):
        return self._binSizeX

    def getBinSizeY(self):
        return self._binSizeY

    def getBinningOperation(self):
        return self._binningOperation

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

        if (self.getBinSizeX() > 1) or (self.getBinSizeY() > 1):
            image.bin(self.getBinSizeX(), self.getBinSizeY(), self.getBinningOperation())

        return image