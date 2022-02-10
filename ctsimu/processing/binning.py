import os    # File and path handling
import numpy
from ..image import *
from ..helpers import *

from .pipeline import Pipeline
from .step import Step

class Step_Binning(Step):
    """ Binning operation for the processing pipeline. """

    def __init__(self, binSizeX=1, binSizeY=1, binningOperation="mean"):
        Step.__init__(self, "Binning")
        self.binSizeX = 1
        self.binSizeY = 1
        self.binningOperation = "mean"

        self.validBinningOperations = ["mean", "min", "max", "sum", None]

        self.setBinning(binSizeX, binSizeY, binningOperation)

    def setBinning(self, x, y, operation='mean'):
        self.setPrepared(False)

        if x is None:
            x = 1

        if y is None:
            y = 1

        if (x >= 0) and (y >= 0):
            if x == 0:
                x = 1
            if y == 0:
                y = 1

            self.binSizeX = int(x)
            self.binSizeY = int(y)
        else:
            raise Exception("The bin size must be >= 1.")

        if operation in self.validBinningOperations:
            self.binningOperation = operation
        else:
            raise Exception("'{}' is not a valid binning operation. Options are: {}".format(operation, self.validBinningOperations))

    def getBinSizeX(self):
        return self.binSizeX

    def getBinSizeY(self):
        return self.binSizeY

    def getBinningOperation(self):
        return self.binningOperation

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

        if (self.getBinSizeX() > 1) or (self.getBinSizeY() > 1):
            image.bin(self.getBinSizeX(), self.getBinSizeY(), self.getBinningOperation())

        return image