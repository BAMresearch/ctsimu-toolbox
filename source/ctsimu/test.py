# -*- coding: UTF-8 -*-

import os    # File and path handling
import numpy
import copy  # for deepcopy
import math

from .image import ImageFile, Image, ImageROI, ImageStack
from .geometry import Geometry
from .processing import ProcessingPipeline, ProcessingStep
from .general import *

def touchDirectory(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)

class generalTest(ProcessingStep):
    """ General class for test scenario evaluations: get image(s), run and store evaluation. """

    def __init__(self, testName="General Test", name=None, nExpectedRuns=1, resultFileDirectory=".", rawOutput=False):
        ProcessingStep.__init__(self, testName)
       
        self._testName = testName
        self._subtests = []

        self._prepared = False 
        self._currentRun = 0
        self._nExpectedRuns = None  # usually, number of projections to evaluate
        self._resultFileDirectory = None
        self._name = None
        self._rawOutput = None

        self.setName(name)
        self.setExpectedRuns(nExpectedRuns)
        self.setResultFileDirectory(resultFileDirectory)
        self.setRawOutput(rawOutput)

        self.reset()

    def reset(self):
        self._currentRun = 0
        self._prepared = False

    def addSubtest(self, subt):
        self._subtests.append(subt)

    def setName(self, name=None):
        """ Set an individual name for the (sub) test. """
        if name != None:
            self._name = name
        else:
            self._name = self._testName

    def setExpectedRuns(self, n=1):
        self._nExpectedRuns = n

    def setResultFileDirectory(self, resultFileDirectory="."):
        """ Set the location where test results should be saved. """
        self._resultFileDirectory = resultFileDirectory
        touchDirectory(self._resultFileDirectory)

    def setRawOutput(self, rawOutput=False):
        """ Save intermediate projections as RAW instead of TIFF? """
        self._rawOutput = rawOutput

    def plotResults(self):
        """ Plot results of evaluation. """
        # Should be called by step's followUp() function, if needed.
        pass
