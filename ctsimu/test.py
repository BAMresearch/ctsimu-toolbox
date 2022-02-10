# -*- coding: UTF-8 -*-

import os    # File and path handling
import numpy
import copy  # for deepcopy
import math

from .image import ImageFile, Image, ImageROI, ImageStack
from .geometry import Geometry
from .processing.pipeline import Pipeline
from .processing.step import Step
from .helpers import *

def touchDirectory(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)

class generalTest(Step):
    """ General class for test scenario evaluations: get image(s), run and store evaluation. """

    def __init__(self, testName="General Test", name=None, nExpectedRuns=1, resultFileDirectory=".", rawOutput=False):
        Step.__init__(self, testName)
       
        self.testName = testName
        self.subtests = []

        self.prepared = False 
        self.currentRun = 0
        self.nExpectedRuns = None  # usually, number of projections to evaluate
        self.resultFileDirectory = None
        self.name = None
        self.rawOutput = None

        self.setName(name)
        self.setExpectedRuns(nExpectedRuns)
        self.setResultFileDirectory(resultFileDirectory)
        self.setRawOutput(rawOutput)

        self.reset()

    def reset(self):
        self.currentRun = 0
        self.prepared = False

    def addSubtest(self, subt):
        self.subtests.append(subt)

    def setName(self, name=None):
        """ Set an individual name for the (sub) test. """
        if name != None:
            self.name = name
        else:
            self.name = self.testName

    def setExpectedRuns(self, n=1):
        self.nExpectedRuns = n

    def setResultFileDirectory(self, resultFileDirectory="."):
        """ Set the location where test results should be saved. """
        self.resultFileDirectory = resultFileDirectory
        touchDirectory(self.resultFileDirectory)

    def setRawOutput(self, rawOutput=False):
        """ Save intermediate projections as RAW instead of TIFF? """
        self.rawOutput = rawOutput

    def plotResults(self):
        """ Plot results of evaluation. """
        # Should be called by step's followUp() function, if needed.
        pass
