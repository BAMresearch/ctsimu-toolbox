import os    # File and path handling
import numpy
from ..image import *
from ..general import *

from .pipeline import ProcessingPipeline
from .processing_step import ProcessingStep

class ProcessingStep_Transform(ProcessingStep):
    """ Rotations and mirroring (flip). """

    def __init__(self, rotate=0, flipX=False, flipY=False):
        ProcessingStep.__init__(self, "Transform")
        self._rotate = 0
        self._flipX  = False
        self._flipY  = False

        self.setRotation(rotate)
        self.setFlip(flipX, flipY)

    def setRotation(self, rotate=0):
        if rotate is None:
            self._rotate = 0
        else:
            if rotate in [0, 90, 180, 270]:
                self._rotate = rotate
            else:
                raise Exception("Image rotation must be one of the following integers: 0, 90, 180, 270.")

    def setFlip(self, flipX=False, flipY=False):
        self._flipX = flipX
        self._flipY = flipY

    def prepare(self):
        """ Nothing to prepare for the transform module. """
        if isinstance(self._pipe, ProcessingPipeline):
            self._prepared = True
            return

        self._prepared = False
        raise Exception("Step must be part of a processing pipeline before it can prepare.")

    def run(self, image):
        """ Transform given image. """
        self.prepare()

        if self._rotate != 0:
            image.rotate("{}".format(self._rotate))

        image.flip(horizontal=self._flipX, vertical=self._flipY)

        return image