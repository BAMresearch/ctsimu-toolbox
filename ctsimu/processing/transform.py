import os    # File and path handling
import numpy
from ..image import *
from ..helpers import *

from .pipeline import Pipeline
from .step import Step

class Step_Transform(Step):
    """ Rotations and mirroring (flip). """

    def __init__(self, rotate=0, flipX=False, flipY=False):
        Step.__init__(self, "Transform")
        self.rotate = 0
        self.flipX  = False
        self.flipY  = False

        self.setRotation(rotate)
        self.setFlip(flipX, flipY)

    def setRotation(self, rotate=0):
        if rotate is None:
            self.rotate = 0
        else:
            if rotate in [0, 90, 180, 270]:
                self.rotate = rotate
            else:
                raise Exception("Image rotation must be one of the following integers: 0, 90, 180, 270.")

    def setFlip(self, flipX=False, flipY=False):
        self.flipX = flipX
        self.flipY = flipY

    def prepare(self):
        """ Nothing to prepare for the transform module. """
        if isinstance(self.pipe, Pipeline):
            self.prepared = True
            return

        self.prepared = False
        raise Exception("Step must be part of a processing pipeline before it can prepare.")

    def run(self, image):
        """ Transform given image. """
        self.prepare()

        if self.rotate != 0:
            image.rotate("{}".format(self.rotate))

        image.flip(horizontal=self.flipX, vertical=self.flipY)

        return image