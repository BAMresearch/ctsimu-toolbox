# -*- coding: UTF-8 -*-
"""
PDF report generation.
"""

import os    # File and path handling
from reportlab.pdfgen import canvas 
from reportlab.pdfbase.ttfonts import TTFont 
from reportlab.pdfbase import pdfmetrics 
from reportlab.lib import colors 
from PIL import Image
from datetime import datetime
from .version import *

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm


class Report(canvas.Canvas):

    def __init__(self,
                 filename, 
                 pagesize: tuple[float, float] | None = A4):
        
        canvas.Canvas.__init__(
            self,
            filename, 
            pagesize, 
            bottomup=1)
        
        self.pagesize = pagesize
        self.margin = [25.0 * mm, 20.0 * mm, 27.0 * mm, 25.0 * mm]
        self.box = [self.margin[0], self.pagesize[1] - self.margin[2],
                    self.pagesize[0] - self.margin[1], self.margin[3]]
        self.boxWidth = self.box[2] - self.box[0]
        self.boxXCenter = self.box[0] + .5 * self.boxWidth

    def getBox(self):
        return self.box

    def getBoxXCenter(self):
        return self.boxXCenter

    def drawHeader(self,
                   left = '',
                   center = '',
                   right = '',
                   pad = 6.0):
        yLine = self.box[1] + pad * mm
        yText = yLine + 2 * mm
        self.drawString(self.box[0], yText, left)
        self.drawCentredString(self.boxXCenter, yText, center)
        self.drawRightString(self.box[2], yText, right)
        self.line(self.box[0], yLine, self.box[2], yLine)

    def drawFooter(self,
                   left = '',
                   center = '',
                   right = '',
                   pad = 6.0):
        yLine = self.box[3] - pad * mm
        yText = yLine - 6 * mm
        self.drawString(self.box[0], yText, left)
        self.drawCentredString(self.boxXCenter, yText, center)
        self.drawRightString(self.box[2], yText, right)
        self.line(self.box[0], yLine, self.box[2], yLine)
