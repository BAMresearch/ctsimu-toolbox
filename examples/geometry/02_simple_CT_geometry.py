# -*- coding: UTF-8 -*-
# File: examples/geometry/02_simple_CT_geometry.py

from ctsimu.geometry import *

# General CT parameters:
SOD = 250.0 # mm
SDD = 800.0 # mm

pixelSize = 0.2 # mm
pixelColumns = 2000
pixelRows = 1000

# Create a CT geometry object:
myCT = Geometry()

# Stage:
myCT.stage.center.set_x(SOD)

# Detector:
myCT.detector.center.set_x(SDD)
myCT.detector.set_size(
    pixels_u = pixelColumns,
    pixels_v = pixelRows,
    pitch_u = pixelSize,
    pitch_v = pixelSize
)

myCT.update() # calculates derived geometry parameters

print(myCT.info())