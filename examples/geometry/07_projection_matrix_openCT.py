# -*- coding: UTF-8 -*-
from ctsimu.geometry import *

# Set up a quick CT geometry:
myCT = Geometry()
myCT.stage.center.x    = 250  # SOD
myCT.detector.center.x = 800  # SDD

# Set the detector size:
myCT.detector.setSize(
	pixelsU = 2000,
	pixelsV = 1000,
	pitchU  = 0.2,
	pitchV  = 0.2)

# Set up a new image coordinate system
# as a standard coordinate system,
# relative to the detector coordinate system:
image = CoordinateSystem()

# Set up a new volume coordinate system
# as a standard coordinate system,
# relative to the stage coordinate system:
volume = CoordinateSystem()

# Calculate the projection matrix:
P = myCT.projectionMatrix(
	imageCS = image,
	volumeCS = volume)

print("openCT projection matrix:")
print(P)