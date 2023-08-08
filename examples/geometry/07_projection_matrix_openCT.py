# -*- coding: UTF-8 -*-
from ctsimu.geometry import *

# Set up a quick CT geometry:
myCT = Geometry()
myCT.stage.center.set_x(250)    # SOD
myCT.detector.center.set_x(800) # SDD

# Set the detector size:
myCT.detector.set_size(
	pixels_u = 2000,
	pixels_v = 1000,
	pitch_u  = 0.2,
	pitch_v  = 0.2)

# Set up a new image coordinate system
# as a standard coordinate system,
# relative to the detector coordinate system:
image = CoordinateSystem()

# Set up a new volume coordinate system
# as a standard coordinate system,
# relative to the stage coordinate system:
volume = CoordinateSystem()

# Calculate the projection matrix:
P = myCT.projection_matrix(
	imageCS = image,
	volumeCS = volume)

print("openCT projection matrix:")
print(P)