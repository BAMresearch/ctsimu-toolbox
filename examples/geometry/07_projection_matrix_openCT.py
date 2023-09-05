# -*- coding: UTF-8 -*-
# File: examples/geometry/07_projection_matrix_openCT.py

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

# Set up a new volume coordinate system
# as a standard coordinate system,
# relative to the stage coordinate system:
volume = CoordinateSystem()
volume.w.invert()  # mirror reconstruction volume

# Calculate the projection matrix:
P = myCT.projection_matrix(volumeCS = volume)

print("openCT projection matrix:")
print(P)

"""
openCT projection matrix:
[[ 0.     3.2    0.     0.   ]
 [ 0.     0.    -3.2    0.   ]
 [-0.004  0.     0.     1.   ]]
"""