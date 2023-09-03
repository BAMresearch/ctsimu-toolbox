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

# Calculate the projection matrix:
P_openCT = myCT.projection_matrix(mode="OpenCT")
P_CERA   = myCT.projection_matrix(mode="CERA")

print("OpenCT:")
print(P_openCT)

print("CERA:")
print(P_CERA)