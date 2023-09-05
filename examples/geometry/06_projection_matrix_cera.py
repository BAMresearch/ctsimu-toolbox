# -*- coding: UTF-8 -*-
# File: examples/geometry/06_projection_matrix_cera.py

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

# Set up a new image coordinate system,
# relative to the detector coordinate system:
image = CoordinateSystem()

# CERA places the origin of the image CS in the center
# of the lower left pixel of the projection image.
image.center.set_x(-(myCT.detector.phys_width  / 2.0) + 0.5*myCT.detector.pitch_u)
image.center.set_y( (myCT.detector.phys_height / 2.0) - 0.5*myCT.detector.pitch_v)

# CERA's unit of the image CS is in px, so we need to
# scale the image CS basis vectors by the pixel size.
# Also, v points up instead of down. This also flips
# the w axis to keep a right-handed coordinate system.
image.u.scale( myCT.detector.pitch_u)
image.v.scale(-myCT.detector.pitch_v)

# CERA's volume coordinate system is equivalent to the CTSimU stage coordinate
# coordinate system, but flipped vertically. Therefore, we need to
# invert the volume's w axis.
volume = CoordinateSystem()
volume.w.invert()

# Calculate the projection matrix:
P = myCT.projection_matrix(imageCS=image, volumeCS=volume)

print("CERA projection matrix:")
print(P)

"""
CERA projection matrix:
[[-3.998e+00  1.600e+01  0.000e+00  9.995e+02]
 [-1.998e+00  0.000e+00  1.600e+01  4.995e+02]
 [-4.000e-03  0.000e+00  0.000e+00  1.000e+00]]
"""