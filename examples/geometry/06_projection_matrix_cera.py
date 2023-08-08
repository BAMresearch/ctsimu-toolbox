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
image.w.scale(-1.0)

# Calculate the projection matrix:
P = myCT.projection_matrix(imageCS=image)

print("CERA projection matrix:")
print(P)