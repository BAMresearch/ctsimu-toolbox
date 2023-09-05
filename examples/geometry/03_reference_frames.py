# -*- coding: UTF-8 -*-
# File: examples/geometry/03_reference_frames.py

import math
from ctsimu.geometry import * # also contains ctsimu_world

# Set up a quick CT geometry with a tilted stage axis:
myCT = Geometry()
myCT.stage.center.set_x(250)  # SOD
myCT.stage.rotate_around_u(angle = math.radians(2.0))
myCT.detector.center.set_x(800)  # SDD

# Assume a specimen in the (tilted) stage
# coordinate system, shifted 5 mm "upwards"
# along the axis of rotation:
mySpecimen = CoordinateSystem()
mySpecimen.translate_z(5.0)

# Change the specimen's reference frame to
# the world coordinate system:
mySpecimen.change_reference_frame(
	cs_from = myCT.stage,
	cs_to = ctsimu_world
)

print("The specimen's world coordinates:")
print(mySpecimen)

"""
The specimen's world coordinates:
Center: ( 250.0000000, -0.1744975,  4.9969541)
u:      ( 1.0000000,  0.0000000,  0.0000000)
v:      ( 0.0000000,  0.9993908,  0.0348995)
w:      ( 0.0000000, -0.0348995,  0.9993908)
"""