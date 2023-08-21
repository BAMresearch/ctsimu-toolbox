# -*- coding: UTF-8 -*-
import math
from ctsimu.geometry import *

mySpecimen = CoordinateSystem()

# Set position and orientation:
mySpecimen.center = Vector(250, 0, 0)
mySpecimen.u = Vector(0, -1,  0)
mySpecimen.v = Vector(0,  0, -1)
mySpecimen.w = Vector(1,  0,  0)

# Manipulate:
mySpecimen.translate(translation_vector=Vector(5.2, 0, 4.3))
mySpecimen.rotate_around_u(angle=math.radians(2.0))
mySpecimen.rotate(axis=Vector(1, 1, 1), angle=math.radians(5.0))

print("My specimen's new location and orientation:")
print(mySpecimen)

"""
My specimen's new location and orientation:
Center: [255.2   0.    4.3]
u:      [ 0.04905096 -0.99746313 -0.05158783]
v:      [-0.01674544  0.05082147 -0.99856736]
w:      [ 0.99865589  0.04984455 -0.01421012]

"""