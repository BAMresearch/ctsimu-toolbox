# -*- coding: UTF-8 -*-
from ctsimu.geometry import *

# Set up a quick CT geometry:
myCT = Geometry()
myCT.stage.center.set_x(250)    # SOD
myCT.detector.center.set_x(800) # SDD

# Calculate the projection matrix:
P = myCT.projection_matrix()

print("Projection Matrix:")
print(P)

"""
Projection Matrix:
[[ 0.     3.2    0.     0.   ]
 [ 0.     0.     3.2    0.   ]
 [-0.004  0.     0.     1.   ]]
"""