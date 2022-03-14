from ctsimu.geometry import *

# Set up a quick CT geometry:
myCT = Geometry()
myCT.stage.center.x    = 250  # SOD
myCT.detector.center.x = 800  # SDD

# Calculate the projection matrix:
P = myCT.projectionMatrix()

print("Projection Matrix:")
print(P)

"""
Projection Matrix:
[ 0.0000000  -3.2000000   0.0000000   0.0000000  ]
[ 0.0000000   0.0000000   3.2000000   0.0000000  ]
[ 0.0040000   0.0000000   0.0000000   1.0000000  ]
"""