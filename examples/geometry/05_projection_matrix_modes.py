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

# Calculate the projection matrix:
P_openCT = myCT.projectionMatrix(mode="openCT")
P_CERA   = myCT.projectionMatrix(mode="CERA")

print("openCT:")
print(P_openCT)

print("CERA:")
print(P_CERA)