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

# Set up a new image coordinate system,
# relative to the detector coordinate system:
image = CoordinateSystem()

# CERA places the origin of the image CS in the center
# of the lower left pixel of the projection image.
image.center.x = -(myCT.detector.physWidth  / 2.0) + 0.5*myCT.detector.pitchU
image.center.y =  (myCT.detector.physHeight / 2.0) - 0.5*myCT.detector.pitchV

# CERA's unit of the image CS is in px, so we need to
# scale the image CS basis vectors by the pixel size.
# Also, v points up instead of down. This also flips
# the w axis to keep a right-handed coordinate system.
image.u.scale( myCT.detector.pitchU)
image.v.scale(-myCT.detector.pitchV)
image.w.scale(-1.0)

# Calculate the projection matrix:
P = myCT.projectionMatrix(imageCS=image)

print("CERA projection matrix:")
print(P)