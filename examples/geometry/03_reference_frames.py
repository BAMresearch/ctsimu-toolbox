from ctsimu.geometry import *
from ctsimu.helpers  import *  # provides deg2rad()

world = CoordinateSystem()

# Set up a quick CT geometry with a tilted stage axis:
myCT = Geometry()
myCT.stage.center.x = 250  # SOD
myCT.stage.rotateAroundU(angle = deg2rad(2.0))
myCT.detector.center.x = 800  # SDD

# Assume a specimen in the (tilted) stage
# coordinate system, shifted 5 mm "upwards"
# along the axis of rotation:
mySpecimen = CoordinateSystem()
mySpecimen.translateZ(5.0)

# Change the specimen's reference frame to
# the world coordinate system:
mySpecimen.changeReferenceFrame(
	fromCS = myCT.stage,
	toCS = world)

print("The specimen's world coordinates:")
print(mySpecimen)