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

"""
The specimen's world coordinates:
Center: ( 250.0000000, -0.1744975,  4.9969541)
u:      ( 1.0000000,  0.0000000,  0.0000000)
v:      ( 0.0000000,  0.9993908,  0.0348995)
w:      ( 0.0000000, -0.0348995,  0.9993908)
"""