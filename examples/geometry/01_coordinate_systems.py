from ctsimu.geometry import *
from ctsimu.helpers  import *  # provides deg2rad()

mySpecimen = CoordinateSystem()

# Set position and orientation:
mySpecimen.center = Vector(250, 0, 0)
mySpecimen.u = Vector(0, -1,  0)
mySpecimen.v = Vector(0,  0, -1)
mySpecimen.w = Vector(1,  0,  0)

# Manipulate:
mySpecimen.translate(translationVector=Vector(5.2, 0, 4.3))
mySpecimen.rotateAroundU(angle=deg2rad(2.0))
mySpecimen.rotate(axis=Vector(1, 1, 1), angle=deg2rad(5.0))

print("My specimen's new location and orientation:")
print(mySpecimen)