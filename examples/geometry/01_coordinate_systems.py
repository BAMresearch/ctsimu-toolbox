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

"""
My specimen's new location and orientation:
Center: ( 255.2000000,  0.0000000,  4.3000000)
u:      ( 0.0490510, -0.9974631, -0.0515878)
v:      (-0.0167454,  0.0508215, -0.9985674)
w:      ( 0.9986559,  0.0498445, -0.0142101)
"""