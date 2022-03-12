from ctsimu.geometry import *

# General CT parameters:
SOD = 250.0 # mm
SDD = 800.0 # mm

pixelSize = 0.2 # mm
pixelColumns = 2000
pixelRows = 1000

# Create a CT geometry object:
myCT = Geometry()

# Stage:
myCT.stage.center.x = SOD

# Detector:
myCT.detector.center.x = SDD
myCT.detector.setSize(
    pixelsU = pixelColumns,
    pixelsV = pixelRows,
    pitchU = pixelSize,
    pitchV = pixelSize
    )

myCT.update() # calculates derived geometry parameters

print(myCT.info())