from ctsimu.geometry import *

# Set up a quick CT geometry:
myCT = Geometry()
myCT.stage.center.x    = 250  # SOD in mm
myCT.detector.center.x = 800  # SDD in mm

# Set the detector size:
myCT.detector.setSize(
	pixelsU = 2000,
	pixelsV = 1000,
	pitchU  = 0.2,
	pitchV  = 0.2)

# Define the size of the reconstruction volume:
volume_size_x = 2000 # voxels
volume_size_y = 2000 # voxels
volume_size_z = 1000 # voxels
voxel_size = 0.0625  # mm/voxel

# Set up a new image coordinate system,
# relative to the detector coordinate system:
image = CoordinateSystem()

# Move the image origin to the upper left corner
# of the detector coordinate system:
image.center.x = -myCT.detector.physWidth  / 2.0
image.center.y = -myCT.detector.physHeight / 2.0

# Our unit of the image CS is in px, so we need to
# scale the image CS basis vectors by the pixel size.
image.u.scale(myCT.detector.pitchU)
image.v.scale(myCT.detector.pitchV)

# Set up a new volume coordinate system,
# relative to the stage coordinate system:
volume = CoordinateSystem()

# Move the volume origin to the front lower right
# corner of the reconstruction volume:
volume.center.x = -volume_size_x * voxel_size / 2.0
volume.center.y = -volume_size_y * voxel_size / 2.0
volume.center.z = -volume_size_z * voxel_size / 2.0

# Our unit of the volume CS is in voxels, so we need to
# scale the volume CS basis vectors by the voxel size.
volume.u.scale(voxel_size)
volume.v.scale(voxel_size)
volume.w.scale(voxel_size)

# Calculate the projection matrix:
P = myCT.projectionMatrix(imageCS=image, volumeCS=volume)

print("My projection matrix:")
print(P)