# -*- coding: UTF-8 -*-
# File: examples/geometry/08_projection_matrix_example3.py

from ctsimu.geometry import *

# Set up a quick CT geometry:
myCT = Geometry()
myCT.stage.center.set_x(250)     # SOD in mm
myCT.detector.center.set_x(800)  # SDD in mm

# Set the detector size:
myCT.detector.set_size(
	pixels_u = 2000,
	pixels_v = 1000,
	pitch_u  = 0.2,
	pitch_v  = 0.2)

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
image.center.set_x(-myCT.detector.phys_width  / 2.0)
image.center.set_y(-myCT.detector.phys_height / 2.0)

# Our unit of the image CS is in px, so we need to
# scale the image CS basis vectors by the pixel size.
image.u.scale(myCT.detector.pitch_u)
image.v.scale(myCT.detector.pitch_v)

# Set up a new volume coordinate system,
# relative to the stage coordinate system:
volume = CoordinateSystem()

# Move the volume origin to the front lower right
# corner of the reconstruction volume:
volume.center.set_x(-volume_size_x * voxel_size / 2.0)
volume.center.set_y(-volume_size_y * voxel_size / 2.0)
volume.center.set_z(-volume_size_z * voxel_size / 2.0)

# Our unit of the volume CS is in voxels, so we need to
# scale the volume CS basis vectors by the voxel size.
volume.u.scale(voxel_size)
volume.v.scale(voxel_size)
volume.w.scale(voxel_size)

# Calculate the projection matrix:
P = myCT.projection_matrix(imageCS=image, volumeCS=volume)

print("My projection matrix:")
print(P)

"""
My projection matrix:
[[-3.33333333e-01  1.33333333e+00  0.00000000e+00  2.33333333e+03]
 [-1.66666667e-01  0.00000000e+00  1.33333333e+00  1.16666667e+03]
 [-3.33333333e-04  0.00000000e+00  0.00000000e+00  1.00000000e+00]]
"""