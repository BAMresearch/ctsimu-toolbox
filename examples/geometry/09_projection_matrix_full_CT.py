# -*- coding: UTF-8 -*-
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

myCT.update() # signals that we made manual changes
myCT.store()  # backups the initial configuration

# Scan configuration:
projections = 3000   # number of projections or angular steps
scan_range  = 360.0  # degrees. One full CT rotation.

# We assume that the projections are stored in single TIFF image files,
# sequentially numbered with four digits, starting at "img_0000.tif".
projectionFilename    = "img_{:04d}.tif" # for openCT
projectionFilePattern = "img_%04d.tif"   # for CERA

# The following two lists will store the projection matrices
# for openCT and for CERA:
matrices_openCT = []
matrices_CERA   = []

# For openCT, we also need to create a list of projection file names:
projectionFilenames = []

# Loop over each frame:
for p in range(projections):
	# Restore the initial configuration from the backup,
	# i.e. the situation before the stage was rotated:
	myCT.restore()

	# Rotate the stage to its current angle:
	current_angle = float(p) * float(scan_range) / float(projections)
	myCT.stage.rotate_around_w(angle = deg2rad(current_angle))
	myCT.update()

	# Calculate a projection matrix for this frame:
	P_openCT = myCT.projection_matrix(mode="openCT")
	P_CERA   = myCT.projection_matrix(mode="CERA")

	# Add to list of projection matrices:
	matrices_openCT.append(P_openCT)
	matrices_CERA.append(P_CERA)

	# Store the current projection filename for openCT:
	projectionFilenames.append(projectionFilename.format(p))


# openCT configuration:
# ----------------------
# We need the bounding box dimensions of the reconstruction volume
# in mm:

voxelSize = 0.0625
bounding_box_x = voxelSize * myCT.detector.pixels_u
bounding_box_y = voxelSize * myCT.detector.pixels_u
bounding_box_z = voxelSize * myCT.detector.pixels_v

# Write the openCT configuration file, including the projection matrices:
write_openCT_file(
	geo=myCT,
	totalAngle=scan_range,
	boundingBoxX=bounding_box_x,
	boundingBoxY=bounding_box_y,
	boundingBoxZ=bounding_box_z,
	matrices=matrices_openCT,
	volumename="recon_openCT",
	filename="recon_openCT.json",
    projectionFilenames=projectionFilenames
)

# CERA configuration:
# -------------------
write_cera_config(
	geo=myCT,
	totalAngle=scan_range,
	projectionFilePattern=projectionFilePattern,
	matrices=matrices_CERA,
	basename="recon_CERA",
	voxelsX=myCT.detector.pixels_u,
    voxelsY=myCT.detector.pixels_u,
    voxelsZ=myCT.detector.pixels_v,
    i0max=44000  # the average free-beam intensity
)