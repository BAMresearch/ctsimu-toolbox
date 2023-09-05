# -*- coding: UTF-8 -*-
# File: examples/geometry/09_projection_matrix_full_CT.py

import math
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

myCT.update() # signal that we made manual changes
myCT.store()  # backup the initial configuration

# Scan configuration:
projections = 3000   # number of projections or angular steps
scan_range  = 360.0  # degrees. One full CT rotation.

# We assume that the projections are stored in single TIFF image files,
# sequentially numbered with four digits, starting at "img_0000.tif".
projection_filename = "img_{:04d}.tif"   # for OpenCT
projection_file_pattern = "img_%04d.tif" # for CERA

# The following two lists will store the projection matrices
# for openCT and for CERA:
matrices_openCT = []
matrices_CERA   = []

# For openCT, we also need to create a list of projection file names:
projection_filenames = []

# Loop over each frame:
for p in range(projections):
	# Restore the initial configuration from the backup,
	# i.e. the situation before the stage was rotated:
	myCT.restore()

	# Rotate the stage to its current angle:
	current_angle = float(p) * float(scan_range) / float(projections)
	myCT.stage.rotate_around_w(angle=math.radians(current_angle))
	myCT.update()

	# Calculate a projection matrix for this frame:
	P_openCT = myCT.projection_matrix(mode="OpenCT")
	P_CERA   = myCT.projection_matrix(mode="CERA")

	# Add to list of projection matrices:
	matrices_openCT.append(P_openCT)
	matrices_CERA.append(P_CERA)

	# Store the current projection filename for openCT:
	projection_filenames.append(projection_filename.format(p))

# Restore CT setup for frame zero:
myCT.restore()

# openCT configuration:
# ----------------------
# Write the openCT configuration file, including the projection matrices:
create_OpenCT_config(
	geo=myCT,
	filename="example_09/recon_openCT.json",
	projection_files=projection_filenames,
	matrices=matrices_openCT,
	volumename="recon_openCT"
)

# CERA configuration:
# -------------------
# Write the CERA configuration file, including the projection matrices:
create_CERA_config(
	geo=myCT,
	total_angle=scan_range,
	projection_file_pattern=projection_file_pattern,
	matrices=matrices_CERA,
	basename="recon_CERA",
	save_dir="example_09",
    i0max=44000  # maximum free-beam intensity
)