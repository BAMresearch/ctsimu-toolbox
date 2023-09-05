# -*- coding: UTF-8 -*-
# File: examples/scenario/06_reconstruction_metadata

from ctsimu.scenario import Scenario
s = Scenario("example.json")

# Information about RAW projection images:
s.metadata.output.projections.filename.set("example_%04d.raw")
s.metadata.output.projections.datatype.set("float32")
s.metadata.output.projections.byteorder.set("little")
s.metadata.output.projections.headersize.file.set(1024)

# ...or just TIFF files:
s.metadata.output.projections.filename.set("example_%04d.tif")

# Projection image size:
s.metadata.output.projections.dimensions.x.set(1500) # pixels
s.metadata.output.projections.dimensions.y.set(1000) # pixels
s.metadata.output.projections.pixelsize.x.set(0.1) # mm
s.metadata.output.projections.pixelsize.y.set(0.1) # mm

# Maximum free-beam intensity gray value:
s.metadata.output.projections.max_intensity.set(44000)

# No flat and dark field images; let's assume
# projections are already corrected:
s.metadata.output.projections.dark_field.filename.set(None)
s.metadata.output.projections.flat_field.filename.set(None)

# Information about the tomogram:
s.metadata.output.tomogram.filename.set("example_recon.raw")
s.metadata.output.tomogram.datatype.set("uint16")
s.metadata.output.tomogram.byteorder.set("little")

# Tomogram size:
s.metadata.output.tomogram.dimensions.x.set(1500) # voxels
s.metadata.output.tomogram.dimensions.y.set(1500) # voxels
s.metadata.output.tomogram.dimensions.z.set(1000) # voxels
s.metadata.output.tomogram.voxelsize.x.set(0.05) # mm
s.metadata.output.tomogram.voxelsize.y.set(0.05) # mm
s.metadata.output.tomogram.voxelsize.z.set(0.05) # mm

s.write_CERA_config(
    save_dir="cera_recon",
    basename="example",
    create_vgi=True
)

s.write_OpenCT_config(
    save_dir="openct_recon",
    basename="example",
    create_vgi=True,
    variant="free",
    abspaths=True
)