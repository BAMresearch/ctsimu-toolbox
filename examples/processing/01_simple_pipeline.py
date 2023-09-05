# -*- coding: UTF-8 -*-
# File: examples/processing/01_simple_pipeline.py

from ctsimu.image import ImageStack
from ctsimu.processing.pipeline import Pipeline

# All projection images are in one raw chunk:
projectionFiles = ImageStack(
    filePattern="tetrahedron/projections/tetra_projections.raw",
    width=150,
    height=150,
    slices=21,
    dataType="uint16",
    byteOrder="little"
)

# Save as individual TIFF images:
outputFiles = ImageStack(
    filePattern="tetrahedron/projections/tetra_%04d.tif",
    dataType="float32"
)

# Create and run pipeline:
pipe = Pipeline(
    inputFileStack=projectionFiles,
    outputFileStack=outputFiles
)
pipe.run()