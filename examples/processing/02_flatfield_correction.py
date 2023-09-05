# -*- coding: UTF-8 -*-
# File: examples/processing/02_flatfield_correction.py

from ctsimu.image import ImageStack
from ctsimu.processing.pipeline import Pipeline
from ctsimu.processing.flat_field import Step_FlatFieldCorrection

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
    filePattern="tetrahedron/projections/corrected/tetra_%04d.tif",
    dataType="float32"
)

# Create a pipeline:
pipe = Pipeline(
    inputFileStack=projectionFiles,
    outputFileStack=outputFiles
)

# 3 bright free-beam images for flat-field correction:
flatField = ImageStack(
    "tetrahedron/projections/tetra_flat_%04d.raw",
    width=150,
    height=150,
    dataType="uint16",
    byteOrder="little"
)

# One dark image:
darkField = ImageStack(
    "tetrahedron/projections/tetra_dark.raw",
    width=150,
    height=150,
    dataType="uint16",
    byteOrder="little"
)

# Create flat-field correction processing step:
ffCorrection = Step_FlatFieldCorrection(
    flatFileStack=flatField,
    darkFileStack=darkField,
    rescaleFactor=50000,
    offsetAfterRescale=0
)

# Add the processing step to the pipeline
# and run the processing:
pipe.addStep(ffCorrection)
pipe.run()