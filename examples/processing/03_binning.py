# -*- coding: UTF-8 -*-
# File: examples/processing/03_binning.py

from ctsimu.image import ImageStack
from ctsimu.processing.pipeline import Pipeline
from ctsimu.processing.binning import Step_Binning

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
    filePattern="tetrahedron/projections/binned/tetra_%04d.tif",
    dataType="float32"
)

# Create a pipeline:
pipe = Pipeline(
    inputFileStack=projectionFiles,
    outputFileStack=outputFiles
)

# Create binning step:
binning = Step_Binning(
    binSizeX=2,
    binSizeY=2,
    binningOperation='mean'
)

# Add the processing step to the pipeline
# and run the processing:
pipe.addStep(binning)
pipe.run()