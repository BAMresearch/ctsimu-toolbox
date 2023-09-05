# Image conversion

## Pipeline concept

![Example for a processing pipeline.](../pictures/pipeline.png "Example for a processing pipeline.")

The required minimum set of parameters are an input file pattern and an output file pattern. As illustrated in the following picture, a file pattern is a string that specifies a single file name (for a single 2D image or a raw chunk of multiple images), or a stack of 2D images using the `%[N]d` placeholder for a number that contains `[N]` digits, e.g. `%04d` for the numbers from `0000` to `9999`.

![File stacks and their file patterns.](../pictures/filestacks.png "File stacks and their file patterns.")

## Input and output files

Any kind of image file (input, output, flat fields) can be specified by creating `ctsimu.image.ImageStack` objects. These objects will later be passed to the processing pipeline or to processing step configurations.

In the following example, a 16 bit unsigned integer raw data volume file that contains 21 slices of 150x150 pixels is converted into a stack of 21 TIFF files, each with a 32 bit float data type.

```python
.. include:: ../../examples/processing/01_simple_pipeline.py
```

For TIFF files, the data type is determined from the header information when reading existing files, but the given value is obeyed when they are saved. This allows to convert the data type from one to another. For RAW files, this information must be provided for all input. Upon saving, RAW and TIFF output files automatically get the input file parameters if no output parameters are provided.

The minimum set of parameters required to set up a `ctsimu.processing.pipeline.Pipeline` are the input and output images. This way, the pipeline can be used to convert images from one format to another. Any parameters that are not specified for the output images (such as data type and byte order) will be taken from the input file parameters.


# Image processing

## Flat-field correction

Images can be processed by adding a processing `ctsimu.processing.step.Step` (or
several to create a processing chain) to the `ctsimu.processing.pipeline.Pipeline`.

The following example performs a dark-field and flat-field correction using
existing dark and bright free-beam images. It uses the processing step called
`ctsimu.processing.flat_field.Step_FlatFieldCorrection`.

```python
.. include:: ../../examples/processing/02_flatfield_correction.py
```

If a stack of multiple images is provided either for flat or dark images, the mean image of this stack will be used for the flat/dark field correction.

A `rescaleFactor` can be provided to rescale the projection gray values after they have been normalized. This is especially necessary if your output files will have an integer data type, as the gray values will typically be in the interval [0, 1] after a flat field correction.

## Binning

For binning, the bin size has to be specified in x and y direction,
as well as the binning operation (possible values: `"mean"`, `"min"`, `"max"`, `"sum"`).

```python
.. include:: ../../examples/processing/03_binning.py
```