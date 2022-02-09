Features
========

In its current version, the toolbox can be used for the following tasks:

Image conversion
----------------
For image conversion, and "empty" `ctsimu.processing.pipeline` can be used. The following properties can be converted:

* Image format conversion (RAW, TIFF),
* Data type conversion (int, float),
* Representation conversion (multi-file slice stacks, single-file RAW volume chunks).

 For details, refer to the documentation of the `ctsimu.processing` module.

Image processing
----------------
For image processing, processing steps can be added to the pipeline. Currently, there are processing steps for

* Flat-field and dark-field correction with provided gain/offset files,
* Pixel binning,
* Median filter,
* Transformations: rotations, flipping.

You can find out more in the documentation of the `ctsimu.processing` module.

Image quality assessment
------------------------
The `ctsimu.image_analysis` module currently provides functions for:

* Line profiles,
* `ctsimu.image_analysis.mtf`: Calculation of the modulation transfer function (MTF) from an edge image,
* `ctsimu.image_analysis.isrb`: Determination of the interpolated basic spatial resolution (iSRb) from a duplex wire image.

Geometry tools
--------------
The `ctsimu.geometry` module provides tools for

* Virtual representation and manipulation of a CT scene in Python,
* Importing CT geometries from [CTSimU scenario description files],
* Creation of reconstruction configurations for openCT (VGSTUDIO) and CERA,
* Calculation of projection matrices.

[CTSimU scenario description files]: https://bamresearch.github.io/ctsimu-scenarios


Getting started
===============

Requirements
------------
A **Python 3** environment (Python 3.9 or higher) is required. The Python distribution that was used to develop and test the toolbox is *Anaconda 3* (available for Windows, macOS and Linux), and therefore a recommendation if you don't know where to start. The following Python packages are required as well. They usually come with a Python environment, or can be installed there easily:

* NumPy,
* SciPy subpackages *ndimage, optimize, stats, signal* and
* Matplotlib (only required if you want plots of evaluation results).


Installation
------------

The file :guilabel:`ctsimu_1.4.4.zip` contains the complete Python package.

To use it **without installation,** copy the folder `ctsimu` from the zip file to the directory where you want to run the Python scripts that need to import the package modules.

If you choose to **install** the package in your Python environment, you can use [pip]. For example, you can run the following command in your *Anaconda Prompt*:

	pip install ctsimu_1.4.4.zip

[pip]: https://pip.pypa.io


About
=====

The CTSimU toolbox was developed for the [CTSimU] project for its needs in reproducible image processing and serves as the project's reference software implementation. It provides useful tools for many CT-related tasks, beyond the aims of the CTSimU project.

The software is released under the **Apache 2.0 license,** its source code is available [on GitHub].

[on GitHub]: https://github.com/BAMresearch/ctsimu-toolbox

Acknowledgements
----------------
This work was funded through the project [CTSimU] (*Radiographic Computed Tomography Simulation for Measurement Uncertainty Evaluation*, WIPANO project 03TNH026A).

[WIPANO] projects are financed by the German [Federal Ministry for Economic Affairs and Climate Action] and managed by [Project Management Jülich].

[CTSimU]: https://www.ctsimu.forschung.fau.de/
[WIPANO]: https://www.innovation-beratung-foerderung.de/INNO/Navigation/DE/WIPANO/wipano.html
[Federal Ministry for Economic Affairs and Climate Action]: https://www.bmwi.de/
[Project Management Jülich]: https://www.ptj.de/

Contributors
------------
The following people contributed to code and documentation of the toolbox.

* David Plotzki
* Bendix Hartlaub
* Florian Wohlgemuth
* Tamara Reuter
* Fabrício Borges de Oliveira
* Carsten Bellon