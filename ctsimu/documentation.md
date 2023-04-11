Features
========

In its current version, the toolbox can be used for the following tasks:

Image conversion
----------------
For image conversion, an "empty" `ctsimu.processing.pipeline` can be used. The following properties can be converted:

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

If you want to **install** the package in your Python environment, you can use [pip]. For example, you can run the following command in your *Anaconda Prompt* to make the toolbox available:

	pip install ctsimu

[pip]: https://pip.pypa.io

To use the package **without installation**, you need to download the package manually. You have the following three options:

* Download the package [from PyPi]. You will get a zipped file called `ctsimu-X.X.X.tar.gz` (where X.X.X is the current version number). Open or unpack this file to get to its contents.
* Download the repository [from GitHub]: press the *clone* button to download a ZIP file. Unpack it or open it to see its contents.
* You can also clone the repository from GitHub via the terminal:

	`git clone https://github.com/BAMresearch/ctsimu-toolbox.git`

From any of these three options, you will get the complete package source code. It contains a folder called `ctsimu`. If you want to use the toolbox from your Python scripts without installing it, you need to copy the folder `ctsimu` to the location of your Python scripts to make the package available.


[from GitHub]: https://github.com/BAMresearch/ctsimu-toolbox
[from PyPi]: https://pypi.org/project/ctsimu/

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