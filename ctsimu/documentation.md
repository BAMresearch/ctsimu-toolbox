Features
========

The toolbox provides the following **main submodules,** each designed
for a different task.

Quick access
------------
The **`ctsimu.toolbox`** sub-module provides a simplified interface to the most
commonly used functionality. It reads [CTSimU scenario description files]
and [metadata files] to carry out the following tasks:

* Flat-field correction,
* Create reconstruction configuration files,
* Standardize or update scenario description files to the latest file format version,
* Run recursive post-processing on whole directory trees (batch processing).

Image conversion & processing
------------------------------
The **`ctsimu.processing`** sub-module provides image processing pipelines that
can be used to convert stacks of projection images into different formats.
The following properties can be converted:

* Image format conversion (RAW, TIFF),
* Data type conversion (int, float),
* Representation conversion (multi-file slice stacks, single-file RAW volume chunks).

Additional processing steps can be added to a pipeline: this will alter the
images before they are saved. Currently, there are processing steps for

* Flat-field and dark-field correction with provided gain/offset files,
* Pixel binning,
* Median filter,
* Transformations: rotations, flipping.

Image quality assessment
------------------------
The **`ctsimu.image_analysis`** sub-module currently provides functions for:

* Line profiles,
* `ctsimu.image_analysis.mtf`: Calculation of the modulation transfer function (MTF) from an edge image,
* `ctsimu.image_analysis.isrb`: Determination of the interpolated basic spatial resolution (iSRb) from a duplex wire image.

Geometry tools
--------------
The **`ctsimu.geometry`** sub-module provides tools for

* Geometric representation and manipulation of a CT scan,
* Calculation of projection matrices,
* Creation of reconstruction configurations for openCT (VGSTUDIO) and CERA,
  even though it is easier to use the `ctsimu.scenario` sub-module for this.

[CTSimU scenario description files]: https://bamresearch.github.io/ctsimu-scenarios/
[metadata files]: https://bamresearch.github.io/ctsimu-scenarios/metadata.html

Scenario handling
-----------------
The **`ctsimu.scenario`** sub-module is used to handle complete
CTSimU scenarios:

* Read, write, explore and manipulate [CTSimU scenario description files],
* Calculate coordinate systems for each frame,
* Generate reconstruction configuration files.

[CTSimU scenario description files]: https://bamresearch.github.io/ctsimu-scenarios/

2D test evaluations
-------------------
The **`ctsimu.evaluation`** sub-module provides scripts and documentation
for CTSimU 2D test evaluations.


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
* Anthony Orth