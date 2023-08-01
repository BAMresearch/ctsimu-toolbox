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
A **Python 3** environment (Python 3.8 or higher) is required. The Python distribution that was used to develop and test the toolbox is *Anaconda 3* (available for Windows, macOS and Linux), and therefore a recommendation if you don't know where to start. The following Python packages are required as well. They usually come with a Python environment, or can be installed there easily:

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

# Usage - Example of a test implementation

In the CTSimU toolbox, the tests are performed using JSON files in Python. There are two files (where X. is the current test, detector, subtest, and the date of the last modification of the JSON file):

- `X..json`: defines environment, geometry, materials, source, detector properties, among others, and - references files in the same folder for the description of the object and spectrum: STL file, CSV and XRS files. You can find the `X..json` file in the path: `ctsimu-toolbox\ctsimu\ctsimu_evaluations\scenarios` or on the [website](https://github.com/BAMresearch/ctsimu-toolbox/tree/main/ctsimu/ctsimu_evaluations/scenarios).
- `X._metadata.json`: Contains information about the simulation results (type and location) and is created via a simulation program (e.g. aRTist). With this result file, the test can be started and executed in the CTSimU toolbox.

The following is a step-by-step guide using the example of the double wire test to determine the detector unsharpness.

## Simulation for CTSimU in a*RT*ist

In this example, a double wire bar measurement is presented to determine the detector unsharpness (the test name is 2D-DW-1). Theoretically, any other simulation programme could be used. We will use the simulation software aRTist to create the necessary projection images for the test. To do this, you need to download and install the appropriate CTSimU module for a*RT*ist. Instructions on how to install and operate the module, as well as the relevant ARTP installation file, can be found on the [BAM research website](https://github.com/BAMresearch/ctsimu-artist-module). 

### Step-by-Step Guide

Step 1:	The prerequisite is that the scenario in this example is configured correctly in a*RT*ist.

- For this the following files should be in one folder:
    - `2D-DW-1_Detektor1_2021-05-26v01r01dp.json` file, the STL of the double wire path and the two files (TXT and XRS) for the description of the spectrum. 

Step 2:	Open aRTist and load the JSON file:
- Load into to the CTSimU-Module the 2D-DW-1_Detector1_2021-05-26v01r01dp.json file, which contains the necessary simulation details, including geometry and spectrum descriptions.

Step 3: Run the scenario and generate projections:

- Start the scenario by clicking on **Run scenario** in a*RT*ist
- Once the simulation is complete, the finished projection images will be visible in the Image Viewer (see figure below).
- In the folder where the JSON file is located, a new folder named `2D-DW-1_Detector1_2021-05-26v01r01dp` will be created, containing a subfolder named `projections`.

![aRTist scene for DW_1-test](C:/Users/jjanczyn/Documents/ctsimu-toolbox/docs/pictures/scene_DW_1.png "aRTist scene for DW_1-test")

Step 4: Review the contents of the `projections` folder. The folder will contain the following files:

- Two TIF files: a projection image of the double wire bar according to EN462-5 and the corresponding flatfield image (background noise image). 
-  Python program: `2D-DW-1_Detektor1_2021-05-26v01r01dp_flat.py` which includes the code for the flatfield correction.
- CTSimU metadata file: 2D-DW-1_Detector1_2021-05-26v01r01dp_metadata.json, which contains all important information about the two projections and is used for the test execution.

## CTSimU Toolbox Test Execution from the Python Console

The test "Double Wire Ridge - Detector Unsharpness" can be found at [BAM research website in CTSimu Evaluation](https://github.com/BAMresearch/ctsimu-toolbox/tree/main/ctsimu/ctsimu_evaluations) with the name `test2D_DW_1.py`. In this Python file, the test name "2D-DW-1" and thus also input name for the Python console can be found in line 43 (as testName="2D-DW-1"):

```python
36 class Test2D_DW_1(generalTest):
37	""" CTSimU test 2D-DW-1: Detector Unsharpness. """
38
39	def __init__(self, resultFileDirectory=".", name=None, rawOutput=False):
40
41		generalTest.__init__(
42			self,
43			testName="2D-DW-1",
44			name=name,
45			nExpectedRuns=2,
46			resultFileDirectory=resultFileDirectory,
47          rawOutput=rawOutput)
```

The names of the sub-tests can be read in line 80 with "SR75" and in line 82 with "SR150":

```python
77 def prepareRun(self, i):
78    if i < len(self.subtests):
79        self.jsonScenarioFile = "scenarios/2D-DW-1_Detektor1_2021-05-26v01r01dp.json"
80        if self.subtests[i] == "SR75":
81            self.jsonScenarioFile = "scenarios/2D-DW-1_Detektor1_2021-05-26v01r01dp.json"
82        elif self.subtests[i] == "SR150":
83            self.jsonScenarioFile = "scenarios/2D-DW-1_Detektor2_2021-05-26v01r01dp.json"
```

This is similarly structured in the other tests. In some of them there are no sub-tests.

The tests SR75 and SR150 are the same only they are stored according to their names, here according to the expected spatial resolution. Since we know in advance in which range the resolution is, we can directly choose the test with the appropriate name. In this way, a suitable name is created for the resulting data. Detector1 should have a spatial resolution of iSRb = 75μm and Detector2 should have iSRb = 150μm.

### Step-by-Step Guide

Running the 2D-DW-1 test in a Python console:

To run the 2D-DW-1 test, the generated X._metadata.json file from a*RT*ist is used. The test commands are then executed in the Python console in the following way: 

Step 1:	Open the Python Console

Step 2: Import the CTSimU Toolbox Library:

- In the Python console, import the CTSimU Toolbox library with the following command:

```python
from ctsimu.toolbox import Toolbox
```

Step 3:	Run a Partial Test:

- To execute a partial test for Detector1 with a spatial resolution of 75 microns (SR75), use the following command:

```python
Toolbox("2D-DW-1", SR75 = "\\...\\ctsimu-toolbox\\ctsimu\\ctsimu_evaluations\\scenarios\\2D-DW-1_Detektor1_2021-05-26v01r01dp\\projections\\2D-DW-1_Detektor1_2021-05-26v01r01dp_metadata.json")
```

Step 4:	Run Both Partial Tests Sequentially:

- To execute both partial tests for Detector1 (SR75) and Detector2 (SR150) sequentially, use the following command:

```python
Toolbox("2D-DW-1", SR75 = "\\...\\ctsimu-toolbox\\ctsimu\\ctsimu_evaluations\\scenarios\\2D-DW-1_Detektor1_2021-05-26v01r01dp\\projections\\2D-DW-1_Detektor1_2021-05-26v01r01dp_metadata.json", SR150 = "\\...\\ctsimu-toolbox\\ctsimu\\ctsimu_evaluations\\scenarios\\2D-DW-1_Detektor2_2021-05-26v01r01dp\\projections\\2D-DW-1_Detektor2_2021-05-26v01r01dp_metadata.json")
```

Step 5:	After executing the test, two folders will be saved in the same path as the `X._metadata.json` file:

- Folder `corrected`: 
  - Contains the flatfield-corrected image.
- Folder `2D-DW-1-results`: Contains grey value profiles and grey value distribution data for each detector's basic spatial resolution (iSRb) with corresponding measurement values (wire spacing [mm] and modulation depth [%]).
  - The values of the interpolation function are saved in TXT files.
  - Corresponding plots are displayed as PNG images (see figures below), including a histogram (grey value distribution over the horizontal distance in pixels) and a plot of the spatial resolution of the investigated detector with the interpolation function and the determined iSRb value (modulation depth [%] over the wire distance [mm]).

![2D-DW-1 SR75 results](C:/Users/jjanczyn/Documents/ctsimu-toolbox/docs/pictures/2D-DW-1_SR75_results.png "2D-DW-1 SR75 results")

![2D-DW-1 SR150 results](C:/Users/jjanczyn/Documents/ctsimu-toolbox/docs/pictures/2D-DW-1_SR150_results.png "2D-DW-1 SR150 results")

By using these steps, you can run the 2D DW-1 test successfully in the CTSimU Toolbox from the Python console and get the results and data you want.

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