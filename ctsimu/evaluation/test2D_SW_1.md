In this scenario, the imaging characteristics of specific scintillators is tested for different X-ray source spectra and two different object materials. Two spherical step wedges are used, one made of aluminum, the other one of iron. Each wedge features nine steps of different thickness. The test compares their gray values (as well as the free beam gray value) for each combination of scintillator, material and spectrum (both monochromatic and polychromatic) to a corresponding ideal image taken with an *ideal* detector, i.e., one whose gray values correspond linearly to the intensity of incident radiation.

The following table lists all the scenarios that need to be simulated for this test, along with the keyword that needs to be used to correctly identify each projection metadata file for the toolbox.

| Keyword       | Material | Spectrum       | Scintillator                           |
| :------------ | :------- | :------------- | :------------------------------------- |
| `D1_Al_50`    | Al       | 50 keV mono    | ideal                                  |
| `D1_Al_120`   | Al       | 120 kV poly    | ideal                                  |
| `D1_Fe_100`   | Fe       | 100 keV mono   | ideal                                  |
| `D1_Fe_200`   | Fe       | 200 kV poly    | ideal                                  |
| `D2_Al_50`    | Al       | 50 keV mono    | 200 &mu;m CsI                          |
| `D2_Al_120`   | Al       | 120 kV poly    | 200 &mu;m CsI                          |
| `D2_Fe_100`   | Fe       | 100 keV mono   | 200 &mu;m CsI                          |
| `D2_Fe_200`   | Fe       | 200 kV poly    | 200 &mu;m CsI                          |
| `D3_Al_50`    | Al       | 50 keV mono    | 500 &mu;m CsI                          |
| `D3_Al_120`   | Al       | 120 kV poly    | 500 &mu;m CsI                          |
| `D3_Fe_100`   | Fe       | 100 keV mono   | 500 &mu;m CsI                          |
| `D3_Fe_200`   | Fe       | 200 kV poly    | 500 &mu;m CsI                          |
| `D4_Al_50`    | Al       | 50 keV mono    | 500 &mu;m Gd<sub>2</sub>O<sub>2</sub>S |
| `D4_Al_120`   | Al       | 120 kV poly    | 500 &mu;m Gd<sub>2</sub>O<sub>2</sub>S |
| `D4_Fe_100`   | Fe       | 100 keV mono   | 500 &mu;m Gd<sub>2</sub>O<sub>2</sub>S |
| `D4_Fe_200`   | Fe       | 200 kV poly    | 500 &mu;m Gd<sub>2</sub>O<sub>2</sub>S |

An example call for this test scenario looks like this:

```python
from ctsimu.toolbox import Toolbox
Toolbox("2D-SW-1",
    D1_Al_50  = "2D-SW-1_Al_Detector1_50keV-mono_metadata.json",
    D1_Al_120 = "2D-SW-1_Al_Detector1_120kV-poly_metadata.json",
    D1_Fe_100 = "2D-SW-1_Fe_Detector1_100keV-mono_metadata.json",
    D1_Fe_200 = "2D-SW-1_Fe_Detector1_200kV-poly_metadata.json",
    D2_Al_50  = "2D-SW-1_Al_Detector2_50keV-mono_metadata.json",
    D2_Al_120 = "2D-SW-1_Al_Detector2_120kV-poly_metadata.json",
    D2_Fe_100 = "2D-SW-1_Fe_Detector2_100keV-mono_metadata.json",
    D2_Fe_200 = "2D-SW-1_Fe_Detector2_200kV-poly_metadata.json",
    D3_Al_50  = "2D-SW-1_Al_Detector3_50keV-mono_metadata.json",
    D3_Al_120 = "2D-SW-1_Al_Detector3_120kV-poly_metadata.json",
    D3_Fe_100 = "2D-SW-1_Fe_Detector3_100keV-mono_metadata.json",
    D3_Fe_200 = "2D-SW-1_Fe_Detector3_200kV-poly_metadata.json",
    D4_Al_50  = "2D-SW-1_Al_Detector4_50keV-mono_metadata.json",
    D4_Al_120 = "2D-SW-1_Al_Detector4_120kV-poly_metadata.json",
    D4_Fe_100 = "2D-SW-1_Fe_Detector4_100keV-mono_metadata.json",
    D4_Fe_200 = "2D-SW-1_Fe_Detector4_200kV-poly_metadata.json"
)
```

If only a sub-sample should be tested, it is not necessary to provide all simulation results at once. However, if a scintillator scenario misses its corresponding ideal detector scenario, it cannot be evaluated.

For the evaluation, the ratio of scintillator gray value to ideal detector gray value is calculated for each step of the spherical wedge, as well as for a free beam region. The grey values are measured as the pixel mean within the regions of interest (ROI) given in the table below and illustrated in the following image, defined as rectangles with an upper left corner coordinate (x<sub>0</sub>, y<sub>0</sub>) in the pixel coordinate system which is included in the ROI, and a lower right corner coordinate (x<sub>1</sub>, y<sub>1</sub>) which is not included anymore in the ROI.

| Step      | x<sub>0</sub> | y<sub>0</sub> | x<sub>1</sub> | y<sub>1</sub> |
| :-------- | :------------ | :------------ | :------------ | :------------ |
| 1         | 510           | 854           | 530           | 875           |
| 2         | 510           | 763           | 530           | 784           |
| 3         | 510           | 672           | 530           | 693           |
| 4         | 510           | 581           | 530           | 602           |
| 5         | 510           | 490           | 530           | 511           |
| 6         | 510           | 399           | 530           | 420           |
| 7         | 510           | 308           | 530           | 329           |
| 8         | 510           | 217           | 530           | 238           |
| 9         | 510           | 126           | 530           | 147           |
| free beam | 510           | 35            | 530           | 56            |

![Regions of interest on the step wedge](../pictures/wedge_ROIs.png "Regions of interest on the step wedge")

The gray value ratios are compared to reference ratios obtained from particle transport Monte-Carlo simulations. The gray value results of such simulations are subject to stochastic noise. For a given mean Monte-Carlo gray value intensity &mu;<sub>scint</sub> in a scintillator scenario, and a corresponding gray value intensity for the ideal detector &mu;<sub>ideal</sub>, Gaussian error propagation gives the error *u* of the Monte-Carlo gray value ratio, given the standard deviations \(\sigma_\text{scint}\) and \(\sigma_\text{ideal}\) of the grey values inside the ROI under consideration:

$$ u\left(\frac{\mu_\text{scint}}{\mu_\text{ideal}}\right) = \sqrt{\frac{\sigma_\text{scint}^2}{\mu_\text{ideal}^2} + \frac{\sigma_\text{ideal}^2 \mu_\text{scint}^2}{\mu_\text{ideal}^4}} $$

During the Monte-Carlo simulations, the effects of primary radiation and scatter radiation can be separated. Accordingly, the evaluation provides two sets of results: one where only primary radiation is considered (file names contain `primary`), and one where both scattered radiation and primary radiation are considered (file names contain `scatter`). Depending on which kind of software is tested, the correct comparison set should be chosen from the two options.

An evaluation run will produce plots such as the example shown below. These plots display the measured gray value ratios and show them in relation to the results from the Monte-Carlo simulations with their corresponding uncertainty *u*.

The following example evaluation result plot displays the gray value ratios between detector 2 (200 &mu;m CsI) and detector 1 (the ideal detector). Small crosses represent the gray value ratios calculated from the Monte-Carlo simulations. Their error bars cover the ratio's uncertainty *u* in both directions. Solid shapes without error bars show the measured ratios from the projections of the simulation software.

![2D-SW-1 example evaluation result](../pictures/2D-SW-1_D2.png "2D-SW-1")

Additionally, the following data files are created:


+ `2D-SW-1_all_GV_means.txt`

    lists all the mean grey values in one file.

+ `2D-SW-1_all_GV_ratios_to_D1.txt`

    lists the grey value ratios for detectors D2, D3 and D4 to the scenarios from detector D1. It is self-evident that the grey value ratios of detector D1 to itself are all 1, but they are listed as well for clarity.

+ `2D-SW-1_all_MCprimary_reference_ratios.txt`

    lists all the grey value ratios and their uncertainties *u* calculated from the Monte-Carlo simulations for the case of only primary radiation simulated (i.e., no scatter radiation taken into account).

+ `2D-SW-1_all_MCscatter_reference_ratios.txt`

    lists all the grey value ratios and their uncertainties *u* calculated from the Monte-Carlo simulations for the case of both primary and scatter radiation simulated.

+ For each sub test, a file such as `2D-SW-1_D1_Al_50_grey_values.txt` is generated which contains the ROI information and mean grey values.