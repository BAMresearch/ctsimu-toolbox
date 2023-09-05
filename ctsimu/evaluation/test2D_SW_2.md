In this test scenario, spectral filtering is tested. Two spherical wedges with 9 steps of different thicknesses are simulated: one made of aluminum, the other one made of iron. For both wedges, we define two scenarios each: one at an X-ray spectrum of 200 kV without any additional filters (apart from the tube window), the other one with an additional copper filter (2 mm) by which the spectrum shall be filtered. This will result in grey value ratios between the wedge steps which are characteristic to each scenario. Those ratios are compared to the ratios obtained from particle-transport Monte-Carlo simulations of the corresponding scenarios.

When running the test evaluation, please take care to use the correct keywords for each sub-test, as shown in the code example below. It is not necessary to run all tests at once, but possible.

```python
from ctsimu.toolbox import Toolbox
Toolbox("2D-SW-2",
  Al   = "2D-SW-2_Al_200kV_noFilter_metadata.json",
  AlCu = "2D-SW-2_Al_200kV_Filter2mmCu_metadata.json",
  Fe   = "2D-SW-2_Fe_200kV_noFilter_metadata.json",
  FeCu = "2D-SW-2_Fe_200kV_Filter2mmCu_metadata.json"
)
```

The regions of interest (ROI) where the gray values are measured are the same as shown for scenario `ctsimu.evaluation.test2D_SW_1`. The errors of the gray value ratios between the steps of the wedges in the Monte-Carlo simulations are obtained by Gaussian error propagation in analogy to the procedure shown for 2D-SW-1.

The result files list the measured grey values and their corresponding ROIs, as well as the measured grey value ratios between subsequent steps. For the ratios, it also lists the Monte-Carlo reference values and their uncertainties. Just like for scenario 2D-SW-1, values are given for simulations considering only `primary` radiation (without scatter radiation), and simulations which consider both (labelled `scatter` in the result files). Plots such as the example shown below are produced for each sub-test.

![2D-SW-2 example evaluation result](../pictures/2D-SW-2_FeCu.png "2D-SW-2")

The example evaluation result displays the gray value ratios between each pair of neighboring steps of the iron wedge. Small crosses represent the grey value ratios calculated from the Monte-Carlo simulations. Their error bars cover the ratio's uncertainty *u* in both directions. Solid circles without error bars show the measured ratios from the projections of the simulation software.