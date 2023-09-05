The term *Boolean models* refers to defined, reproducible behavior in regions where models overlap. In these regions, an object of higher priority is supposed to completely replace any other object of lower priority. This test employs the step wedge known from the former two scenarios, but splits it into a left-hand piece that shows straight edges in the projection image and a right-hand piece that is beveled (see the following picture). When positioned correctly, both pieces touch at the central three steps (contiguous surfaces). Between the upper three steps, there is a small gap of decreasing width, and for the lower three steps, both pieces overlap. The right-hand piece is supposed to completely replace the left-hand piece in those overlap regions

![Boolean models](../pictures/boolean_models.png "Boolean models")

The picture shows example projection images for the two tests for 2D-SW-3. In each of the two sub-scenarios, the beveled right-hand piece is supposed to replace the left-hand piece in the region of the lower three steps where both models overlap. For the upper three steps, there remains a gap of decreasing width and at the central three steps, both surfaces touch.

The test considers two sub-scenarios. In the first scenario, the two pieces are both assumed to be of aluminum. In the projection image, there should be no detectable gray value change at the object boundaries for the lower six steps. For the second scenario, the left-hand piece is assumed to be of aluminum and the right-hand piece of titanium. There must be a detectable gray value change at pixels where the expected material boundaries are located. In all other areas, the gray values within each step should not change.

Both scenarios can be evaluated by the toolbox at once, but single evaluations of each sub-scenario are possible as well. Use the keywords `Al_Al` and `Al_Ti` as shown in the following example to correctly identify each sub-scenario.

```python
from ctsimu.toolbox import Toolbox
Toolbox("2D-SW-3",
    Al_Al = "2D-SW-3_Al-Al_metadata.json",
    Al_Ti = "2D-SW-3_Al-Ti_metadata.json"
)
```

To detect locations of gray value changes and possible anomalies, the first step in the evaluation process is to calculate the gray value derivatives for each horizontal line of pixels. The scenario specifies a 1000x1000 pixel projection image. We assume a pixel coordinate system where *x* denotes the index within a row and *y* denotes the index within a column. The upper left pixel has index (0, 0), the lower right pixel is located at index (999, 999). In the derivative image, a pixel's gray value \(p^\prime\) is calculated as the difference of the original gray value \(p\) of its successor at position *x+1* and its own gray value:

$$ p^\prime (x, y) = p(x+1, y) - p(x, y) ~~~~ \forall~x \in [0, 998],~y \in [0, 999]. $$

The last pixel column at \(x=999\) is omitted as it has no successor. We therefore get a derivative image with a size of 999x1000 pixels. Numerical limitations may lead to minor deviations in regions where constant gray values are expected. We ignore any gray value changes below 0.04% of the free beam intensity, which corresponds to 24 gray values for the required free beam gray value of 60000. This leads to cleaned gray value derivatives

$$
    p^\prime_{\text{cleaned}}(x, y) = \begin{cases}
        0 & p^\prime(x, y) < 24\,\text{GV} \\
        p^\prime(x, y) & p^\prime(x, y) \geqslant 24\,\text{GV}
    \end{cases}.
$$

The scenario requires monochromatic radiation with a photon energy of 80 keV. The specified densities of \(\rho_\text{Al}=2.6989\,\text{g}/\text{cm}^3\) and \(\rho_\text{Ti}=4.506\,\text{g}/\text{cm}^3\) and the mass attenuation coefficients \(\mu_\text{Al}/\rho=0.2018\,\text{cm}^2/\text{g}\) and \(\mu_\text{Ti}/\rho=0.4052\,\text{cm}^2/\text{g}\) ([Hubbell et al.]) lead to the specific attenuation coefficients \(\mu_\text{Al}=0.0544638/\text{mm}\) and \(\mu_\text{Ti}=0.18234/\text{mm}\) for the materials used in this scenario. The Beer-Lambert law of attenuation, the thickness of each step and the requirement for a free beam intensity of 60000 gray values lead to the expected gray values (after flat-field correction) listed in the following table.

[Hubbell et al.]: https://doi.org/10.18434/T4D01F


| Step                | 1     | 2     | 3     | 4     | 5     | 6     | 7     | 8     | 9     |
| :------------------ | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: |
| Thickness in mm     | 41.00 | 35.89 | 30.79 | 25.68 | 20.58 | 15.47 | 10.36 | 5.26  |  0.15 |
| Al gray value       | 6432  | 8496  | 11217 | 14816 | 19560 | 25837 | 34127 | 45054 | 59512 |
| Ti gray value       | 34    | 86    | 219   | 555   | 1407  | 3573  | 9073  | 22994 | 58381 |
| &Delta;             | 6398  | 8410  | 10998 | 14261 | 18153 | 22264 | 25054 | 22060 | 1131  |
| y<sub>0</sub> in px | 819   | 728   | 637   | 546   | 455   | 364   | 273   | 182   | 91    |
| y<sub>1</sub> in px | 908   | 817   | 726   | 635   | 544   | 453   | 362   | 271   | 180   |


The table lists the total material thickness of each step of the spherical wedge that a ray from the source encounters. The expected analytical gray values from the law of radiation attenuation are listed, as well as the difference &Delta; between aluminum and titanium. In the projection image, each step is fully present in their respective region \(y\in[y_0, y_1]\). The rows in between are *forbidden* lines where transitions between steps take place.}

Each step of the wedge occupies 1/11 of the height of the projection image. The locations \(y_\text{step}\) of their boundaries are listed on the right side of the graph. \(\lfloor y_\text{step} \rfloor\) gives the floor function value, i.e., the pixel index of the *forbidden rows* where the transitions from one step to another reside. The maximum widths of the material gap and the overlap region are shown.

![Boolean models](../pictures/wedge_step-y-locations.png "Boolean models")

The algorithm now treats each horizontal row of the image of gray value derivatives separately. For each row, the following procedure is applied.

1. The pixels with indices \(x\in[80, 99]\) and \(x\in[870, 919]\) are set to zero. This ignores any gray value transitions and pixel anomalies in the left and right border regions where gray values transition from free beam intensities to material intensities (i.e., where the wedge begins and ends).

2. **Gap region.** For each pixel row, the width of the air gap is calculated. For continuous coordinates \(x\) and \(y\) in the pixel coordinate system, the width \(w_\text{gap}\) (in pixels) is given by the linear function
    $$w_{\text{gap}}(y) = -\frac{5\,\text{px}}{3h_{\text{step}}} \left( y - 4h_{\text{step}} \right),$$
    with the height of a step \(h_\text{step} = \frac{1000\,\text{px}}{11}\). See the above picture for reference. If we only allow discrete integer coordinates for the pixels, each row \(r\) therefore has an upper and a lower gap width:
    $$w_{\text{gap,upper}}(r) = -\frac{5\,\text{px}}{3h_{\text{step}}} \left( r - 4h_{\text{step}} \right),$$
    $$w_{\text{gap,lower}}(r) = -\frac{5\,\text{px}}{3h_{\text{step}}} \left( (r+1) - 4h_{\text{step}} \right).$$
    The left boundary of the air gap is marked by the pixel at \(x_{\text{gap,left}}=500\), the first one to contain free-beam contributions in the intensity. The location of the right boundary depends on the gap width. We can calculate two pixel indices for the interval in which the boundary should be found, tolerances included:
    $$x_{\text{gap,right,0}} = x_{\text{gap,left}} + \lfloor w_{\text{gap,lower}} \rfloor - 1,$$
    $$x_{\text{gap,right,1}} = x_{\text{gap,left}} + \lfloor w_{\text{gap,upper}} \rfloor.$$
    Now we need to consider the following cases, ordered from highest to lowest priority.

    - If \(w_{\text{gap,lower}}\leq 0\) we are not in a region where the gap is found (step number is lower than 7). In this case, we can continue with step~3 of the algorithm.

    - If \(0 < w_\text{gap,lower}\leq 1\) we assume that the air gap is too small to be fully resolved in between the two wedges. The sum of gray value changes across the complete gap region is calculated. It must match (i.e., be within 5% of) the expected gray value transition from the left material (Al) to the right material (Al or Ti). If this is the case, we set the gray value differences in the area of the gap to zero for all \(x\in\left[ 498, x_{\text{gap,right,1}} \right]\). If the current row \(y\) is a *forbidden row*, we ignore the gray value differences and set the pixels to zero as well. If both conditions are not fulfilled, we set the value of the pixel at \(x=499\) to the expected gray value difference and all other pixels in the considered range to zero. This puts the whole weight of the anomaly onto one pixel and makes counting anomalies easier afterwards.

    + If \(w_{\text{gap,lower}} > 1\) it should be possible to resolve the air gap (i.e., there should be at least one pixel that contains the free beam intensity in the original projection).

        - At first, we treat the boundary between the material of the left wedge and the environment. In the difference image, the sum over the boundary pixels at \(x\in\left[ 498, 499 \right]\) should match (within 5%) the expected gray value change from the Al step to the free beam intensity. If \(y\) is a *forbidden row*, we only check if the slope of the profile points upwards towards higher gray values, i.e., we check for a positive gray value difference for both pixels under consideration. If any of these two conditions is fulfilled, we set the two boundary pixels to zero in the difference image. Otherwise we set the value of the pixel at \(x=499\) to the expected gray value difference, and all other pixels in the considered range to zero.

        - We repeat this procedure for the right boundary, from environment to the right wedge, for the pixels at \(x\in \left[x_{\text{gap,right,0}}, x_{\text{gap,right,1}} \right]\). It is checked whether there is the expected (within 5%) drop in gray values to lower intensities. For *forbidden rows*, it is only necessary that there is a negative difference (decreasing slope). Otherwise, an anomaly is introduced to the difference image, as described above.

3. **Contact region.** The pixel rows that belong to steps 4, 5 and 6 are considered. Only the three pixels at the material intersection are treated: \(x\in\left[ 498, 500 \right]\). For the Al/Al scenario, nothing is done here. There is no expected gray value difference, and any existing gray value anomaly will persist for the final evaluation. For the Al/Ti scenario, the sum of gray value differences of these three pixels must match (within 5%) the expected gray value difference between Al and Ti. If the current row \(y\) is a *forbidden row*, the gray value differences for these three pixels must be zero or negative, i.e. the gray value must not rise at the transition from Al to Ti. If any of these two conditions is fulfilled, we set the value of these three pixels to zero in the difference image. Otherwise, we set the value of the pixel at \(x=498\) to the expected gray value difference, and all other pixels in the considered range to zero.

4. **Displacement region.** The pixel rows for steps 1, 2 and 3 are treated very similarly. For the Al/Al scenario, no gray value difference is expected and nothing is done here. For the Al/Ti scenario, the leftmost and rightmost pixel coordinate \(x\) is calculated that should be considered when checking the gray value changes:
    $$x_{\text{boundary,0}}(r) = 498-\frac{10\,\text{px}}{3h_{\text{step}}} \left( r - 7h_{\text{step}} \right),$$
    $$x_{\text{boundary,1}}(r) = 500-\frac{10\,\text{px}}{3h_{\text{step}}} \left( (r+1) - 7h_{\text{step}} \right).$$
    The sum of gray value changes for all pixels in \(x\in\left[x_{\text{boundary,0}}, x_{\text{boundary,1}} \right]\) must match (within 5%) the expected gray value difference between Al and Ti. If the current row \(y\) is a *forbidden row*, the gray value differences for these three pixels must be zero or negative, i.e., the gray value must not rise at the transition from Al to Ti. If any of these two conditions is fulfilled, we set the value of these three pixels to zero in the difference image. Otherwise, we set the value of the pixel at \(x_{\text{boundary,1}}\) to the expected gray value difference, and all other pixels in the considered range to zero.

After this algorithm has been applied to each pixel row, all pixels in the difference image should have a value of zero. Pixels with other values are counted as anomalies. The mean gray value difference of all pixel anomalies is calculated.