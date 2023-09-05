Any spatially extended focal spot leads to a spatially extended point spread function (PSF) on the detector and therefore to image blurring. In this test scenario, we check if Gaussian spot intensity profiles are correctly simulated by the software to lead to the expected image blurring, represented by certain expected MTF20 values that we calculate as described below.

This scenario defines two sub-tests: one with a Gaussian spot intensity profile with a standard deviation of &sigma;<sub>1</sub> = 8 µm (Spot 1), and a second one with &sigma;<sub>2</sub> = 24 µm (Spot 2). The following code example shows how to identify each subtest correctly when using the toolbox to run the evaluation.

```python
from ctsimu.toolbox import Toolbox
Toolbox("2D-WE-1",
    spot1 = "2D-WE-1_Spot1_metadata.json",
    spot2 = "2D-WE-1_Spot2_metadata.json")
```

Our evaluation approach is to measure the edge spread function (ESF) across the edge of an ideal absorber (the same edge as used for `ctsimu.evaluation.test2D_WE_2`). In subsequent steps, we will calculate the line spread function (LSF) to evaluate the spot size and the modulation transfer function (MTF) to determine where the system's frequency transmission drops to 20% (MTF20 value). [Rossmann et al., 1969], [Jones et al., 1958].

[Rossmann et al., 1969]: https://doi.org/10.1148/93.2.257
[Jones et al., 1958]: https://doi.org/10.1364/JOSA.48.000934

The edge spread function (ESF) is simply the line profile across the edge. The center point *c* of the line profile is placed exactly in the center of the vertical part of the edge (see picture below). The starting point p<sub>0</sub> and the end point p<sub>1</sub> of the profile result from the line profile parameters: we choose a profile length of *l* = 100.1 px and a bin size (resolution) of *r* = 0.1 px. This results in a profile of *N*=1001 bins along the line across the edge. We choose a profile region that is *w* = 200 px wide.

![Line profile region of interest](../pictures/edge_lineprofile.png "Line profile region of interest")

Due to the pixel rasterization of the projection and the finite sampling of the line profile, we will lose frequency information even in the ideal case of a point source. We therefore use the ideal edge image from `ctsimu.evaluation.test2D_WE_2` to determine the sampled LSF and MTF under the ideal assumption of a point source. This will give us the best LSF and MTF that are achievable with the given detector, edge tilt angle, and the chosen evaluation method, which means they represent our scenario's fundamental spread and transfer functions. These results are shown in the following picture.

![ESF_LSF_MTF](../pictures/sampled_pointsource_ESF_LSF_MTF.png "Fundamental edge spread function (ESF), normalized fundamental line spread function (LSF) and normalized fundamental modulation transfer function (MTF) from the ideal edge image assuming a point source.")

A convolution of the fundamental ESF and LSF with the profile of a Gaussian spot will lead to the expected ideal ESF and LSF at finite spot sizes. Generally, this convolution has to be done using the point spread function (PSF), which, apart from amplitude normalization, is identical to the LSF when an isotropic Gaussian distribution of the spot's intensity is assumed. ([Gopala Rao et al., 1967])

[Gopala Rao et al., 1967]: https://doi.org/10.1364/JOSA.57.001159

To calculate the analytical sampled LSF of an ideal non-rotated edge under a Gaussian spot with a standard deviation &sigma;, we integrate the spot's intensity function along the full width of each bin of the line profile, assuming that the distribution is centered at \(s_\text{c}=\frac{N}{2}\) (with the number of bins *N*):
$$
  \text{LSF}_\text{Gaussian, j} = A\cdot \int_{j-\frac{N}{2}}^{j+1-\frac{N}{2}} \text{e}^{-s^2/(2\sigma^2)} \text{d} s = A \cdot \sqrt{\frac{2}{\pi}}\sigma \cdot \left[\text{erf}\left(\frac{j+1-\frac{N}{2}}{\sqrt{2}\sigma} \right) - \text{erf}\left(\frac{j-\frac{N}{2}}{\sqrt{2}\sigma} \right) \right].
$$
The parameter *A* is chosen such that the maximum of the LSF is normalized to a value of 1. The ideal expected LSF is then computed by the discrete convolution (assuming the number of bins *N* is an odd number)
$$
  \text{LSF}_{\text{ideal}, j} = (\text{LSF}_\text{fundamental} \ast \text{LSF}_\text{Gaussian})_j = \sum_{k=0}^{N-1} \text{LSF}_{\text{fundamental}, k} \cdot \text{LSF}_{\text{Gaussian}, j-k+\frac{N-1}{2}}.
$$
In the scenario, the edge is located exactly halfway between source and detector. Therefore, no magnification of the point spread function (PSF) or line spread function (LSF) is to be expected, and the Gaussian spot size can be determined directly from the measured LSF without the need for scaling. In theory, it would be possible to use a deconvolution to remove the fundamental LSF from the measured LSF of a given simulated image such that only the shape of the spot intensity profile remains. However, deconvolution is an ill-posed problem, and it is only possible if all of our ideal assumptions are met by the simulation software. To avoid instabilities and the introduction of any additional errors, we therefore determine the expected spot sizes by Gaussian fits to the ideal LSF, i.e., from the convolution of fundamental LSF and Gaussian LSF that results from the above equation. Using a least-squares fit, we determine the expected ideal \(\tilde{\sigma}\) that are listed in the following table.


| Spot size                      | Expected LSF std. dev.                | Expected MTF20 frequency             |
| :----------------------------- | :------------------------------------ | :----------------------------------- |
| &sigma;<sub>1</sub> = 8 &mu;m  | \(\tilde{\sigma}_1\) =  9.579(2)&mu;m | MTF20<sub>1</sub> = 29.852 cycles/mm |
| &sigma;<sub>2</sub> = 24 &mu;m | \(\tilde{\sigma}_2\) = 24.552(1)&mu;m | MTF20<sub>2</sub> = 11.770 cycles/mm |

The ideal MTF can be calculated either by multiplication of the fundamental MTF and the Gaussian MTF (as calculated from the Gaussian LSF),
$$
  \text{MTF}_\text{ideal} = \text{MTF}_\text{fundamental} \cdot \text{MTF}_\text{Gaussian},
$$
or by applying a Fourier transform to the ideal LSF. Apart from numerical artifacts, both results are identical; here, we use the second approach.

For each image to be measured, the MTF20 value is determined and compared to the MTF20 value of the ideal MTF. The MTF20 value is the frequency where the MTF drops to 20% of its maximum value. The expected values from our ideal MTFs are listed in the table above.