This test scenario consists of two subtests. For a free-beam, projection images with two different noise levels are to be simulated: a signal-to-noise ratio (SNR) of 100, and a second one with an SNR of 250, both of which refer to the point of maximum intensity.

Once simulated, the metadata file for each subtest that shall be evaluated must be passed to the toolbox, identified by the correct argument keyword `SNR100` or `SNR250`. You don't have to run both test evaluations at once, but it is possible.

```python
from ctsimu.toolbox import Toolbox
Toolbox("2D-FB-1",
  SNR100="2D-FB-1_Detektor1_SNR-100_metadata.json",
  SNR250="2D-FB-1_Detektor2_SNR-250_metadata.json"
)
```

The detector is small and very far away from the source, leading to a very homogeneous illumination. The simulation software must provide a noise-free projection image to be used for the flat-field correction of the noisy projection image. Both files must be correctly referenced in the metadata file that is passed to the toolbox. The SNR is then evaluated from the complete pixel ensemble of the flat-field corrected image:

$$ \text{SNR} = \frac{\left< I \right>}{\sigma}, $$

with \(\left< I \right>\) being the mean grey value and \(\sigma\) the root mean square grey value deviation (RMSD) of the pixel ensemble.