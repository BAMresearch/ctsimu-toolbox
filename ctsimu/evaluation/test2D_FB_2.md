This test scenario checks if the intensity profile is correctly rendered under the assumption of an isotropic point source. The inverse square law and the law of radiation incidence must be obeyed:

$$I\sim \cos(\alpha)/r^2$$

The test compares the central horizontal row of pixels at index `y=250` with the analytical intensity profile.

To run this test evaluation, you can simply pass the name of the test and the metadata file which describes the projection that you want to evaluate to the `ctsimu.toolbox.Toolbox`:

```python
from ctsimu.toolbox import Toolbox
Toolbox("2D-FB-2", "2D-FB-2_metadata.json")
```