This test checks for correct positioning of sample and detector. An STL file that represents a thin foil with ten holes must be placed according to the correct bounding box position and orientation as described in the scenario. The STL file's native coordinates do not match the coordinates where it has to be placed. Additionally, the correct detector position and orientation must be set up in order to get the correct projection image. The evaluation procedure will calculate translation vectors, rotation angles and scale factors relative to the ideal pixel coordinates of the holes.

```python
from ctsimu.toolbox import Toolbox
Toolbox("2D-HS-1", "2D-HS-1_metadata.json")
```

The test tries to find the ten holes in the projection image. If it cannot identify exactly ten holes, the test will be aborted completely. The centers of the holes are determined by fitting a circle to each hole boundary ([Coope et al., 1993]) after application of an edge filter. Hole number 0 has a bigger diameter than the other holes and serves as a symmetry breaker and reference anchor from which all other holes are identified.

[Coope et al., 1993]: https://doi.org/10.1007/BF00939613

![Hole sheet standard](../pictures/hole_sheet.png "Hole sheet standard with holes labeled")