If a detector pixel is partly covered by an ideal absorber, its gray value should (approximately) scale inversely proportionally with the area fraction that is covered. This effect is tested here using a very thin (0.01 mm) edge that is made of a high-density material. Its upper right corner is placed at the center of the detector's central pixel, and then tilted around this point by 3&deg;, as illustrated in the figure below.

To run the evaluation for this test, simply pass its identifier and the name of your metadata file (for projection and flat field) to the toolbox:

```python
from ctsimu.toolbox import Toolbox
Toolbox("2D-WE-2", "2D-WE-2_metadata.json")
```

![Partial pixel coverage](../pictures/edge.png "A thin edge partly covers pixels on the detector.")

To calculate the ideal projection image of an ideal edge, we use the spherical trigonometry approach (see internal toolbox documentation: "How it works under the hood") to obtain an analytical value for each pixel's intensity. For pixels which are fully exposed or fully covered, this is straightforward. Pixels that are partly covered must be separated into an exposed region and a covered region. To treat the covered regions, we regard each pixel as a square polygon in the detector plane. Each of these polygons is then clipped using the shadow of the edge as a clipping polygon. By employing the [Sutherland-Hodgman polygon clipping] algorithm we obtain a clipped polygon for which we can calculate the solid angle in the same way as we do for a full pixel. The edges of these clipped polygons are illustrated by dotted red lines in the figure above.

[Sutherland-Hodgman polygon clipping]: https://doi.org/10.1145/360767.360802

The intensity of a pixel is therefore proportional to the solid angle of the uncovered area of the full pixel polygon:

$$ I(\text{pixel}) \sim \Omega(\text{full pixel}) - \Omega(\text{clipped pixel}). $$

The toolbox regards any pixel with a gray value above zero and below 60000 as *partly covered*. For the analytical image, it identifies the set *A* of partly covered pixels, and for the provided projection image to be evaluated, it identifies a set *M* of partly covered pixels.

The following results are then calculated.

+ A difference image, showing the gray value difference \(\Delta = \text{GV}_\text{measured}-\text{GV}_\text{analytical}\) for each pixel.
+ The number of partly covered pixels in the analytical image, \(\left|A\right|=1053\), and in the projection image to be evaluated, \(\left|M\right|\).
+ The number of pixels that are partly covered in both images: \(\left|A \cap M\right|\), i.e., the number of pixels in the intersection of both sets. This number should match \(\left|M\right|\) if the edge is positioned correctly. If \(\left|M\right|\) is higher, this means that partly covered pixels have been identified in the image to be evaluated which should be either fully exposed or fully covered.
+ The ratio *r* of the number of pixels that are covered in both vs. the number of pixels covered in the analytical image:
    $$
    r = \frac{\left|A \cap M\right|}{\left|A\right|}
    $$
    Ideally, this ratio is 1 (i.e., 100%). However, it only makes a statement about the extent to which *M* covers the ideal set *A*, and does not take into account any pixels that might be incorrectly partly covered outside of the region of analytically partly-covered pixels.
+ The root mean square difference between the analytical and the measured image for all pixels in the analytically ideal set *A*, i.e., for all pixels that are partly covered in the analytical edge image:
    $$
    \text{RMSD} = \sqrt{\frac{1}{\left|A\right|} \sum_{p \in A} \left[ I_\text{measured}(p) - I_\text{analytical}(p) \right]^2 }.
    $$
+ A list of the coordinates of *partly covered* pixels in the analytical edge image (i.e., all pixels from set \(\left|A\right|\), their gray value in the analytical and in the measured image, and their respective difference.
