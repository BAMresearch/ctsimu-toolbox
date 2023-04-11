# Coordinate Systems

The geometry subpackage provides a `CoordinateSystem` class that lets you create and manipulate objects in a virtual CT scene. Such an oject has a position (`center`) and three basis vectors (`u`, `v`, `w`). These basis vectors are assumed to be orthogonal, but they do not have to be unit vectors. The `center` acts as the pivot point for rotations.

Example for creating and manipulating an object:

```python
.. include:: ../examples/geometry/01_coordinate_systems.py
```

# Full CT Geometry

For a full CT, we need an X-ray source, a stage for the specimens, and a detector. The `Geometry` class bundles three coordinate systems (one for each of those components), and additional information about the detector (using the `Detector` class, an extension of a regular `CoordinateSystem`). The following figure shows their standard orientations when a CT geometry is initialized.

![Standard coordinate system](pictures/geometry.png "Standard coordinate system")

The orientation of the coordinate system and all components can be changed by rotations or by manually setting the object basis vectors. However, it is important to keep the following conventions.

## Detector convention

* The detector's `u` vector is its row vector.
* The detector's `v` vector is its column vector.
* The detector's `w` vector has no special meaning. It is a planar normal that must be chosen such that the detector's coordinate system remains right-handed.

## Stage convention

* The stage's `w` vector is its axis of CT rotation.
* The stage's `center` typically refers to the center of the reconstructed volume (possibly the center of the specimen) and is *not* meant to describe the location of the turntable (which would normally be at a lower position).

## Source convention

There is currently no restriction on the source coordinate system. We usually assume its `w` axis to be the direction of the principal ray, but this is not a necessity.

In the following examples, the source will be located at the origin `(0, 0, 0)` of the world coordinate system, whereas stage and detector are placed in positive *x* direction (see figure above).

## Example Setup

In the following example, we set up a standard CT geometry.

```python
.. include:: ../examples/geometry/02_simple_CT_geometry.py
```

# Reference Frames

Implicitly, each coordinate system has a reference coordinate system (its *reference frame*) in which its `center` and `u`, `v`, `w` basis vectors are located and described. Typically, we assume that this is a right-handed standard coordinate system. It does not necessarily have to be the world coordinate system. For example, you might want to attach a specimen to the sample stage by implicitly making the stage coordinate system its reference coordinate system. You, the programmer, have to know the reference coordinate system of your objects, as this information is not explicitly stored by the toolbox.

Any new `CoordinateSystem` object is initialized to be a right-handed standard coordinate system with its center at `(0, 0, 0)`:

```python
from ctsimu.geometry import *

myWorld = CoordinateSystem()

print("My World:")
print(myWorld)

"""
My World:
Center: ( 0.0000000,  0.0000000,  0.0000000)
u:      ( 1.0000000,  0.0000000,  0.0000000)
v:      ( 0.0000000,  1.0000000,  0.0000000)
w:      ( 0.0000000,  0.0000000,  1.0000000)
"""
```

You can change the reference frame of a `CoordinateSystem`. In the following example, we set up a CT geometry with a stage that is tilted by 2°. We place a specimen object in the stage coordinate system and move it "upwards" by 5 mm along the (now tilted) axis of rotation. Afterwards, we change the specimen's reference frame to see where it is actually located in the world coordinate system. Refer to the image of the standard orientations above to see what is going on with the coordinate systems.

```python
.. include:: ../examples/geometry/03_reference_frames.py
```

**Note:** When changing reference frames, the original (`csFrom`) and the target reference frame (`csTo`) must both have the same common reference frame for themselves. In the example above, we change the reference frame from the stage coordinate system to the world coordinate system. Both of them have the same reference frame: the world coordinate system (which is special, because it is also a reference for itself).


# Projection Matrices

A projection matrix maps a 3D point coordinate `(x, y, z)` from the stage coordinate system to a 2D point coordinate `(u, v)` in the detector coordinate system. They are used by some reconstruction softwares to describe arbitrary scan trajectories. For such a reconstruction, we need one projection matrix for each projection image.

![Euclidean Mapping](pictures/pmatrix_mapping_euclidean.png "Mapping a coordinate from the stage coordinate system to a coordinate in the detector coordinate system")


## Mathematical Background

We operate in [homogeneous coordinates](https://en.wikipedia.org/wiki/Homogeneous_coordinates) because we are in a projective geometry and this concept allows us to describe translations in space by a matrix. Homogeneous coordinates describe rays in a space and are therefore only defined up to a scale factor, which is carried as an additional coordinate. This way, a Euclidean 3-vector turns into a homogeneous 4-vector. Our mapping becomes:

![Homogeneous Mapping](pictures/pmatrix_mapping_homogeneous.png "Mapping in homogeneous coordinates")

The following picture illustrates a 1D projective geometry. The *h* axis is our scale factor for the homogeneous coordinates. In this geometry, all points on a ray are equivalent. The points at *h*=1 are the normalized homogeneous coordinates.

![Projective geometry for a 1D Euclidean space](pictures/hom_coords.png "Projective geometry for a 1D Euclidean space")

We can use this concept to describe translations in space with a matrix multiplication:

![Translation Matrix](pictures/pmatrix_translation.png "Translation Matrix")

We now consider two coordinate systems: the stage coordinate system (our origin) and the detector coordinate system (our target). The world coordinate system is not important in this context. In the following picture, 3D (Euclidean) coordinates are expressed in terms of the stage coordinate system, and 2D (Euclidean) coordinates are expressed in terms of the detector coordinate system (ignoring its *w* axis).

![Projective system](pictures/geometry_stage_projection.png "Projective system")

To calculate a projection matrix for the current geometry, we have to consider five subsequent transformations. Each transformation is expressed by a matrix. The final projection matrix is then the product of these five transformation matrices.

1. We shift the origin of the coordinate system from the stage to the source, which is our center of projection.

    ![Shifting the origin to the source](pictures/pmatrix_F.png "Shifting the origin to the source")

2. We perform a basis transformation to express the 3D coordinates in terms of the axes of the detector coordinate system. The origin remains at the source. You can also think of this transformation as a rotation into the detector coordinate system.

    ![Rotation into detector coordinate system](pictures/pmatrix_R.png "Rotation into detector coordinate system")

    All basis vectors in this matrix are assumed to be unit vectors.

    This matrix takes care of any stage or detector tilts.

    After this transformation, the third (*"z"*) coordinate in a vector now refers to its position on the detector normal (its *w* axis). Therefore, this third coordinate now contains something similar to what we would normally call the SOD (source-object distance) of that point. The fourth coordinate of our homogeneous vector has not been scaled so far (α=1), which means we have not left the projective plane which we call home (our real world). This is important to keep in mind for the next step.

3. We use a matrix that reduces the dimension of our vector by one (from a homogeneous 4-vector to a homogeneous 3-vector). This step is sometimes called the actual *projection*.

    ![Projection reduces dimension](pictures/pmatrix_D.png "Projection reduces dimension")

    In the previous step, the third component of the 4-vector used to be something similar to the SOD. This has now become the *scale component* β of our homogeneous 3-vector (because a multiplication with this matrix throws away the fourth vector component, which has still been α=1). This means we are now in a projective plane β=SOD, away from the detector plane of our home world (which would be at β=1).

    This problem is solved in the end by a simple renormalization of the matrix. Stay tuned!

4. We take care of the magnification and any additional scaling.

    ![Scaling](pictures/pmatrix_S.png "Scaling")

    The SDD (source-detector distance) in this case means the length of the principal ray from source to detector (i.e., the ray that is parallel to the detector normal *w* and orthogonally hits the detector plane).

    This matrix simply scales any image at the projective plane β=1 such that its *u* and *v* component will obey the magnification by the SDD (source-detector distance). Sometimes, the stage coordinate system is expressed in a different unit than the detector coordinate system (e.g. mm vs. px). In this case, we can introduce scale factors s<sub>u</sub> and s<sub>v</sub> that take care of further scaling, e.g. to handle the pixel size.

    Note that the final renormalization will turn out to be a division by the SOD (as mentioned in the previous step). This will convert the SDD-factors of this matrix into the actual magnification: M=SDD/SOD. We do not incorporate this here because the SOD as a parameter is not well-defined and might lead to confusion in a non-standard geometry.

5. The origin of the detector coordinate system might not be where the principal ray hits the detector (i.e., the center of projection projected onto the detector). We need to take care of this additional shift:

    ![Translation on detector](pictures/pmatrix_T.png "Translation on detector")

The final projection matrix is a 3×4 matrix that results from a multiplication of these five matrices and a renormalization by the lower-right component (p<sub>23</sub>) to get back to the projective plane of our home world (see step 3).

![Final projection matrix](pictures/pmatrix_P.png "Final projection matrix")


## Generating Projection Matrices

You can call the function `Geometry.projectionMatrix()` to get a projection matrix for the geometry's current configuration.

```python
.. include:: ../examples/geometry/04_projection_matrix.py
```

### openCT & CERA

The toolbox provides two pre-configured modes to calculate projection matrices for openCT (which can be used in VGSTUDIO MAX) and for SIEMENS CERA. Each software needs slightly different projection matrices, because they define their detector coordinate system in different ways. See the next section about the image and volume coordinate system for details.

In the following example, we calculate a projection matrix for each software by defining the `mode` when calling the `Geometry.projectionMatrix()` function.

```python
.. include:: ../examples/geometry/05_projection_matrix_modes.py
```

openCT's definition of the detector coordinate system matches our definition. In this case, we get the same projection matrix as without the `mode` parameter.


### Image & Volume Coordinate Systems

Depending on the reconstruction software, the **image coordinate system** of the projection image does not have to match our standard detector coordinate system. Also, the **volume coordinate system** of the reconstructed volume does not have to match the stage coordinate system.

In the following three examples, we will show how to use the parameter `imageCS` and `volumeCS` to define our own image and volume coordinate systems.

![Image and volume coordinate system](pictures/image_stage_cs.png "Image and volume coordinate system")

**Note:** The image coordinate system is expressed in terms of the detector coordinate system (its reference coordinate system). Similarly, the volume coordinate system is expressed in terms of the stage coordinate system. To set the scale factor for the image or volume coordinate system, we set the lengths of their basis vectors to the correct conversion factor (e.g., the pixel size in mm/px or the voxel size in mm/voxel).

#### Example 1: CERA

Even though we have a pre-defined more for CERA, we will use its image coordinate system (illustrated above) to show how to set up an image coordinate system for CERA manually.

CERA's volume coordinate system matches our stage coordinate system, so we won't have to create our own volume coordinate system.

CERA's image coordinate system has its origin in the center of the lower left pixel of the detector. This means we have to move its origin by half the detector's physical width to the left and half the detector's physical height downwards from the origin of the detector coordinate system, and then back by half a physical pixel size (the detector pitch). We can use the attributes `physWidth` and `physHeight`, which are automatically calculated when calling `setSize()`.

```python
.. include:: ../examples/geometry/06_projection_matrix_cera.py
```

#### Example 2: openCT

In the case of openCT, the image coordinate system matches our detector coordinate system. Also, the volume coordinate system matches our definition of the stage coordinate system. In this case, we can simply call the `projectionMatrix()` function without any parameters.

However, for the sake of clarity, we create a standard coordinate system for both the image and volume coordinate system. Because they are expressed in terms of their respective reference coordinate systems (detector and stage), they exactly represent their respective reference CS as seen from the outside world.

```python
.. include:: ../examples/geometry/07_projection_matrix_openCT.py
```

#### Example 3

For the third example (see illustration), we will have to move the origin of the image coordinate system to the upper left corner of the detector and scale its basis vectors by the pixel size because its units are pixels.

The volume coordinate system has its origin at the front lower right corner of the reconstruction volume. Because it is no longer at the stage's center, we will actually have to define the volume's physical size in order to correctly calculate the corner coordinate in terms of the stage coordinate system.

We also scale the basis vectors of the volume coordinate system by the voxel size, because we assume that our reconstruction software expresses its volume coordinates in voxel units instead of world units (mm).

```python
.. include:: ../examples/geometry/08_projection_matrix_example3.py
```

## Simulating a complete CT scan

A single projection matrix is not enough to describe a full CT scan. We need one projection matrix for each frame (i.e., for each projection image).

We can use a loop to set up each frame and collect the projection matrices in a list. Afterwards, we can pass this list of matrices to the function `writeOpenCTFile()` or `writeCERAconfig()` to create specific reconstruction configuration files for each reconstruction software.

In the loop, it is advisable not to rotate the stage incrementally for each frame by a certain angular increment. This could lead to the accumulation of small floating-point rounding inaccuracies. Instead, we create a backup of the initial setup (at frame zero) using the `Geometry.store()` function. In each step of the loop, we restore this initial configuration by calling `Geometry.restore()` and then rotate the stage to its current absolute angle. This approach of parameterizing the whole CT trajectory as a deterministic function that only depends on the initial configuration and the current frame number is preferred over incremental changes in a loop, but might not always be feasible.

The following example shows how to simulate a simple CT scan (one full stage rotation with 3000 equidistant projection images) and how to create configuration files for the reconstruction.

```python
.. include:: ../examples/geometry/09_projection_matrix_full_CT.py
```