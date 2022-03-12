# -*- coding: UTF-8 -*-
"""
Coordinate systems, transformations and projection matrix functionality.

Coordinate Systems
==================

The geometry subpackage provides a `CoordinateSystem` class that lets you create and manipulate objects in a virtual CT scene. Such an oject has a position (`center`) and three basis vectors (`u`, `v`, `w`). These basis vectors are assumed to be orthogonal, but they do not have to be unit vectors. The `center` acts as the pivot point for rotations.

Example for creating and manipulating an object:

```python
from ctsimu.geometry import *
from ctsimu.helpers  import *  # provides deg2rad()

mySpecimen = CoordinateSystem()

# Set position and orientation:
mySpecimen.center = Vector(250, 0, 0)
mySpecimen.u = Vector(0, -1,  0)
mySpecimen.v = Vector(0,  0, -1)
mySpecimen.w = Vector(1,  0,  0)

# Manipulate:
mySpecimen.translate(translationVector = Vector(5.2, 0, 4.3))
mySpecimen.rotateAroundU(angle = deg2rad(2.0))
mySpecimen.rotate(axis = Vector(1, 1, 1), angle = deg2rad(5.0))

print("My specimen's new location and orientation:")
print(mySpecimen)
```

Full CT Geometry
================

For a full CT, we need an X-ray source, a stage for the specimens, and a detector. The `Geometry` class bundles three coordinate systems (one for each of those components), and additional information about the detector (using the `Detector` class, an extension of a regular `CoordinateSystem`). The following figure shows their standard orientations when a CT geometry is initialized.

![Standard coordinate system](pictures/geometry.png "Standard coordinate system")

The orientation of the coordinate system and all components can be changed by rotations or by manually setting the object basis vectors. However, it is important to keep the following conventions.

Detector convention
-------------------

* The detector's `u` vector is its row vector.
* The detector's `v` vector is its column vector.
* The detector's `w` vector has no special meaning. It is a planar normal that must be chosen such that the detector's coordinate system remains right-handed.

Stage convention
----------------

* The stage's `w` vector is its axis of CT rotation.

Source convention
-----------------

There is currently no restriction on the source coordinate system. We usually assume its `w` axis to be the direction of the principal ray, but this is not a necessity.

Example Setup
-------------

In the following example, we set up a standard CT geometry.

```python
from ctsimu.geometry import *

# General CT parameters:
SOD = 250.0 # mm
SDD = 800.0 # mm

pixelSize = 0.2 # mm
pixelColumns = 2000
pixelRows = 1000

# Create a CT geometry object:
myCT = Geometry()

# Stage:
myCT.stage.center.x = SOD

# Detector:
myCT.detector.center.x = SDD
myCT.detector.setSize(
    pixelsU = pixelColumns,
    pixelsV = pixelRows,
    pitchU = pixelSize,
    pitchV = pixelSize
    )

myCT.update() # calculates derived geometry parameters

print(myCT.info())
```

Reference Frames
================

Implicitly, each coordinate system has a reference coordinate system (its *reference frame*) in which its `center` and `u`, `v`, `w` basis vectors are located and described. Typically, we assume that this is a standard coordinate system called the **world coordinate system**. Any new `CoordinateSystem` object is initialized to be aligned with the world coordinate system:

```python
from ctsimu.geometry import *

world = CoordinateSystem()

print("The World:")
print(world)

\"\"\"
The World:
Center: ( 0.0000000,  0.0000000,  0.0000000)
u:      ( 1.0000000,  0.0000000,  0.0000000)
v:      ( 0.0000000,  1.0000000,  0.0000000)
w:      ( 0.0000000,  0.0000000,  1.0000000)
\"\"\"
```

You can change the reference frame of a `CoordinateSystem`. In the following example, we set up a CT geometry with a stage that is tilted by 2°. We place a specimen object in the stage coordinate system and move it "upwards" by 5 mm along the (now tilted) axis of rotation. Afterwards, we change the specimen's reference frame to see where it is actually located in the world coordinate system.

```python
from ctsimu.geometry import *
from ctsimu.helpers  import *  # provides deg2rad()

world = CoordinateSystem()

# Set up a quick CT geometry with a tilted stage axis:
myCT = Geometry()
myCT.stage.center.x = 250  # SOD
myCT.stage.rotateAroundU(angle = deg2rad(2.0))
myCT.detector.center.x = 800  # SDD

# Assume a specimen in the (tilted) stage
# coordinate system, shifted 5 mm "upwards"
# along the axis of rotation:
mySpecimen = CoordinateSystem()
mySpecimen.translateZ(5.0)

# Change the specimen's reference frame to
# the world coordinate system:
mySpecimen.changeReferenceFrame(
    fromCS = myCT.stage,
    toCS = world)

print("The specimen's world coordinates:")
print(mySpecimen)

\"\"\"
The specimen's world coordinates:
Center: ( 250.0000000, -0.1744975,  4.9969541)
u:      ( 1.0000000,  0.0000000,  0.0000000)
v:      ( 0.0000000,  0.9993908,  0.0348995)
w:      ( 0.0000000, -0.0348995,  0.9993908)
\"\"\"
```

**Note:** when changing reference frames, the original and the target reference frame must both have a common reference frame for themselves. In the example above, we change the reference frame from the stage coordinate system to the world coordinate system. Both of them have the same reference frame: the world coordinate system (which is special, because it is also a reference for itself).


Projection Matrices
===================

A projection matrix maps a 3D point coordinate `(x, y, z)` from the stage coordinate system to a 2D point coordinate `(u, v)` in the detector coordinate system. They are used by some reconstruction softwares to describe arbitrary scan trajectories. For such a reconstruction, we need one projection matrix for each projection image.

![Euclidean Mapping](pictures/pmatrix_mapping_euclidean.png "Mapping a coordinate from the stage coordinate system to a coordinate in the detector coordinate system")

Mathematical Background
-----------------------

We operate in [homogeneous coordinates](https://en.wikipedia.org/wiki/Homogeneous_coordinates), because we are in a projective geometry and this concept allows us to describe translations in space by a matrix. Homogeneous coordinates describe rays in a space and are therefore only defined up to a scale factor, which is carried as an additional coordinate. This way, a Euclidean 3-vector turns into a homogeneous 4-vector. Our mapping becomes:

![Homogeneous Mapping](pictures/pmatrix_mapping_homogeneous.png "Mapping in homogeneous coordinates")

The following picture illustrates a 1D projective geometry. The *h* axis is our scale factor for the homogeneous coordinates. In this geometry, all points on a ray are equivalent. The points at *h*=1 are the normalized homogeneous coordinates.

![Projective geometry for a 1D Euclidean space](pictures/hom_coords.png "Projective geometry for a 1D Euclidean space")

We can use this concept to describe translations in space with a matrix multiplication:

![Translation Matrix](pictures/pmatrix_translation.png "Translation Matrix")

We now consider two coordinate systems: the stage coordinate system (our origin) and the detector coordinate system (our target). The world coordinate system is not important in this context. In the following picture, 3D coordinates are expressed in terms of the stage coordinate system, and 2D coordinates are expressed in terms of the detector coordinate system (its *w* axis is irrelevant here).

![Projective system](pictures/geometry_stage_projection.png "Projective system")

To calculate a projection matrix for the current geometry, we have to consider five subsequent transformations. Each transformation is expressed by a matrix. The final projection matrix is then the product of these five transformation matrices.

1. We shift the origin of the coordinate system from the stage to the source, which is our center of projection.

    ![Shifting the origin to the source](pictures/pmatrix_F.png "Shifting the origin to the source")

2. We perform a basis transformation to express the 3D coordinates in terms of the axes of the detector coordinate system. The origin remains at the source. You can also think of this transformation as a rotation into the detector coordinate system.

    ![Rotation into detector coordinate system](pictures/pmatrix_R.png "Rotation into detector coordinate system")

    All basis vectors in this matrix are assumed to be unit vectors.

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

"""

import numpy
import os    # File and path handling
import json
import math
import copy
import pkgutil
from datetime import datetime

from .primitives import *
from .image import Image  # To create detector flat field

def basisTransformMatrix(fromCS:'CoordinateSystem', toCS:'CoordinateSystem') -> Matrix:
    """Calculate a matrix that transforms coordinates from `fromCS` to `toCS`.

    `fromCS` and `toCS` must have the same common reference frame
    (e.g. the world coordinate system).

    Parameters
    ----------
    fromCS : CoordinateSystem
        The origin coordinate system.

    toCS : CoordinateSystem
        The target coordinate system.

    Returns
    -------
    T : Matrix
        The 3x3 basis transformation matrix.

    References
    ----------
    * S. Widnall: [Lecture L3 - Vectors, Matrices and Coordinate Transformations]
    [Lecture L3 - Vectors, Matrices and Coordinate Transformations]: https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-07-dynamics-fall-2009/lecture-notes/MIT16_07F09_Lec03.pdf
    """

    T = Matrix(3, 3)

    # Row 1:
    T.value[0][0] = toCS.u.unitVector().dot(fromCS.u.unitVector())
    T.value[0][1] = toCS.u.unitVector().dot(fromCS.v.unitVector())
    T.value[0][2] = toCS.u.unitVector().dot(fromCS.w.unitVector())

    # Row 2:
    T.value[1][0] = toCS.v.unitVector().dot(fromCS.u.unitVector())
    T.value[1][1] = toCS.v.unitVector().dot(fromCS.v.unitVector())
    T.value[1][2] = toCS.v.unitVector().dot(fromCS.w.unitVector())

    # Row 3:
    T.value[2][0] = toCS.w.unitVector().dot(fromCS.u.unitVector())
    T.value[2][1] = toCS.w.unitVector().dot(fromCS.v.unitVector())
    T.value[2][2] = toCS.w.unitVector().dot(fromCS.w.unitVector())

    return T

class CoordinateSystem:
    """Coordinate system: center point and axis vectors.

    An object according to the CTSimU scenario specification,
    containing a center coordinate and an orientation in 3D space.
    
    The center and axis vectors are expressed in terms of the
    object's reference coordinate system, which must be known implicitly
    when objects of this class are used.

    Geometrical objects could be source, stage or detector.
    Samples would need additional attention due to possible attachment
    to stage coordinate system (instead of world).

    Attributes
    ----------
    center : Vector
        The location of the center point in a reference
        coordinate system (usually world or stage).

    u : Vector
        Basis vector for the u axis.

    v : Vector
        Basis vector for the v axis.

    w : Vector
        Basis vector for the w axis.
    """

    def __init__(self):
        """Initialize as a standard world coordinate system."""
        self.center  = Vector(0, 0, 0)
        self.u = Vector(1, 0, 0)
        self.v = Vector(0, 1, 0)
        self.w = Vector(0, 0, 1)

    def __str__(self):
        """Information string for easy printing."""

        txt  = "Center: {}\n".format(self.center)
        txt += "u:      {}\n".format(self.u)
        txt += "v:      {}\n".format(self.v)
        txt += "w:      {}\n".format(self.w)
        return txt

    def json_import(self, geometry: dict):
        """Set up geometry from a JSON dictionary.

        Parameters
        ----------
        geometry : dict
            A parsed JSON dictionary from a [CTSimU scenario description] file.

        [CTSimU scenario description]: https://bamresearch.github.io/ctsimu-scenarios/

        Raises
        ------
        KeyError
            When expected JSON keys for center and vector x, y, z
            components are not found in the dictionary.
        """

        # Get center position from JSON dict:
        if "center" in geometry:
            if "x" in geometry["center"]:
                cx = in_mm(geometry["center"]["x"])
            else:
                raise KeyError("No \"x\" coordinate found for the center.")

            if "y" in geometry["center"]:
                cy = in_mm(geometry["center"]["y"])
            else:
                raise KeyError("No \"y\" coordinate found for the center.")
            
            if "z" in geometry["center"]:
                cz = in_mm(geometry["center"]["z"])
            else:
                raise KeyError("No \"z\" coordinate found for the center.")

            self.center.set(cx, cy, cz)
        # Try old British spelling (up to file format v0.9)
        elif "centre" in geometry:
            if "x" in geometry["centre"]:
                cx = in_mm(geometry["centre"]["x"])
            else:
                raise KeyError("No \"x\" coordinate found for the center.")

            if "y" in geometry["centre"]:
                cy = in_mm(geometry["centre"]["y"])
            else:
                raise KeyError("No \"y\" coordinate found for the center.")
            
            if "z" in geometry["centre"]:
                cz = in_mm(geometry["centre"]["z"])
            else:
                raise KeyError("No \"z\" coordinate found for the center.")

            self.center.set(cx, cy, cz)
        else:
            raise KeyError("JSON file is missing a geometry \"center\" section.")

        # Get vector u from JSON dict:
        if "vector_u" in geometry:
            if "x" in geometry["vector_u"]:
                ux = geometry["vector_u"]["x"]
            else:
                raise KeyError("No \"x\" component found for vector u.")

            if "y" in geometry["vector_u"]:
                uy = geometry["vector_u"]["y"]
            else:
                raise KeyError("No \"y\" component found for vector u.")

            if "z" in geometry["vector_u"]:
                uz = geometry["vector_u"]["z"]
            else:
                raise KeyError("No \"z\" component found for vector u.")
        else:
            raise KeyError("JSON file is missing a geometry \"vector_u\" section.")

        # Get vector w from JSON dict:
        if "vector_w" in geometry:
            if "x" in geometry["vector_w"]:
                wx = geometry["vector_w"]["x"]
            else:
                raise KeyError("No \"x\" component found for vector w.")

            if "y" in geometry["vector_w"]:
                wy = geometry["vector_w"]["y"]
            else:
                raise KeyError("No \"y\" component found for vector w.")

            if "z" in geometry["vector_w"]:
                wz = geometry["vector_w"]["z"]
            else:
                raise KeyError("No \"z\" component found for vector w.")
        else:
            raise KeyError("JSON file is missing a geometry \"vector_w\" section.")

        # Set up the geometry from the information given in the JSON file:
        c = Vector(cx, cy, cz)  # center
        u = Vector(ux, uy, uz)  # u basis vector
        w = Vector(wx, wy, wz)  # w basis vector
        v = w.cross(u)
        self.setup(c, u, v, w)
        self.makeUnitCS()

        # Apply deviations from the now-ideal geometry:
        if "deviation" in geometry:
            # Positional deviations:
            if "position" in geometry["deviation"]:
                if "x" in geometry["deviation"]["position"]:
                    if geometry["deviation"]["position"]["x"] != None:
                        translationX = in_mm(geometry["deviation"]["position"]["x"])
                        self.translateX(translationX)
    
                if "y" in geometry["deviation"]["position"]:
                    if geometry["deviation"]["position"]["y"] != None:
                        translationY = in_mm(geometry["deviation"]["position"]["y"])
                        self.translateY(translationY)

                if "z" in geometry["deviation"]["position"]:
                    if geometry["deviation"]["position"]["z"] != None:
                        translationZ = in_mm(geometry["deviation"]["position"]["z"])
                        self.translateZ(translationZ)

            # Rotations according to w''v'u convention:
            if "rotation" in geometry["deviation"]:
                if "w" in geometry["deviation"]["rotation"]:
                    if geometry["deviation"]["rotation"]["w"] != None:
                        angleAroundW = in_rad(geometry["deviation"]["rotation"]["w"])
                        self.rotateAroundW(angleAroundW)

                if "v" in geometry["deviation"]["rotation"]:
                    if geometry["deviation"]["rotation"]["v"] != None:
                        angleAroundV = in_rad(geometry["deviation"]["rotation"]["v"])
                        self.rotateAroundV(angleAroundV)

                if "u" in geometry["deviation"]["rotation"]:
                    if geometry["deviation"]["rotation"]["u"] != None:
                        angleAroundU = in_rad(geometry["deviation"]["rotation"]["u"])
                        self.rotateAroundU(angleAroundU)

    def setup(self, center:Vector, u:Vector, v:Vector, w:Vector):
        """Set up center and orientation manually.

        Parameters
        ----------
        center : Vector
            Object's center point in reference coordinate system,
            origin of local {u,v,w} coordinate system.

        u : Vector
            Basis vector u in terms of reference coordinate system.

        v : Vector
            Basis vector v in terms of reference coordinate system.

        w : Vector
            Basis vector w in terms of reference coordinate system.

        Notes
        -----
        All basis vectors must be orthogonal.
        """

        # Create new vectors from given components:
        self.center  = center
        self.u = u
        self.v = v
        self.w = w

    def update(self):
        """Signal a manual update to the center position or orientation vectors."""
        self.center.update()
        self.u.update()
        self.v.update()
        self.w.update()

    def makeUnitCS(self):
        """Convert all basis vectors to unit vectors."""
        self.u.makeUnitVector()
        self.v.makeUnitVector()
        self.w.makeUnitVector()

    def translate(self, translationVector: Vector):
        """Move object in space.

        Parameters
        ----------
        translationVector : Vector
            Vector by which the object's center point should be shifted.
            Its components are added to the center's components.
        """
        self.center.add(translationVector)

    def translateX(self, dx: float):
        """Move object in x direction.

        Parameters
        ----------
        dx : float
            Shift amount in x direction.
        """
        self.center.setx(self.center.x + float(dx))

    def translateY(self, dy: float):
        """Move object in y direction.

        Parameters
        ----------
        dy : float
            Shift amount in y direction.
        """
        self.center.sety(self.center.y + float(dy))

    def translateZ(self, dz: float):
        """Move object in z direction.

        Parameters
        ----------
        dz : float
            Shift amount in z direction.
        """
        self.center.setz(self.center.z + float(dz))

    def rotateAroundU(self, angle: float):
        """Rotate object around its u axis by given angle [rad].
        
        Parameters
        ----------
        angle : float
            Rotation angle in rad, mathematically positive direction (right-hand rule).
        """
        self.v.rotate(self.u, angle)
        self.w.rotate(self.u, angle)

    def rotateAroundV(self, angle: float):
        """Rotate object around its v axis by given angle [rad].
        
        Parameters
        ----------
        angle : float
            Rotation angle in rad, mathematically positive direction (right-hand rule).
        """
        self.u.rotate(self.v, angle)
        self.w.rotate(self.v, angle)

    def rotateAroundW(self, angle: float):
        """Rotate object around its w axis by given angle [rad].
        
        Parameters
        ----------
        angle : float
            Rotation angle in rad, mathematically positive direction (right-hand rule).
        """
        self.u.rotate(self.w, angle)
        self.v.rotate(self.w, angle)

    def rotate(self, axis: Vector, angle: float):
        """Rotate object around a given axis by the given angle [rad].
        
        Parameters
        ----------
        axis : Vector
            The axis of rotation, in terms of the object's
            reference coordinate system (e.g. world).
        
        angle : float
            Rotation angle in rad, mathematically positive direction (right-hand rule).
        """
        self.u.rotate(axis, angle)
        self.v.rotate(axis, angle)
        self.w.rotate(axis, angle)

    def changeReferenceFrame(self, fromCS:'CoordinateSystem', toCS:'CoordinateSystem'):
        """Change the object's reference coordinate system.
        
        Parameters
        ----------
        fromCS : CoordinateSystem
            Current reference coordinate system.
        
        toCS : CoordinateSystem
            New reference coordinate system.

        Notes
        -----
        Both fromCS and toCS must be in the same reference coordinate system
        (e.g., the world coordinate system).
        """

        # Rotate basis vectors into toCS:
        T = basisTransformMatrix(fromCS, toCS)
        self.u = T * self.u
        self.v = T * self.v
        self.w = T * self.w

        world = CoordinateSystem()

        # Move center to toCS:
        # 1. Translate self.center by difference of toCS and fromCS
        #    -> Origins are "superimposed".
        # 2. Rotate self.center from fromCS to toCS.

        # Translation vector in world coordinates:
        translator = fromCS.center - toCS.center  # in world coordinates
        # Translation vector in fromCS coordinates:
        M = basisTransformMatrix(world, fromCS)
        translator = M*translator
        relCenter = self.center + translator
        self.center = T*relCenter

class Detector(CoordinateSystem):
    """Detector as geometrical object.

    With additional attributes for the spatial extension and
    the pixel coordinate system.

    Attributes
    ----------
    pixelsU : int
        Number of pixels in u direction.
    
    pixelsV : int
        Number of pixels in v direction.
    
    pitchU : float
        Size of a pixel in u direction.
        In units of the reference coordinate system.
    
    pitchV : float
        Size of a pixel in v direction.
        In units of the reference coordinate system.
    
    physWidth : float
        Physical size in u direction.
        In units of the reference coordinate system.
        Computed automatically after calling `setSize()`.
    
    physHeight : float
        Physical size in v direction.
        In units of the reference coordinate system.
        Computed automatically after calling `setSize()`.
 
    pixelOrigin : Vector
        Origin of the pixel coordinate system in terms of the reference
        coordinate system. This is the outermost corner of the
        (0,0) pixel of the detector (often the "upper left" corner).
        Computed automatically after calling `setSize()`.

    Notes
    -----
    Use `setSize()` to set the size of the detector, given its number of pixels
    and the pitch. This function automatically computes the physical dimensions
    `physWidth` and `physHeight` and the origin of the pixel coordinate system.
    """

    def __init__(self):
        """Initialize as a standard CoordinateSystem.

        Orientation, position and size must be set up manually afterwards.
        """

        # Call init from parent class:
        CoordinateSystem.__init__(self)

        self.pixelsU     = None  # Detector pixels in u direction
        self.pixelsV     = None  # Detector pixels in v direction
        self.pitchU = None  # Size of a pixel in u direction in units of reference coordinate system
        self.pitchV = None  # Size of a pixel in v direction in units of reference coordinate system
        self.physWidth   = 0    # Physical width in units of reference coordinate system
        self.physHeight  = 0    # Physical height in units of reference coordinate system

        self.pixelOrigin = Vector()  # origin of pixel coordinate system in terms of reference coordinate system

    def setSize(self, pixelsU:int = None, pixelsV:int = None, pitchU:float = None, pitchV:float = None):
        """Set the physical size of the detector.

        From the given parameters (number of pixels and pitch), the physical
        size of the detector and the position of the origin of the pixel
        coordinate system will be calculated. Make sure that the orientation
        vectors and position of the detector are correct before calling
        `setSize()`, or call `computeGeometryParameters()` if you update
        the detector orientation or position later on.

        Parameters
        ----------
        pixelsU : int
            Number of pixels in u direction.

        pixelsV : int
            Number of pixels in v direction.

        pitchU : float
            Pixel pitch in u direction.

        pitchV : float
            Pixel pitch in v direction.
        """

        self.pixelsU = int(pixelsU)
        self.pixelsV = int(pixelsV)
        self.pitchU = float(pitchU)
        self.pitchV = float(pitchV)

        self.computeGeometryParameters()

    def computeGeometryParameters(self):
        """Calculate the physical width and height, and the position of the
        pixel coordinate system origin.

        These calculations assume that the size, position and
        orientation of the detector are correctly set up.

        Results are assigned to their member variables (attributes).
        """

        if (self.pixelsU is not None) and (self.pixelsV is not None) and (self.pitchU is not None) and (self.pitchV is not None):
            # Physical width and height:
            self.physWidth  = self.pixelsU * self.pitchU
            self.physHeight = self.pixelsV * self.pitchV

            # Vectors of the detector coordinate system:
            ux = self.u.unitVector().x
            uy = self.u.unitVector().y
            uz = self.u.unitVector().z
            vx = self.v.unitVector().x
            vy = self.v.unitVector().y
            vz = self.v.unitVector().z

            # World coordinates of origin (0,0) of detector's pixel coordinate system:
            self.pixelOrigin.x = self.center.x - 0.5*(ux*self.physWidth + vx*self.physHeight)
            self.pixelOrigin.y = self.center.y - 0.5*(uy*self.physWidth + vy*self.physHeight)
            self.pixelOrigin.z = self.center.z - 0.5*(uz*self.physWidth + vz*self.physHeight)

    def cols(self) -> int:
        """Returns the number of detector columns (i.e., pixels in u direction).

        Returns
        -------
        pixelsU : int
            Number of detector columns (i.e., pixels in u direction).
        """
        return self.pixelsU

    def rows(self) -> int:
        """Returns the number of detector rows (i.e., pixels in v direction).

        Returns
        -------
        pixelsV : int
            Number of detector rows (i.e., pixels in v direction).
        """
        return self.pixelsV

    def pixelVector(self, x: float, y: float) -> Vector:
        """World position vector for given pixel coordinate.

        The pixel coordinate system has its origin at the detector corner with
        the lowest coordinate in terms of its u and v basis vectors. Typically,
        this is the upper left corner, but your arrangement may differ.

        Integer coordinates always refer to the pixel corner that is closest to
        the origin of the pixel coordinate system, whereas the center of a pixel
        therefore has a ".5" coordinate in the pixel coordinate system.
        For example, the first pixel (0, 0) would have center coordinates
        (0.5, 0.5).

        To get the center coordinates for a given integer pixel location,
        `pixelVectorCenter()` may be used.

        Parameters
        ----------
        x : float
            x position in pixel coordinate system.

        y : float
            y position in pixel coordinate system.

        Returns
        -------
        pixelVector : Vector
            Pixel position in reference coordinate system (usually world)
            as a 3D vector.
        """

        # x, y are coordinates in pixel coordinates system
        px = self.pixelOrigin.x + self.u.x*x*self.pitchU + self.v.x*y*self.pitchV
        py = self.pixelOrigin.y + self.u.y*x*self.pitchU + self.v.y*y*self.pitchV
        pz = self.pixelOrigin.z + self.u.z*x*self.pitchU + self.v.z*y*self.pitchV
        pixelVector = Vector(px, py, pz)
        return pixelVector

    def pixelVectorCenter(self, x: float, y: float) -> Vector:
        """World position vector of pixel center, for a pixel given in integer coordinates.

        Parameters
        ----------
        x : float
            Integer x coordinate, specifies a pixel in the pixel coordinate system.

        y : float
            Integer y coordinate, specifies a pixel in the pixel coordinate system.

        Returns
        -------
        pixelVector : Vector
            Position of the pixel center in the reference coordinate system
            (usually world) as a 3D vector.

        Notes
        -----
        If `float` coordinates are passed (non-integer),
        they are converted to integers using `math.floor`.
        """
        return self.pixelVector(float(math.floor(x))+0.5, float(math.floor(y))+0.5)


class Geometry:
    """Bundles geometry information about the complete CT setup.

    Keeps the source, stage and detector as a set and provides methods
    to calculate geometry parameters and projection matrices.
    
    Attributes
    ----------
    detector : Detector
        The detector geometry.

    source : CoordinateSystem
        The source geometry.

    stage : CoordinateSystem
        The stage geometry.

    SDD : float
        Shortest distance between source center and detector plane.
        Calculated automatically by `update()`.

    SOD : float
        Distance between source center and stage center.
        Calculated automatically by `update()`.

    ODD : float
        Shortest distance between stage center and detector plane.
        Calculated automatically by `update()`.

    brightestSpotWorld : Vector
        Location of the intensity maximum on the detector, in world coordinates.
        Assuming an isotropically radiating source.
        Calculated automatically by `update()`.

    brightestSpotDetector : Vector
        Location of the intensity maximum on the detector, in terms of
        detector coordinate system. Assuming an isotropically radiating source.
        Calculated automatically by `update()`.
    """

    def __init__(self, jsonFile:str = None, jsonFileFromPkg:str = None):
        """Initialize using the provided JSON geometry specification.

        Parameters
        ----------
        jsonFile : str, optional
            Location of a CTSimU scenario description file to
            import the geometry.

        jsonFileFromPkg : str, optional
            Load the geometry from a JSON file included in the package,
            usually for internal purposes only.

        """
        self.detector    = Detector()
        self.source      = CoordinateSystem()
        self.stage       = CoordinateSystem()

        # Initialize source and detector to standard CTSimU orientation:
        self.detector.u = Vector(0, -1,  0)
        self.detector.v = Vector(0,  0, -1)
        self.detector.w = Vector(1,  0,  0)
        self.source.u   = Vector(0, -1,  0)
        self.source.v   = Vector(0,  0, -1)
        self.source.w   = Vector(1,  0,  0)

        self.SDD = None
        self.SOD = None
        self.ODD = None
        self.brightestSpotWorld = None
        self.brightestSpotDetector = None

        jsonText = None
        if jsonFileFromPkg is not None:  # from package
            jsonFile = jsonFileFromPkg
            jsonText = pkgutil.get_data(__name__, jsonFileFromPkg).decode()
        elif jsonFile is not None:  # from file
            if os.path.isfile(jsonFile):
                log("JSON File: {}".format(jsonFile))
                jsonFilePtr = open(jsonFile, "r")
                jsonText = jsonFilePtr.read()
                jsonFilePtr.close()
            else:
                raise Exception("Can't find " + jsonFile)
        else:
            return

        if(jsonText is not None):
            try:
                jsonDict = json.loads(jsonText)
            except:
                raise Exception("Error parsing JSON file: {}".format(jsonFile))

            # Detector size and pixel pitch:
            pixelsU = getFieldOrNone(jsonDict, "detector", "columns", "value")
            pixelsV = getFieldOrNone(jsonDict, "detector", "rows", "value")
            pitchU = getFieldOrNone(jsonDict, "detector", "pixel_pitch", "u", "value")
            pitchV = getFieldOrNone(jsonDict, "detector", "pixel_pitch", "v", "value")

            try:
                detectorGeometry = getFieldOrNone(jsonDict, "geometry", "detector")
                if detectorGeometry != None:
                    self.detector.json_import(detectorGeometry)
                    self.detector.setSize(pixelsU, pixelsV, pitchU, pitchV)
                else:
                    raise Exception("JSON file does not contain a valid \"detector\" section in \"geometry\".")
            except Exception as e:
                log("Something went wrong when setting up the detector geometry using the JSON file description.")
                raise Exception(e)

            try:
                sourceGeometry = getFieldOrNone(jsonDict, "geometry", "source")
                if sourceGeometry != None:
                    self.source.json_import(sourceGeometry)
                else:
                    raise Exception("JSON file does not contain a valid \"source\" section in \"geometry\".")
            except Exception as e:
                log("Something went wrong when setting up the source geometry using the JSON file description.")
                raise Exception(e)

            try:
                stageGeometry = getFieldOrNone(jsonDict, "geometry", "stage")
                if stageGeometry != None:
                    self.stage.json_import(stageGeometry)
                else:
                    raise Exception("JSON file does not contain a valid \"stage\" section in \"geometry\".")
            except Exception as e:
                log("Something went wrong when setting up the stage geometry using the JSON file description.")
                raise Exception(e)

            self.update()            
        else:
            raise Exception("JSON scenario file not available.")

    def __str__(self):
        return self.info()
        
    def update(self):
        """Calculate derived geometry parameters.

        Calculates the SOD, SDD, ODD, and location of the intensity maximum
        on the detector (in world and detector coordinates) for the
        curent geometry. Results are stored in the following member variables
        (attributes).

        SDD: Shortest distance between source center and detector plane.

        SOD: Distance between source center and stage center.

        ODD: Shortest distance between stage center and detector plane.

        brightestSpotWorld: Location of the intensity maximum on the detector,
            in world coordinates.  Assuming an isotropically radiating source.

        brightestSpotDetector: Location of the intensity maximum on the
            detector, in terms of detector coordinate system.
            Assuming an isotropically radiating source.
        """

        self.source.update()
        self.stage.update()
        self.detector.update()

        # SOD, SDD, ODD
        world = CoordinateSystem()
        source_from_image = copy.deepcopy(self.source)
        stage_from_detector  = copy.deepcopy(self.stage)

        source_from_image.changeReferenceFrame(world, self.detector)
        stage_from_detector.changeReferenceFrame(world, self.detector)

        self.SDD = abs(source_from_image.center.z)
        self.ODD = abs(stage_from_detector.center.z)
        self.SOD = self.source.center.distance(self.stage.center)

        ## Brightest Spot in World Coordinate System:
        self.brightestSpotWorld = copy.deepcopy(self.detector.w)
        self.brightestSpotWorld.scale(self.SDD)
        self.brightestSpotWorld.add(self.source.center)

        ## Brightest Spot in Detector Coordinate System:
        self.brightestSpotDetector = copy.deepcopy(self.brightestSpotWorld)
        self.brightestSpotDetector.subtract(self.detector.center)
        
        pxU = self.brightestSpotDetector.dot(self.detector.u) / self.detector.pitchU + self.detector.cols()/2.0
        pxV = self.brightestSpotDetector.dot(self.detector.v) / self.detector.pitchV + self.detector.rows()/2.0

        self.brightestSpotDetector = Vector(pxU, pxV, 0)


    def info(self) -> str:
        """Generate an information string about the current geometry.
    
        Returns
        -------
        txt : string
            Information string for humans.
        """

        self.update()

        txt  = "Detector\n"
        txt += "===========================================================\n"
        txt += "Center:          {}\n".format(self.detector.center)
        txt += "u:               {}\n".format(self.detector.u)
        txt += "v:               {}\n".format(self.detector.v)
        txt += "w:               {}\n".format(self.detector.w)
        txt += "Pixels:          {cols} x {rows}\n".format(cols=self.detector.cols(), rows=self.detector.rows())
        txt += "Pitch:           {pitchU} x {pitchV}\n".format(pitchU=self.detector.pitchU, pitchV=self.detector.pitchV)
        txt += "Physical Size:   {width} x {height}\n".format(width=self.detector.physWidth, height=self.detector.physHeight)

        txt += "Brightest Spot:\n"
        txt += "  World:         {}\n".format(self.brightestSpotWorld)
        txt += "  Pixels:        {}\n".format(self.brightestSpotDetector)

        txt += "\n"
        txt += "Source:\n"
        txt += "===========================================================\n"
        txt += "Center:          {}\n".format(self.source.center)
        txt += "u:               {}\n".format(self.source.u)
        txt += "v:               {}\n".format(self.source.v)
        txt += "w:               {}\n".format(self.source.w)

        txt += "\n"
        txt += "Stage:\n"
        txt += "===========================================================\n"
        txt += "Center:          {}\n".format(self.stage.center)
        txt += "u:               {}\n".format(self.stage.u)
        txt += "v:               {}\n".format(self.stage.v)
        txt += "w:               {}\n".format(self.stage.w)

        txt += "\n"
        txt += "Geometry Parameters:\n"
        txt += "===========================================================\n"
        # Source - Detector distance (SDD) defined by shortest distance between source and detector:
        txt += "SDD:             {}\n".format(self.SDD)
        txt += "ODD:             {}\n".format(self.ODD)
        txt += "SOD:             {}\n".format(self.SOD)

        return txt

    def projectionMatrix(self,
                         volumeCS:CoordinateSystem=None,
                         imageCS:CoordinateSystem=None,
                         mode:str=None,
                         mirror:bool=True):
        """Calculate a projection matrix for the current geometry.

        Parameters
        ----------
        volumeCS : CoordinateSystem, optional
            Position of the volume coordinate system in terms of the
            stage coordinate system. If `None` is given, the volume
            coordinate system is assumed to be the stage coordinate system.
            See notes for details.

        imageCS : CoordinateSystem, optional
            Position of the image coordinate system in terms of the
            detector coordinate system. If `None` is given, the image
            coordinate system is assumed to be the detector coordinate system.
            See notes for details.

        mode : str, optional
            Pre-defined modes. Either "openCT" or "CERA" are supported.
            They override the `volumeCS` and `imageCS`, which can be set
            to `None` when using one of the pre-defined modes.

        mirror : bool, optional
            Whether or not a mirror operation should be applied to the
            reconstruction volume. For many 3D processing softwares,
            this parameter should be set to `True` to avoid having to
            mirror the volume after loading it into the software.

        Returns
        -------
        P : Matrix
            Projection matrix.

        Notes
        -----
        The image coordinate system (`imageCS`) should match the notation
        used by the reconstruction software, and is expressed in terms of
        the detector coordinate system.

        The detector coordinate system has its origin at the detector center,
        the u unit vector points in the row vector direction, and the
        v unit vector points in column vector direction (they are always assumed
        to be unit vectors).

        The center (origin) of the `imageCS` should be where the reconstruction
        software places the origin of its own projection image coordinate
        system. For example, CERA places it at the center of the lower-left pixel
        of the projection image.
        """

        validModes = ["openCT", "CERA"]

        if mode is not None:
            if mode in validModes:  # Override imageCS
                image = CoordinateSystem()

                if mode == "openCT":
                    """openCT places the origin of the image CS at the detector 
                    center. The constructor places it at (0,0,0) automatically,
                    so there is nothing to do. Comments for illustration."""
                    # image.center.x = 0
                    # image.center.y = 0
                    # image.center.z = 0

                    """openCT's image CS is in mm units. We assume that all
                    other coordinate systems are in mm as well here (at least
                    when imported from JSON file). No scaling of the basis vectors is necessary."""
                    # image.u.scale(1.0)
                    # image.v.scale(1.0)
                    # image.w.scale(1.0)

                elif mode == "CERA":
                    """CERA places the origin of the image CS in the center
                    of the lower left pixel of the projection image."""
                    image.center.x = -self.detector.physWidth  / 2.0 + 0.5*self.detector.pitchU
                    image.center.y =  self.detector.physHeight / 2.0 - 0.5*self.detector.pitchV
                    # image.center.z = 0

                    """CERA's unit of the image CS is in px, so we need to
                    scale the image CS basis vectors by the pixel size.
                    Also, v points up instead of down."""
                    image.u.scale( self.detector.pitchU)
                    image.v.scale(-self.detector.pitchV)
                    image.w.scale(-1.0)
            else:
                raise RuntimeError("Unsupported mode for projection matrix: \"{}\"".format(mode))
        elif imageCS is not None:
            image = copy.deepcopy(imageCS)
        else:
             # Set a standard coordinate system. Results in pure
             # detector coordinate system after transformation.
            image = CoordinateSystem()

        world    = CoordinateSystem()
        source   = copy.deepcopy(self.source)

        # The 3D volume (reconstruction space).
        volume = None
        if volumeCS is not None:
            volume = copy.deepcopy(volumeCS)
            volume = volume.changeReferenceFrame(self.stage, world)
        else:
            volume = copy.deepcopy(self.stage)

        """The scale factors are derived from the lengths of the basis
        vectors of the image CS compared to the detector CS unit vectors."""
        scale_volume_u = self.stage.u.unitVector().dot(volume.u)
        scale_volume_v = self.stage.v.unitVector().dot(volume.v)
        scale_volume_w = self.stage.w.unitVector().dot(volume.w)

        # Detach the image CS from the detector CS and
        # express it in terms of the world CS:
        image.changeReferenceFrame(self.detector, world)

        """The scale factors are derived from the lengths of the basis
        vectors of the image CS compared to the detector CS unit vectors."""
        scale_detector_u = self.detector.u.unitVector().dot(image.u)
        scale_detector_v = self.detector.v.unitVector().dot(image.v)
        scale_detector_w = self.detector.w.unitVector().dot(image.w)

        # Save a source CS as seen from the detector CS. This is convenient to
        # later get the SDD, ufoc and vfoc:
        source_from_image = copy.deepcopy(self.source)
        source_from_image.changeReferenceFrame(world, image)

        # Make the volume CS the new world CS:
        source.changeReferenceFrame(world, volume)
        image.changeReferenceFrame(world, volume)
        volume.changeReferenceFrame(world, volume)

        # Translation vector from volume to source:
        rfoc = source.center - volume.center
        xfoc = rfoc.x
        yfoc = rfoc.y
        zfoc = rfoc.z
        SOD = rfoc.length()

        # Focus point on detector: principal, perpendicular ray.
        # In the detector coordinate system, ufoc and vfoc are the u and v coordinates
        # of the source center; SDD (perpendicular to detector plane) is source w coordinate.
        ufoc = source_from_image.center.x / scale_u
        vfoc = source_from_image.center.y / scale_v
        wfoc = source_from_image.center.z / scale_w
        SDD  = abs(source_from_image.center.z)

        # Mirror volume:
        if mirror:
            M = Matrix(values=[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]])
        else:
            M = Matrix(values=[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0,  1, 0], [0, 0, 0, 1]])

        # Scale: volume units -> world units
        A = Matrix(values=[[scale_volume_u, 0, 0, 0], [0, scale_volume_v, 0, 0], [0, 0, scale_volume_w, 0], [0, 0, 0, 1]])

        # Move origin to source (the origin of the camera CS)
        F = Matrix(values=[[1, 0, 0, -xfoc], [0, 1, 0, -yfoc], [0, 0, 1, -zfoc]])

        # Rotations:
        R = basisTransformMatrix(volume, image)

        # Projection onto detector and scaling:
        S = Matrix(values=[[SDD/scale_detector_u, 0, 0], [0, SDD/scale_detector_v, 0], [0, 0, 1.0/scale_detector_w]])

        # Shift in detector CS: (ufoc and vfoc must be in scaled units)
        T = Matrix(values=[[1, 0, ufoc], [0, 1, vfoc], [0, 0, 1]])

        # Multiply all together:
        P = T * (S * (R * (F * (A*M))))

        # Renormalize:
        lower_right = P.get(col=3, row=2)
        print("lower right: {} SOD: {}".format(lower_right, SOD))
        if lower_right != 0:
            P.scale(1.0/lower_right)

        return P

    def createDetectorFlatField_rays(self):
        """ Calculate an analytical free beam intensity distribution
            picture for the given detector, to be used for an
            ideal flat field correction. """
        width      = self.detector.cols()
        height     = self.detector.rows()
        pixelSizeU = self.detector.pitchU
        pixelSizeV = self.detector.pitchV

        if(width is None):
            raise Exception("The detector width (in pixels) must be provided through a valid CTSimU JSON file.")
        if(height is None):
            raise Exception("The detector height (in pixels) must be provided through a valid CTSimU JSON file.")
        if(pixelSizeU is None):
            raise Exception("The pixel size (in mm) in u direction must be provided through a valid CTSimU JSON file.")
        if(pixelSizeV is None):
            raise Exception("The pixel size (in mm) in v direction must be provided through a valid CTSimU JSON file.")

        flatField = Image()
        flatField.shape(width, height, 0, flatField.getInternalDataType())

        # Positions of detector and source center:
        dx = self.detector.center.x
        dy = self.detector.center.y
        dz = self.detector.center.z

        sx = self.source.center.x
        sy = self.source.center.y
        sz = self.source.center.z

        # Vectors of the detector coordinate system:
        ux = self.detector.u.x
        uy = self.detector.u.y
        uz = self.detector.u.z

        vx = self.detector.v.x
        vy = self.detector.v.y
        vz = self.detector.v.z

        wx = self.detector.w.x
        wy = self.detector.w.y
        wz = self.detector.w.z


        # Angle 'alpha' between detector normal and connection line [detector center -- source]:
        connectionLine = Vector(dx-sx, dy-sy, dz-sz)

        alpha = abs(self.detector.w.angle(connectionLine))
        if alpha > (math.pi/2):
            alpha = math.pi - alpha

        # Distance between source center and detector center:
        dist = self.detector.center.distance(self.source.center)

        # Source - Detector distance (SDD) defined by shortest distance between source and detector:
        SDD = dist * math.cos(alpha)

        log("Geometry definition from JSON file:\n\
   Detector Angle:    {}\n\
   Detector Distance: {}\n\
   SDD:               {}\n\
   Pixels U:          {}\n\
   Pixels V:          {}\n\
   Pitch U:           {}\n\
   Pitch V:           {}\n\
   Source:            {}, {}, {}\n\
   Detector:          {}, {}, {}\n\
   Connection Vector: {}, {}, {}\n\
   Detector Vector U: {}, {}, {}\n\
   Detector Vector V: {}, {}, {}\n\
   Detector Vector W: {}, {}, {}".format(alpha, dist, SDD, width, height, pixelSizeU, pixelSizeV, sx, sy, sz, dx, dy, dz, connectionLine.x, connectionLine.y, connectionLine.z, ux, uy, uz, vx, vy, vz, wx, wy, wz))

        maxIntensity = 0
        maxX = 0
        maxY = 0
        minDistToSource = 0
        brightestIncidenceAngle = 0

        gridSize = 3
        gridSizeSq = gridSize*gridSize

        for x in range(width):
            for y in range(height):
                factorSum = 0
                for gx in range(gridSize):
                    for gy in range(gridSize):
                        # Calculate coordinates of pixel center in mm:
                        # Grid with margin:
                        stepSize = 1.0 / (gridSize+1)
                        pixel = self.detector.pixelVector(x+(gx+1)*stepSize, y+(gy+1)*stepSize)

                        # Grid with no margin:
                        #if gridSize > 1:
                        #    stepSize = 1.0 / (gridSize-1)
                        #    pixel = self.detector.pixelVector(x+gx*stepSize, y+gy*stepSize)
                        #else:
                        #    pixel = self.detector.pixelVectorCenter(x, y)

                        distToSource = self.source.center.distance(pixel)

                        # Angle of incident rays:
                        vecSourceToPixel = Vector(pixel.x-sx, pixel.y-sy, pixel.z-sz)
                        incidenceAngle = abs(self.detector.w.angle(vecSourceToPixel))
                        if incidenceAngle > (math.pi/2):
                            incidenceAngle = math.pi - incidenceAngle

                        intensityFactor = math.cos(incidenceAngle)*math.pow(SDD/distToSource, 2)
                        factorSum += intensityFactor

                intensityWeight = factorSum / gridSizeSq

                if intensityWeight > maxIntensity:
                    maxIntensity = intensityWeight
                    maxX = x
                    maxY = y
                    minDistToSource = distToSource
                    brightestIncidenceAngle = incidenceAngle

                flatField.setPixel(x, y, intensityWeight)

            progress = 100*(float(x+1)/float(width))
            print("\rCalculating analytical flat field... {:0.1f}%    ".format(progress), end='')

        print("\rCalculating analytical flat field... 100%  ")

        #print("Brightest Pixel: {}, {}".format(maxX, maxY))
        print("  Dist to Source: {}".format(minDistToSource))
        print("  Angle: {} rad = {} deg".format(brightestIncidenceAngle, 180*brightestIncidenceAngle/math.pi))

        return flatField

    def pixelAreaOnUnitSphere(self, A, B, C, D):
        # Source must be at (0, 0, 0) relative to given detector,
        # and A, B, C, D must be vectors pointing to pixel corners in
        # world coordinate system.

        # Define normals of circular planes, pointing into the
        # triangle or out of the triangle:
        ABin  = A.cross(B)
        BCin  = B.cross(C)
        CAin  = C.cross(A)
        ABout = ABin.inverse()
        BCout = BCin.inverse()
        CAout = CAin.inverse()
        
        ACin  = A.cross(C)
        DAin  = D.cross(A)
        CDin  = C.cross(D)
        ACout = ACin.inverse()
        DAout = DAin.inverse()
        CDout = CDin.inverse()

        # Spherical Triangle ABC:
        alpha = ABin.angle(CAout)
        beta  = BCin.angle(ABout)
        gamma = CAin.angle(BCout)

        # areaABC = alpha + beta + gamma - math.pi

        # Spherical Triangle ACD:
        rho   = ACin.angle(DAout)
        sigma = CDin.angle(ACout)
        tau   = DAin.angle(CDout)

        # areaACD = rho + tau + sigma - math.pi

        pxSphericalArea = (alpha + beta + gamma + rho + sigma + tau) - 2*math.pi

        return pxSphericalArea

    def triangleAreaOnUnitSphere(self, A, B, C):
        # Source must be at (0, 0, 0) relative to given detector,
        # and A, B, C must be vectors pointing to ´triangle corners in
        # world coordinate system.

        # Define normals of circular planes, pointing into the
        # triangle or out of the triangle:
        ABin  = A.cross(B)
        BCin  = B.cross(C)
        CAin  = C.cross(A)
        ABout = ABin.inverse()
        BCout = BCin.inverse()
        CAout = CAin.inverse()

        # Spherical Triangle ABC:
        alpha = ABin.angle(CAout)
        beta  = BCin.angle(ABout)
        gamma = CAin.angle(BCout)

        areaABC = alpha + beta + gamma - math.pi

        return areaABC

    def polygonAreaOnUnitSphere(self, polygon):
        # Source must be at (0, 0, 0) relative to given detector,
        # and A, B, C must be vectors pointing to triangle corners in
        # world coordinate system.

        if len(polygon.points) >= 3:
            # Start at first point
            p1 = polygon.points[0]

            area = 0

            for i in range(1, len(polygon.points)-1):
                p2 = polygon.points[i]
                p3 = polygon.points[i+1]

                area += self.triangleAreaOnUnitSphere(p1, p2, p3)

            return area
        else:
            return 0

    """
    def createDetectorFlatField_sphere_old(self, clippingPolygon=None):
        # Positions of detector and source center:
        dx = self.detector.center.x
        dy = self.detector.center.y
        dz = self.detector.center.z

        sx = self.source.center.x
        sy = self.source.center.y
        sz = self.source.center.z

        # Angle 'alpha' between detector normal and connection line [detector center -- source]:
        connectionLine = Vector(dx-sx, dy-sy, dz-sz)

        alpha = abs(self.detector.w.angle(connectionLine))
        if alpha > (math.pi/2):
            alpha = math.pi - alpha

        # Distance between source center and detector center:
        dist = self.detector.center.distance(self.source.center)

        # Source - Detector distance (SDD) defined by shortest distance between source and detector,
        # or distance between source and spot of highest intensity on detector.
        SDD = dist * math.cos(alpha)

        # Create a new detector in a coordinate system where source is at (0, 0, 0):
        det = copy.deepcopy(self.detector)
        translationVector = Vector(-sx, -sy, -sz)
        det.translate(translationVector)
        det.computeGeometryParameters()

        # Calculate the area of the theoretical "brightest pixel" on the unit sphere:
        hpu = 0.5*det.pitchU
        hpv = 0.5*det.pitchV
        A = Vector(SDD,  hpu,  hpv)
        B = Vector(SDD, -hpu,  hpv)
        C = Vector(SDD, -hpu, -hpv)
        D = Vector(SDD,  hpu, -hpv)
        areaOfBrightestPixel = self.pixelAreaOnUnitSphere(A, B, C, D)

        print("SDD: {}".format(SDD))
        print("Pitch: {}, {}".format(hpu, hpv))
        #print("Brightest Pixel Area: {}".format(areaOfBrightestPixel))

        flatField = Image()
        flatField.shape(det.cols(), det.rows(), 0, flatField.getInternalDataType())

        maxArea = 0
        maxX = 0
        maxY = 0
        maxCenter = 0
        # Go through pixels:
        for x in range(det.cols()):
            for y in range(det.rows()):
                # Define Pixel corners:
                A = det.pixelVector(x,   y)
                B = det.pixelVector(x+1, y)
                C = det.pixelVector(x+1, y+1)
                D = det.pixelVector(x,   y+1)

                pxSphericalArea = self.pixelAreaOnUnitSphere(A, B, C, D)
                flatField.setPixel(x, y, pxSphericalArea)

                if pxSphericalArea > maxArea:
                    maxArea = pxSphericalArea
                    maxX = x
                    maxY = y

            progress = 100*(float(x+1)/float(det.cols()))
            print("\rCalculating analytical intensity profile... {:0.1f}%    ".format(progress), end='')

        # Method #1: renormalize to area of theoretically brightest pixel:
        flatField.divide(areaOfBrightestPixel)

        # Method #2: rescale maximum of actual brightest pixel to 1.0:
        #flatField.renormalize(newMin=0, newMax=1.0, currentMin=0)

        #flatField.save("ff.tif", numpy.dtype('float32'))

        print("\rCalculating analytical intensity profile... 100%  ")


        maxCenter = det.pixelVectorCenter(maxX, maxY)
        distToSource = maxCenter.length()
        incidenceAngle = abs(self.detector.w.angle(maxCenter))

        #print("Brightest Pixel: {}, {}".format(maxX, maxY))
        print("  Vector: {}, {}, {}".format(maxCenter.x, maxCenter.y, maxCenter.z))
        print("  Distance to Source: {}".format(distToSource))
        print("  Spherical Area: {}".format(maxArea))
        print("  Angle: {} rad = {} deg".format(incidenceAngle, 180*incidenceAngle/math.pi))

        return flatField
    """

    def createDetectorFlatField_sphere(self, *coverPolygons):
        """ Calculate an analytical free beam intensity distribution
            picture for the given detector, to be used for an
            ideal flat field correction.

            Geometrical approach using spherical geometry. """

        # Change to the detector coordinate system:
        D = copy.deepcopy(self.detector)
        S = copy.deepcopy(self.source)
        world = CoordinateSystem()  # will be initialized as world

        S.changeReferenceFrame(world, D)
        D.changeReferenceFrame(world, D)
        D.computeGeometryParameters()

        # Source - Detector distance (SDD) defined by shortest distance between source and detector,
        # or distance between source and spot of highest intensity on detector.
        SDD = abs(S.center.z)

        # Calculate the area of the theoretical "brightest pixel" on the unit sphere:
        pu = D.pitchU
        pv = D.pitchV
        nRows = D.rows()
        nCols = D.cols()

        hpu = 0.5*pu
        hpv = 0.5*pv
        pA = Vector(SDD,  hpu,  hpv)
        pB = Vector(SDD, -hpu,  hpv)
        pC = Vector(SDD, -hpu, -hpv)
        pD = Vector(SDD,  hpu, -hpv)
        areaOfBrightestPixel = self.pixelAreaOnUnitSphere(pA, pB, pC, pD)

        # Full flat field image (without any clipping bodies):
        flatField = Image()
        flatField.shape(D.cols(), D.rows(), 0, flatField.getInternalDataType())

        # A second image with a clipping body under consideration: (both will be returned)
        clippedFlatField = None
        if len(coverPolygons) > 0:
            clippedFlatField = Image()
            clippedFlatField.shape(D.cols(), D.rows(), 0, flatField.getInternalDataType())

        # Upper left detector corner in world coordinates (remember: world is now the detector CS)
        p00 = D.pixelVector(0, 0)

        stepRight      = Vector(pu, 0,  0)
        stepDown       = Vector(0,  pv, 0)
        stepRightDown  = Vector(pu, pv, 0)

        # Move the clipping polygon to a coordinate system
        # where source is centered:
        for coverPolygon in coverPolygons:
            for p in range(len(coverPolygon.points)):
                coverPolygon.points[p] = coverPolygon.points[p] - S.center

        # Go through pixels:
        for x in range(nCols):
            for y in range(nRows):
                # Pixel in world coordinates (remember: world is now the detector CS)
                shift = Vector(x*pu, y*pv, 0)

                # Define Pixel corners:
                pA = p00 + shift
                pB = pA  + stepRight
                pC = pA  + stepRightDown
                pD = pA  + stepDown

                # Center source at (0, 0, 0):
                pA = pA - S.center
                pB = pB - S.center
                pC = pC - S.center
                pD = pD - S.center

                pixelPolygon = Polygon(pA, pB, pC, pD)
                pxSphericalArea  = self.polygonAreaOnUnitSphere(pixelPolygon)

                flatField.setPixel(x, y, pxSphericalArea)

                if len(coverPolygons) > 0:
                    for coverPolygon in coverPolygons:
                        pixelPolygon = pixelPolygon.clip(coverPolygon)
                        
                    # Remove the intensity covered by the clipping polygon:
                    pixelPolygon.make3D(zComponent=SDD)
                    subarea = self.polygonAreaOnUnitSphere(pixelPolygon)
                    pxSphericalArea -= subarea

                    clippedFlatField.setPixel(x, y, pxSphericalArea)

            progress = 100*(float(x+1)/float(D.cols()))
            print("\rCalculating analytical intensity profile... {:0.1f}%    ".format(progress), end='')

        # Method #1: renormalize to area of theoretically brightest pixel:
        flatField.divide(areaOfBrightestPixel)
        if clippedFlatField != None:
            clippedFlatField.divide(areaOfBrightestPixel)

        # Method #2: rescale maximum of actual brightest pixel to 1.0:
        #flatField.renormalize(newMin=0, newMax=1.0, currentMin=0)

        #flatField.save("ff.tif", numpy.dtype('float32'))

        print("\rCalculating analytical intensity profile... 100%  ")

        return flatField, clippedFlatField

    def solidAngle(self, l, m):
        """ Solid angle helper function for intensity profile. Approach by Florian Wohlgemuth. """
        if l != 0:
            return (l/abs(l)) * math.atan(abs(l)*m/math.sqrt(1.0+l**2+m**2))
        else:
            return 0

    def createDetectorFlatField_analytical(self):
        """ Calculate an analytical free beam intensity distribution
            picture for the given detector, to be used for an
            ideal flat field correction.

            Analytical approach by Florian Wohlgemuth. """

        width  = self.detector.cols()
        height = self.detector.rows()
        pitchU = self.detector.pitchU
        pitchV = self.detector.pitchV

        if(width is None):
            raise Exception("The detector width (in pixels) must be provided through a valid CTSimU JSON file.")
        if(height is None):
            raise Exception("The detector height (in pixels) must be provided through a valid CTSimU JSON file.")
        if(pitchU is None):
            raise Exception("The pixel size (in mm) in u direction must be provided through a valid CTSimU JSON file.")
        if(pitchV is None):
            raise Exception("The pixel size (in mm) in v direction must be provided through a valid CTSimU JSON file.")

        flatField = Image()
        flatField.shape(width, height, 0, flatField.getInternalDataType())

        # Positions of detector and source center:
        dx = self.detector.center.x
        dy = self.detector.center.y
        dz = self.detector.center.z

        sx = self.source.center.x
        sy = self.source.center.y
        sz = self.source.center.z

        # Vectors of the detector coordinate system:
        ux = self.detector.u.x
        uy = self.detector.u.y
        uz = self.detector.u.z

        vx = self.detector.v.x
        vy = self.detector.v.y
        vz = self.detector.v.z

        wx = self.detector.w.x
        wy = self.detector.w.y
        wz = self.detector.w.z


        # Angle 'alpha' between detector normal and connection line [detector center -- source]:
        connectionLine = Vector(dx-sx, dy-sy, dz-sz)

        alpha = abs(self.detector.w.angle(connectionLine))
        if alpha > (math.pi/2):
            alpha = math.pi - alpha

        # Distance between source center and detector center:
        dist = self.detector.center.distance(self.source.center)

        # Source - Detector distance (SDD) defined by shortest distance between source and detector:
        SDD = dist * math.cos(alpha)

        maxIntensity = 0
        maxX = 0
        maxY = 0

        # Create a new detector in a coordinate system where source is at (0, 0, 0):
        det = copy.deepcopy(self.detector)
        translationVector = Vector(-sx, -sy, -sz)
        det.translate(translationVector)
        det.computeGeometryParameters()

        upperLeft_u = det.pixelVector(0, 0).dot(self.detector.u)
        upperLeft_v = det.pixelVector(0, 0).dot(self.detector.v)
        upperLeft_w = det.pixelVector(0, 0).dot(self.detector.w)

        if upperLeft_w != 0:   # check if detector is not facing its edge towards the source
            for x in range(width):
                for y in range(height):
                    nu = x
                    nv = y
                    lambda0 = (upperLeft_u + nu*pitchU) / upperLeft_w
                    lambda1 = (upperLeft_u + (nu+1)*pitchU) / upperLeft_w
                    mu0     = (upperLeft_v + nv*pitchV) / upperLeft_w
                    mu1     = (upperLeft_v + (nv+1)*pitchV) / upperLeft_w

                    omega = self.solidAngle(lambda0, mu0) + self.solidAngle(lambda1, mu1) - self.solidAngle(lambda0, mu1) - self.solidAngle(lambda1, mu0)

                    if omega > maxIntensity:
                        maxIntensity = omega
                        maxX = x
                        maxY = y

                    flatField.setPixel(x, y, omega)

                progress = 100*(float(x+1)/float(width))
                print("\rCalculating analytical flat field... {:0.1f}%    ".format(progress), end='')

        print("\rCalculating analytical flat field... 100%  ")

        #print("Brightest Pixel: {}, {}".format(maxX, maxY))
        # print("  Dist to Source: {}".format(minDistToSource))
        # print("  Angle: {} rad = {} deg".format(brightestIncidenceAngle, 180*brightestIncidenceAngle/math.pi))

        # Method #1: find hypothetical brightest pixel
        # Calculate the area of the theoretical "brightest pixel" on the unit sphere:
        hpu = 0.5*det.pitchU
        hpv = 0.5*det.pitchV
        A = Vector(SDD,  hpu,  hpv)
        B = Vector(SDD, -hpu,  hpv)
        C = Vector(SDD, -hpu, -hpv)
        D = Vector(SDD,  hpu, -hpv)
        areaOfBrightestPixel = self.pixelAreaOnUnitSphere(A, B, C, D)
        flatField.divide(areaOfBrightestPixel)

        # Method #2: rescale actual maximum to 1.
        #flatField.renormalize(newMin=0, newMax=1.0, currentMin=0)

        #flatField.save("ff.tif", dataType="float32")
        return flatField


    def createDetectorFlatField(self):
        return createDetectorFlatField_analytical()

def writeCERAconfig(geo, totalAngle, projectionFilePattern, matrices, basename, voxelsX, voxelsY, voxelsZ, i0max=60000):
    now = datetime.now()

    nProjections = len(matrices)
    projTableString = """projtable.txt version 3
{timestring}

# format: angle / entries of projection matrices
{nProjections}
""".format(
    nProjections=nProjections,
    timestring=now.strftime("%a %b %d %H:%M:%S %Y")
    )

    for i in range(nProjections):
        m=matrices[i]
        projTableString += """@{nr}
0.0 0.0
{c00} {c01} {c02} {c03}
{c10} {c11} {c12} {c13}
{c20} {c21} {c22} {c23}

""".format(
        nr=(i+1),
        c00=m.get(col=0, row=0),
        c01=m.get(col=1, row=0),
        c02=m.get(col=2, row=0),
        c03=m.get(col=3, row=0),
        c10=m.get(col=0, row=1),
        c11=m.get(col=1, row=1),
        c12=m.get(col=2, row=1),
        c13=m.get(col=3, row=1),
        c20=m.get(col=0, row=2),
        c21=m.get(col=1, row=2),
        c22=m.get(col=2, row=2),
        c23=m.get(col=3, row=2)
    )

    with open("{}_projtable.txt".format(basename), 'w') as f:
        f.write(projTableString)
        f.close()

    voxelSizeXY = geo.detector.pitchU * geo.SOD / geo.SDD
    voxelSizeZ  = geo.detector.pitchV * geo.SOD / geo.SDD

    configFileString = """#CERACONFIG

[Projections]
NumChannelsPerRow = {nCols}
NumRows = {nRows}
PixelSizeU = {psu}
PixelSizeV = {psv}
Rotation = None
FlipU = false
FlipV = true
Padding = 0
BigEndian = false
CropBorderRight = 0
CropBorderLeft = 0
CropBorderTop = 0
CropBorderBottom = 0
BinningFactor = None
SkipProjectionInterval = 1
ProjectionDataDomain = Intensity
RawHeaderSize = 0

[Volume]
SizeX = {volx}
SizeY = {voly}
SizeZ = {volz}
# Midpoints are only necessary for reconstructions
# without projection matrices.
MidpointX = {midpointx}
MidpointY = {midpointy}
MidpointZ = {midpointz}
VoxelSizeX = {vsx}
VoxelSizeY = {vsy}
VoxelSizeZ = {vsz}
Datatype = float

[CustomKeys]
NumProjections = {nProjections}    
ProjectionFileType = tiff
VolumeOutputPath = {basename}.raw
ProjectionStartNum = 0
ProjectionFilenameMask = {projFilePattern}

[CustomKeys.ProjectionMatrices]
SourceObjectDistance = {SOD}
SourceImageDistance = {SDD}
DetectorOffsetU = {offu}
DetectorOffsetV = {offv}
StartAngle = {startAngle}
ScanAngle = {scanAngle}
AquisitionDirection = CW
a = {a}
b = {b}
c = {c}
ProjectionMatrixFilename = {basename}_projtable.txt

[Backprojection]
ClearOutOfRegionVoxels = false
InterpolationMode = bilinear
FloatingPointPrecision = half
Enabled = true

[Filtering]
Enabled = true
Kernel = shepp

[I0Log]
Enabled = true
Epsilon = 1.0E-5
GlobalI0Value = {i0max}
""".format(
    basename=basename,
    nCols=int(geo.detector.cols()),
    nRows=int(geo.detector.rows()),
    psu=geo.detector.pitchU,
    psv=geo.detector.pitchV,
    volx=int(voxelsX),
    voly=int(voxelsY),
    volz=int(voxelsZ),
    midpointx=0, #!
    midpointy=0, #!
    midpointz=0, #!
    vsx=voxelSizeXY,
    vsy=voxelSizeXY,
    vsz=voxelSizeZ,
    nProjections=int(nProjections),
    projFilePattern=projectionFilePattern,
    SOD=geo.SOD,
    SDD=geo.SDD,
    offu=0, #!
    offv=0, #!
    startAngle=0, #!
    scanAngle=totalAngle,
    a=0, #!
    b=0, #!
    c=0, #!
    i0max=i0max
    )

    with open("{}.config".format(basename), 'w') as f:
        f.write(configFileString)
        f.close()


def writeOpenCTFile(geo, totalAngle, boundingBoxX, boundingBoxY, boundingBoxZ, matrices, filename, volumename, projectionFilenames):
    nProjections = len(matrices)
    matrixString = ""

    i = 0
    for m in matrices:
        if i>0:
            matrixString += ",\n\n            "

        matrixString += "[ "
        for row in range(m.rows):
            if row > 0:
                matrixString += ",\n              "
            matrixString += "["
            for col in range(m.cols):
                if col > 0:
                    matrixString += ", "

                matrixString += "{}".format(m.get(col=col, row=row))
            matrixString += "]"
        matrixString += " ]"
        i += 1

    filesString = ""
    if len(matrices) == len(projectionFilenames):
        for i in range(len(matrices)):
            if i > 0:
                filesString += ",\n                "
            #filesString += '"{:05d}.tif"'.format(i)
            filesString += '"{}"'.format(projectionFilenames[i])
    else:
        raise Exception("The number of projection matrices ({}) does not match the number of projection file names ({}).".format(len(matrices), len(projectionFilenames)))


    content = """{{
    "version": {{"major":1, "minor":0}},
    "openCTJSON":     {{
        "versionMajor": 1,
        "versionMinor": 0,
        "revisionNumber": 0,
        "variant": "FreeTrajectoryCBCTScan"
    }},
    "units": {{
        "length": "Millimeter"
    }},
    "volumeName":  "{volumeName}",
    "projections": {{
        "numProjections":  {nProjections},
        "intensityDomain": true,
        "images":          {{
            "directory": ".",
            "dataType":  "UInt16",
            "fileType":  "TIFF",
            "files":     [
                {filesString}
            ]
        }},
        "matrices": [
            {matrixString}
        ]
    }},
    "geometry":    {{
        "totalAngle":           {totalAngle},
        "skipAngle":            0,
        "detectorPixel":        [
            {nPixelsX},
            {nPixelsY}
        ],
        "detectorSize":         [
            {detectorSizeX},
            {detectorSizeY}
        ],
        "mirrorDetectorAxis":   "",
        "distanceSourceObject": {SOD},
        "distanceObjectDetector": {ODD},
        "objectBoundingBox":    [
            [
                {bbx},
                0.0,
                0.0,
                0.0
            ],
            [
                0.0,
                {bby},
                0.0,
                0.0
            ],
            [
                0.0,
                0.0,
                {bbz},
                0.0
            ],
            [
                0.0,
                0.0,
                0.0,
                1.0
            ]
        ]
    }},
    "corrections": {{
        "brightImages": {{
            "directory": "",
            "dataType":  "",
            "fileType":  "",
            "files":     []
        }},
        "darkImage":    {{
            "file":     "",
            "dataType": "",
            "fileType": ""
        }},
        "badPixelMask": {{
            "file":     "",
            "dataType": "",
            "fileType": ""
        }},
        "intensities":  []
    }}
}}""".format(
    nProjections=nProjections,
    matrixString=matrixString,
    nPixelsX=int(geo.detector.cols()),
    nPixelsY=int(geo.detector.rows()),
    detectorSizeX=geo.detector.physWidth,
    detectorSizeY=geo.detector.physHeight,
    SOD=geo.SOD,
    ODD=geo.ODD,
    totalAngle=totalAngle,
    bbx=boundingBoxX,
    bby=boundingBoxY,
    bbz=boundingBoxZ,
    filesString=filesString,
    volumeName=volumename
    )

    with open(filename, 'w') as f:
        f.write(content)
        f.close()