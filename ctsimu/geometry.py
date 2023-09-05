# -*- coding: UTF-8 -*-
"""
Coordinate systems, transformations and projection matrix functionality.

.. include:: ./geometry.md
"""

import numpy
import os    # File and path handling
import json
import math
import copy
import pkgutil
import warnings
from datetime import datetime

from .primitives import *
from .image import Image  # To create detector flat field

class CoordinateSystem:
    """Coordinate system: center point and axis vectors.

    The center and axis vectors are expressed in terms of the
    object's reference coordinate system, which must be known implicitly
    when objects of this class are used.

    Attributes
    ----------
    center : ctsimu.primitives.Vector
        The location of the center point in a reference
        coordinate system (usually world or stage).

    u : ctsimu.primitives.Vector
        Basis vector for the u axis.

    v : ctsimu.primitives.Vector
        Basis vector for the v axis.

    w : ctsimu.primitives.Vector
        Basis vector for the w axis.
    """

    def __init__(self):
        """Initialized as a standard world coordinate system."""
        self.center = Vector(0, 0, 0)
        self.u      = Vector(1, 0, 0)
        self.v      = Vector(0, 1, 0)
        self.w      = Vector(0, 0, 1)

    def __str__(self):
        """Information string for easy printing."""
        txt  = "Center: {}\n".format(self.center)
        txt += "u:      {}\n".format(self.u)
        txt += "v:      {}\n".format(self.v)
        txt += "w:      {}\n".format(self.w)
        return txt

    def reset(self):
        """Reset to a standard world coordinate system."""
        self.center = Vector(0, 0, 0)
        self.u      = Vector(1, 0, 0)
        self.v      = Vector(0, 1, 0)
        self.w      = Vector(0, 0, 1)

    def make_unit_coordinate_system(self):
        """Convert all basis vectors into unit vectors."""
        self.u.make_unit_vector()
        self.v.make_unit_vector()
        self.w.make_unit_vector()
        self.update()

    def make_from_vectors(self, center:'Vector', u:'Vector', w:'Vector'):
        """Create a right-handed coordinate system from the `center`,
        `u` vector (first basis vector) and `w` vector (third basis vector).
        The vector `v` will be determined from the cross product `w`×`u`.

        Parameters
        ----------
        center : Vector
            Object's center point in reference coordinate system,
            origin of local {u,v,w} coordinate system.

        u : Vector
            Basis vector u in terms of reference coordinate system.

        w : Vector
            Basis vector w in terms of reference coordinate system.

        Notes
        -----
        Basis vectors must be orthogonal.
        """
        self.center  = center
        self.u = u
        self.w = w
        self.v = self.w.cross(self.u)

        self.update()

    def make(self, cx:float, cy:float, cz:float, ux:float, uy:float, uz:float, wx:float, wy:float, wz:float):
        """Set up the coordinate system from vector components (all floats)
        for the center (`cx`, `cy`, `cz`), the `u` vector (first basis vector,
        `ux`, `uy`, `uz`) and the `w` vector (third basis vector, `wx`, `wy`, `wz`).

        Parameters
        ----------
        cx : float
            Center x coordinate.

        cy : float
            Center y coordinate.

        cz : float
            Center z coordinate.

        ux : float
            `u` vector x component.

        uy : float
            `u` vector y component.

        uz : float
            `u` vector z component.

        wx : float
            `w` vector x component.

        wy : float
            `w` vector y component.

        wz : float
            `w` vector z component.
        """
        self.center = Vector(cx, cy, cz)
        self.u      = Vector(ux, uy, uz)
        self.w      = Vector(wx, wy, wz)
        self.v      = self.w.cross(self.u)

        self.update()

    def set_u_w(self, u:'Vector', w:'Vector'):
        """Set u and w vector, calculate v from cross product (right-handed).

        Parameters
        ----------
        u : ctsimu.primitives.Vector
            Basis vector for u direction.

        v : ctsimu.primitives.Vector
            Basis vector for v direction.
        """
        self.u = u
        self.w = w
        self.v = w.cross(u)

    def get_copy(self) -> 'CoordinateSystem':
        """Get a copy of this coordinate system.

        Returns
        -------
        copy_cs : CoordinateSystem
            Copy of this coordinate system.
        """
        new_cs = CoordinateSystem()
        new_cs.center = Vector(self.center.x(), self.center.y(), self.center.z())
        new_cs.u      = Vector(self.u.x(), self.u.y(), self.u.z())
        new_cs.v      = Vector(self.v.x(), self.v.y(), self.v.z())
        new_cs.w      = Vector(self.w.x(), self.w.y(), self.w.z())
        return new_cs

    def copy_cs(self, other:'CoordinateSystem'):
        """Make this CoordinateSystem a copy of the `other` coordinate system.

        Parameters
        ----------
        other : CoordinateSystem
            Another coordinate system to copy.
        """
        self.center = other.center.get_copy()
        self.u      = other.u.get_copy()
        self.v      = other.v.get_copy()
        self.w      = other.w.get_copy()

    def update(self):
        """Signal a manual update to the center position or orientation vectors."""
        self.center.update()
        self.u.update()
        self.v.update()
        self.w.update()

    def translate(self, translation_vector:'Vector'):
        """Shift center by given translation vector.

        Parameters
        ----------
        translation_vector : Vector
            Vector by which the object's center point should be shifted.
            Its components are added to the center's components.
        """
        self.center.add(translation_vector)

    def translate_in_direction(self, direction:'Vector', distance:float):
        """Shift center in given `direction` by given `distance`.

        Parameters
        ----------
        direction : Vector
            Vector along which the center point should be shifted.
            It must not be a unit vector.

        distance : float
            Distance by which the center point will travel
        """
        t = direction.unit_vector().scaled(factor=distance)
        self.translate(translation_vector=t)

    def translate_x(self, dx:float):
        """Translate coordinate system in x direction of reference
        coordinate system by distance `dx`.

        Parameters
        ----------
        dx : float
            Shift amount in x direction.
        """
        self.center.set_x(self.center.x() + float(dx))

    def translate_y(self, dy: float):
        """Translate coordinate system in y direction of reference
        coordinate system by distance `dy`.

        Parameters
        ----------
        dy : float
            Shift amount in y direction.
        """
        self.center.set_y(self.center.y() + float(dy))

    def translate_z(self, dz: float):
        """Translate coordinate system in z direction of reference
        coordinate system by distance `dz`.

        Parameters
        ----------
        dz : float
            Shift amount in z direction.
        """
        self.center.set_z(self.center.z() + float(dz))

    def translate_u(self, du:float):
        """Translate coordinate system in u direction by distance `du`.

        Parameters
        ----------
        du : float
            Shift amount in u direction.
        """
        self.translate_in_direction(direction=self.u, distance=du)

    def translate_v(self, dv:float):
        """Translate coordinate system in v direction by distance `dv`.

        Parameters
        ----------
        dv : float
            Shift amount in v direction.
        """
        self.translate_in_direction(direction=self.v, distance=dv)

    def translate_w(self, dw:float):
        """Translate coordinate system in w direction by distance `dw`.

        Parameters
        ----------
        dw : float
            Shift amount in w direction.
        """
        self.translate_in_direction(direction=self.w, distance=dw)

    def rotate(self, axis:'Vector', angle:float):
        """Rotate coordinate system around a given axis by the given angle (in rad).

        This does not move the center point, as the axis vector is assumed
        to be attached to the center of the coordinate system.

        Parameters
        ----------
        axis : Vector
            The axis of rotation, in terms of the object's
            reference coordinate system.

        angle : float
            Rotation angle (in rad), mathematically positive direction (right-hand rule).
        """
        R = rotation_matrix(axis, angle)
        self.u.transform(R)
        self.v.transform(R)
        self.w.transform(R)

    def rotate_around_pivot_point(self, axis:'Vector', angle:float, pivot:'Vector'):
        """Rotate coordinate system around a pivot point.
        Generally, this will result in a different center position,
        as the axis of rotation is assumed to be attached to the
        pivot point instead of the center of the coordinate system.

        Parameters
        ----------
        axis : Vector
            Rotation axis, in terms of the object's reference coordinate system.

        angle : float
            Rotation angle (in rad).

        pivot : Vector
            Pivot point, in terms of the object's reference coordinate system.
        """

        # Move coordinate system such that pivot point is at world origin:
        self.center.subtract(pivot)

        # Rotate center point and transform back into
        # world coordinate system:
        self.center.rotate(axis, angle)
        self.center.add(pivot)

        # Rotate the coordinate system itself:
        self.rotate(axis, angle)

    def rotate_around_x(self, angle: float):
        """Rotate object around x axis of its reference coordinate system
        by given angle (in rad).

        Parameters
        ----------
        angle : float
            Rotation angle in rad, mathematically positive direction (right-hand rule).
        """
        if angle != 0:
            x_axis = Vector(1, 0, 0)
            self.rotate(x_axis, angle)

    def rotate_around_y(self, angle: float):
        """Rotate object around y axis of its reference coordinate system
        by given angle (in rad).

        Parameters
        ----------
        angle : float
            Rotation angle in rad, mathematically positive direction (right-hand rule).
        """
        if angle != 0:
            y_axis = Vector(0, 1, 0)
            self.rotate(y_axis, angle)

    def rotate_around_z(self, angle: float):
        """Rotate object around z axis of its reference coordinate system
        by given angle (in rad).

        Parameters
        ----------
        angle : float
            Rotation angle in rad, mathematically positive direction (right-hand rule).
        """
        if angle != 0:
            z_axis = Vector(0, 0, 1)
            self.rotate(z_axis, angle)

    def rotate_around_u(self, angle: float):
        """Rotate object around its u axis by given angle (in rad).

        Parameters
        ----------
        angle : float
            Rotation angle in rad, mathematically positive direction (right-hand rule).
        """
        self.v.rotate(self.u, angle)
        self.w.rotate(self.u, angle)

    def rotate_around_v(self, angle: float):
        """Rotate object around its v axis by given angle (in rad).

        Parameters
        ----------
        angle : float
            Rotation angle in rad, mathematically positive direction (right-hand rule).
        """
        self.u.rotate(self.v, angle)
        self.w.rotate(self.v, angle)

    def rotate_around_w(self, angle: float):
        """Rotate object around its w axis by given angle (in rad).

        Parameters
        ----------
        angle : float
            Rotation angle in rad, mathematically positive direction (right-hand rule).
        """
        self.u.rotate(self.w, angle)
        self.v.rotate(self.w, angle)

    def transform(self, cs_from:'CoordinateSystem', cs_to:'CoordinateSystem'):
        """Relative transformation in world coordinates
        from `cs_from` to `cs_to`, result will be in world coordinates.

        Assuming this CS, `cs_from` and `cs_to`
        all three are independent coordinate systems in a common
        reference coordinate system (e.g. world). This function
        will calculate the necessary translation and rotation that
        would have to be done to superimpose `cs_from` with `cs_to`.
        This translation and rotation will, however, be applied
        to this CS, not to `cs_from`.

        Parameters
        ----------
        cs_from : CoordinateSystem
            Coordinate system before the transformation.

        cs_to : CoordinateSystem
            Coordinate system after the transformation.
        """

        # -- TRANSLATION:
        t = cs_from.center.to(cs_to.center)
        self.translate(t)

        # We need a copy of cs_from and cs_to because later on,
        # we might have to transform them and don't want to
        # affect the original cs_from passed to this function.
        # Also, cs_from or cs_to could simply be pointers to
        # this coordinate system.
        cs_fromCopy = cs_from.get_copy()
        cs_toCopy   = cs_to.get_copy()

        # -- ROTATIONS
        # Rotation to bring w axis from -> to
        wFrom = cs_fromCopy.w
        wTo   = cs_toCopy.w
        rotationAxis = wFrom.cross(wTo)

        if rotationAxis.length() == 0:
            if wTo.dot(wFrom) < 0:
                # 180° flip; vectors point in opposite direction.
                # Rotation axis is another CS basis vector.
                rotationAxis = cs_fromCopy.u.get_copy()
            else:
                # wFrom already points in direction of wTo.
                pass

        if rotationAxis.length() > 0:
            rotationAngle = wFrom.angle(wTo)
            if rotationAngle != 0:
                self.rotate_around_pivot_point(rotationAxis, rotationAngle, cs_toCopy.center)

                # Also rotate `cs_from` to make calculation of
                # rotation around u axis possible (next step):
                cs_fromCopy.rotate(rotationAxis, rotationAngle)

        # Rotation to bring u axis from -> to (around now fixed w axis)
        uFrom = cs_fromCopy.u
        uTo   = cs_toCopy.u
        rotationAxis = uFrom.cross(uTo)

        if rotationAxis.length() == 0:
            if uTo.dot(uFrom) < 0:
                # 180° flip; vectors point in opposite direction.
                # Rotation axis is another CS basis vector.
                rotationAxis = cs_fromCopy.w.get_copy()
            else:
                # uFrom already points in direction of uTo.
                pass

        if rotationAxis.length() > 0:
            rotationAngle = uFrom.angle(uTo)
            if rotationAngle != 0:
                self.rotate_around_pivot_point(rotationAxis, rotationAngle, cs_toCopy.center)

    def change_reference_frame(self, cs_from:'CoordinateSystem', cs_to:'CoordinateSystem'):
        """Change the object's reference coordinate system.

        Parameters
        ----------
        cs_from : CoordinateSystem
            Current reference coordinate system.

        cs_to : CoordinateSystem
            New reference coordinate system.

        Notes
        -----
        Both `cs_from` and `cs_to` must be in the same reference coordinate system
        (e.g., the world coordinate system).
        """

        # Rotate basis vectors into cs_to:
        T = basis_transform_matrix(cs_from, cs_to)
        self.u.transform(T)
        self.v.transform(T)
        self.w.transform(T)

        self.center = change_reference_frame_of_point(self.center, cs_from, cs_to)

ctsimu_world = CoordinateSystem()

def basis_transform_matrix(cs_from:'CoordinateSystem', cs_to:'CoordinateSystem') -> 'Matrix':
    """A matrix that transforms coordinates from `cs_from` to `cs_to`.

    `cs_from` and `cs_to` must have the same common reference frame
    (e.g. the world coordinate system). A shift in origins is not
    taken into account, i.e., their origins are assumed to be at
    the same position.

    Parameters
    ----------
    cs_from : CoordinateSystem
        The origin coordinate system.

    cs_to : CoordinateSystem
        The target coordinate system.

    Returns
    -------
    T : Matrix
        The 3x3 basis transformation matrix.

    References
    ----------
    * S. Widnall: [Lecture L3 - Vectors, Matrices and Coordinate Transformations]
    [Lecture L3 - Vectors, Matrices and Coordinate Transformations]: https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec03/
    """
    T = Matrix(3, 3)

    # Row 1:
    T.value[0][0] = cs_to.u.unit_vector().dot(cs_from.u.unit_vector())
    T.value[0][1] = cs_to.u.unit_vector().dot(cs_from.v.unit_vector())
    T.value[0][2] = cs_to.u.unit_vector().dot(cs_from.w.unit_vector())

    # Row 2:
    T.value[1][0] = cs_to.v.unit_vector().dot(cs_from.u.unit_vector())
    T.value[1][1] = cs_to.v.unit_vector().dot(cs_from.v.unit_vector())
    T.value[1][2] = cs_to.v.unit_vector().dot(cs_from.w.unit_vector())

    # Row 3:
    T.value[2][0] = cs_to.w.unit_vector().dot(cs_from.u.unit_vector())
    T.value[2][1] = cs_to.w.unit_vector().dot(cs_from.v.unit_vector())
    T.value[2][2] = cs_to.w.unit_vector().dot(cs_from.w.unit_vector())

    return T

def change_reference_frame_of_direction(direction:'Vector', cs_from:'CoordinateSystem', cs_to:'CoordinateSystem') -> 'Vector':
    """For a `direction` in `cs_from`, get the new direction in terms of `cs_to`.
    `cs_from` and `cs_to` must be in the same reference coordinate system.

    Parameters
    ----------
    direction : Vector
        Direction in terms of `cs_from`.

    cs_from : CoordinateSystem
        The original coordinate system.

    cs_to : CoordinateSystem
        The target coordinate system, in which the direction should be expressed.

    Returns
    -------
    direction_in_cs_to : Vector
        The direction in terms of cs_to.
    """

    # Rotation matrix to rotate base vectors into cs_to:
    R = basis_transform_matrix(cs_from, cs_to)

    # Perform rotation:
    return (R * direction)

def change_reference_frame_of_point(point:'Vector', cs_from:'CoordinateSystem', cs_to:'CoordinateSystem') -> 'Vector':
    """For a `point` coordinate in `cs_from`, get the new coordinate in terms of
    `cs_to`. `cs_from` and `cs_to` must be in the same reference coordinate system.

    Parameters
    ----------
    point : Vector
        Point coordinates in terms of `cs_from`.

    cs_from : CoordinateSystem
        The original coordinate system.

    cs_to : CoordinateSystem
        The target coordinate system, in which the point coordinates
        should be expressed.

    Returns
    -------
    point_in_cs_to : Vector
        The point coordinates in terms of cs_to.
    """

    # Place the point in the common reference coordinate system
    # (mathematically, this is always the 'world'):
    point_in_to = point.get_copy()
    R_to_world = basis_transform_matrix(cs_from, ctsimu_world)
    point_in_to.transform(R_to_world)
    point_in_to.add(cs_from.center)

    # Move point to the target coordinate system:
    point_in_to.subtract(cs_to.center)
    R_to_to = basis_transform_matrix(ctsimu_world, cs_to)
    point_in_to.transform(R_to_to)

    return point_in_to


class DetectorGeometry(CoordinateSystem):
    """Detector as geometrical object.

    With additional attributes for the spatial extension and
    the pixel coordinate system.

    Attributes
    ----------
    pixels_u : int
        Number of pixels in u direction.

    pixels_v : int
        Number of pixels in v direction.

    pitch_u : float
        Size of a pixel in u direction.
        In units of the reference coordinate system.

    pitch_v : float
        Size of a pixel in v direction.
        In units of the reference coordinate system.

    phys_width : float
        Physical size in u direction.
        In units of the reference coordinate system.
        Computed automatically after calling `set_size()`.

    phys_height : float
        Physical size in v direction.
        In units of the reference coordinate system.
        Computed automatically after calling `set_size()`.

    pixel_origin : Vector
        Origin of the pixel coordinate system in terms of the reference
        coordinate system. This is the outermost corner of the
        (0,0) pixel of the detector (often the "upper left" corner).
        Computed automatically after calling `set_size()`.

    Notes
    -----
    Use `set_size()` to set the size of the detector, given its number of pixels
    and the pitch. This function automatically computes the physical dimensions
    `phys_width` and `phys_height` and the origin of the pixel coordinate system.
    """

    def __init__(self):
        """Initialize as a standard CoordinateSystem.

        Orientation, position and size must be set up manually afterwards.
        """

        # Call init from parent class:
        CoordinateSystem.__init__(self)

        self.pixels_u     = None  # Detector pixels in u direction
        self.pixels_v     = None  # Detector pixels in v direction
        self.pitch_u      = None  # Size of a pixel in u direction in units of reference coordinate system
        self.pitch_v      = None  # Size of a pixel in v direction in units of reference coordinate system
        self.phys_width   = 0     # Physical width in units of reference coordinate system
        self.phys_height  = 0     # Physical height in units of reference coordinate system

        self.pixel_origin = Vector()  # origin of pixel coordinate system in terms of reference coordinate system

    def get_copy(self):
        new_detector = DetectorGeometry()

        new_detector.center = Vector(self.center.x(), self.center.y(), self.center.z())
        new_detector.u      = Vector(self.u.x(), self.u.y(), self.u.z())
        new_detector.v      = Vector(self.v.x(), self.v.y(), self.v.z())
        new_detector.w      = Vector(self.w.x(), self.w.y(), self.w.z())

        new_detector.pixels_u     = self.pixels_u
        new_detector.pixels_v     = self.pixels_v
        new_detector.pitch_u      = self.pitch_u
        new_detector.pitch_v      = self.pitch_v
        new_detector.phys_width   = self.phys_width
        new_detector.phys_height  = self.phys_height
        new_detector.pixel_origin = self.pixel_origin.get_copy()

        return new_detector

    def size_is_set(self):
        if (self.pixels_u is None) or (self.pixels_v is None) or (self.pitch_u is None) or (self.pitch_v is None):
            return False

        return True

    def set_size(self, pixels_u:int = None, pixels_v:int = None, pitch_u:float = None, pitch_v:float = None):
        """Set the physical size of the detector.

        From the given parameters (number of pixels and pitch), the physical
        size of the detector and the position of the origin of the pixel
        coordinate system will be calculated. Make sure that the orientation
        vectors and position of the detector are correct before calling
        `set_size()`, or call `compute_geometry_parameters()` if you update
        the detector orientation or position later on.

        Parameters
        ----------
        pixels_u : int
            Number of pixels in u direction.

        pixels_v : int
            Number of pixels in v direction.

        pitch_u : float
            Pixel pitch in u direction.

        pitch_v : float
            Pixel pitch in v direction.
        """

        self.pixels_u = int(pixels_u)
        self.pixels_v = int(pixels_v)
        self.pitch_u = float(pitch_u)
        self.pitch_v = float(pitch_v)

        self.compute_geometry_parameters()

    def compute_geometry_parameters(self):
        """Calculate the physical width and height, and the position of the
        pixel coordinate system origin.

        These calculations assume that the size, position and
        orientation of the detector are correctly set up.

        Results are assigned to their member variables (attributes).
        """

        if self.size_is_set():
            # Physical width and height:
            self.phys_width  = self.pixels_u * self.pitch_u
            self.phys_height = self.pixels_v * self.pitch_v

            # Vectors of the detector coordinate system:
            ux = self.u.unit_vector().x()
            uy = self.u.unit_vector().y()
            uz = self.u.unit_vector().z()
            vx = self.v.unit_vector().x()
            vy = self.v.unit_vector().y()
            vz = self.v.unit_vector().z()

            # World coordinates of origin (0,0) of detector's pixel coordinate system:
            self.pixel_origin.set_x(self.center.x() - 0.5*(ux*self.phys_width + vx*self.phys_height))
            self.pixel_origin.set_y(self.center.y() - 0.5*(uy*self.phys_width + vy*self.phys_height))
            self.pixel_origin.set_z(self.center.z() - 0.5*(uz*self.phys_width + vz*self.phys_height))

    def cols(self) -> int:
        """Returns the number of detector columns (i.e., pixels in u direction).

        Returns
        -------
        pixels_u : int
            Number of detector columns (i.e., pixels in u direction).
        """
        return self.pixels_u

    def rows(self) -> int:
        """Returns the number of detector rows (i.e., pixels in v direction).

        Returns
        -------
        pixels_v : int
            Number of detector rows (i.e., pixels in v direction).
        """
        return self.pixels_v

    def pixel_vector(self, x: float, y: float) -> Vector:
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
        `pixel_vector_center()` may be used.

        Parameters
        ----------
        x : float
            x position in pixel coordinate system.

        y : float
            y position in pixel coordinate system.

        Returns
        -------
        pixel_vector : Vector
            Pixel position in reference coordinate system (usually world)
            as a 3D vector.
        """

        # x, y are coordinates in pixel coordinates system
        px = self.pixel_origin.x() + self.u.x()*x*self.pitch_u + self.v.x()*y*self.pitch_v
        py = self.pixel_origin.y() + self.u.y()*x*self.pitch_u + self.v.y()*y*self.pitch_v
        pz = self.pixel_origin.z() + self.u.z()*x*self.pitch_u + self.v.z()*y*self.pitch_v
        pixel_vector = Vector(px, py, pz)
        return pixel_vector

    def pixel_vector_center(self, x: float, y: float) -> Vector:
        """World position vector of pixel center, for a pixel given in integer coordinates.

        Parameters
        ----------
        x : float
            Integer x coordinate, specifies a pixel in the pixel coordinate system.

        y : float
            Integer y coordinate, specifies a pixel in the pixel coordinate system.

        Returns
        -------
        pixel_vector : Vector
            Position of the pixel center in the reference coordinate system
            (usually world) as a 3D vector.

        Notes
        -----
        If `float` coordinates are passed (non-integer),
        they are converted to integers using `math.floor`.
        """
        return self.pixel_vector(float(math.floor(x))+0.5, float(math.floor(y))+0.5)


class Geometry:
    """Geometry information about the complete CT setup.

    Keeps the source, stage and detector in one bundle and provides methods
    to calculate geometry parameters and projection matrices.

    Attributes
    ----------
    detector : DetectorGeometry
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

    brightest_spot_world : Vector
        Location of the intensity maximum on the detector, in world coordinates.
        Assuming an isotropically radiating source.
        Calculated automatically by `update()`.

    brightest_spot_detector : Vector
        Location of the intensity maximum on the detector, in terms of
        detector coordinate system. Assuming an isotropically radiating source.
        Calculated automatically by `update()`.
    """

    def __init__(self):
        self.detector    = DetectorGeometry()
        self.source      = CoordinateSystem()
        self.stage       = CoordinateSystem()

        # Backup geometry after calling store():
        self._detector_stored = None
        self._source_stored   = None
        self._stage_stored    = None

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
        self.brightest_spot_world = None
        self.brightest_spot_detector = None

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

        brightest_spot_world: Location of the intensity maximum on the detector,
            in world coordinates.  Assuming an isotropically radiating source.

        brightest_spot_detector: Location of the intensity maximum on the
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

        source_from_image.change_reference_frame(world, self.detector)
        stage_from_detector.change_reference_frame(world, self.detector)

        self.SDD = abs(source_from_image.center.z())
        self.ODD = abs(stage_from_detector.center.z())
        self.SOD = self.source.center.distance(self.stage.center)

        ## Brightest Spot in World Coordinate System:
        self.brightest_spot_world = copy.deepcopy(self.detector.w)
        self.brightest_spot_world.scale(self.SDD)
        self.brightest_spot_world.add(self.source.center)

        ## Brightest Spot in Detector Coordinate System:
        self.brightest_spot_detector = copy.deepcopy(self.brightest_spot_world)
        self.brightest_spot_detector.subtract(self.detector.center)

        pxU = self.brightest_spot_detector.dot(self.detector.u) / self.detector.pitch_u + self.detector.cols()/2.0
        pxV = self.brightest_spot_detector.dot(self.detector.v) / self.detector.pitch_v + self.detector.rows()/2.0

        self.brightest_spot_detector = Vector(pxU, pxV, 0)

        self.detector.compute_geometry_parameters()

    def store(self):
        """Store the current configuration in a backup buffer.

        The primary purpose of this function is to create a backup
        of the initial configuration, which can then always be recovered
        by a call of `Geometry.restore()`. This allows the simulation of a
        parameterized scan trajectory where each step's (or frame's)
        configuration is deterministically calculated from the initial state,
        rather than using incremental changes which could lead to the
        accumulation of rounding inaccuracies.
        """
        self._source_stored   = copy.deepcopy(self.source)
        self._detector_stored = copy.deepcopy(self.detector)
        self._stage_stored    = copy.deepcopy(self.stage)

    def restore(self):
        """Restore the configuration that has been saved by `Geometry.store()`."""

        if self._source_stored is not None:
            self.source = copy.deepcopy(self._source_stored)

        if self._detector_stored is not None:
            self.detector = copy.deepcopy(self._detector_stored)

        if self._stage_stored is not None:
            self.stage = copy.deepcopy(self._stage_stored)

        self.update()

    def info(self) -> str:
        """Generate an information string about the current geometry.

        Returns
        -------
        txt : string
            Information string for humans.
        """

        self.update()

        txt  = "Detector\n"
        txt += "========\n"
        txt += "Center:          {}\n".format(self.detector.center)
        txt += "u:               {}\n".format(self.detector.u)
        txt += "v:               {}\n".format(self.detector.v)
        txt += "w:               {}\n".format(self.detector.w)
        txt += "Pixels:          {cols} x {rows}\n".format(cols=self.detector.cols(), rows=self.detector.rows())
        txt += "Pitch:           {pitch_u} x {pitch_v}\n".format(pitch_u=self.detector.pitch_u, pitch_v=self.detector.pitch_v)
        txt += "Physical Size:   {width} x {height}\n".format(width=self.detector.phys_width, height=self.detector.phys_height)

        txt += "Brightest Spot:\n"
        txt += "  World:         {}\n".format(self.brightest_spot_world)
        txt += "  Pixels:        {}\n".format(self.brightest_spot_detector)

        txt += "\n"
        txt += "Source\n"
        txt += "======\n"
        txt += "Center:          {}\n".format(self.source.center)
        txt += "u:               {}\n".format(self.source.u)
        txt += "v:               {}\n".format(self.source.v)
        txt += "w:               {}\n".format(self.source.w)

        txt += "\n"
        txt += "Stage\n"
        txt += "=====\n"
        txt += "Center:          {}\n".format(self.stage.center)
        txt += "u:               {}\n".format(self.stage.u)
        txt += "v:               {}\n".format(self.stage.v)
        txt += "w:               {}\n".format(self.stage.w)

        txt += "\n"
        txt += "Geometry Parameters\n"
        txt += "===================\n"
        # Source - Detector distance (SDD) defined by shortest distance between source and detector:
        txt += "SDD:             {}\n".format(self.SDD)
        txt += "ODD:             {}\n".format(self.ODD)
        txt += "SOD:             {}\n".format(self.SOD)

        return txt

    def get_CERA_standard_circular_parameters(self, start_angle:float=0) -> dict:
        """Calculate all parameters for an ideal circular trajectory
        reconstruction in CERA without projection matrices.

        These can be added to the reconstruction config file for CERA.

        Parameters
        ----------
        start_angle : float
            Reconstruction start angle (in degrees). The start angle can be
            tuned to change the in-plane rotation of the reconstruction images.
            Depending on the current stage rotation in this geometry, the
            start angle will be adjusted. Consider this parameter more like
            an offset to the start angle of stage rotation.

        Returns
        -------
        cera_parameters : dict
            The dictionary contains the following keys:

            + `"R"`: CERA's source-object distance (SOD)

            + `"D"`: CERA's source-detector distance (SDD)

            + `"ODD"`: object-detector distance (ODD = SDD - SOD)

            + `"a"`: CERA's a tilt

            + `"b"`: CERA's b tilt

            + `"c"`: CERA's c tilt

            + `"u0"`: Detector u offset (px)

            + `"v0"`: Detector v offset (px)

            + `"start_angle"`

            + `"volume_midpoint"`: dict

                - `"x"`, `"y"` and `"z"`

            + `"voxelsize"`: dict

                - `"x"`, `"y"` and `"z"`
        """

        cera_detector = self.detector.get_copy()

        # Number and size of pixels:
        nu  = cera_detector.pixels_u
        nv  = cera_detector.pixels_v
        psu = cera_detector.pitch_u
        psv = cera_detector.pitch_v

        # Default number of voxels for the reconstruction volume
        # is based on the detector size:
        n_voxels_x = nu
        n_voxels_y = nu
        n_voxels_z = nv

        # CERA's detector CS has its origin in the lower left corner instead of the center.
        # Let's move there:
        half_width  = psu*nu / 2.0
        half_height = psv*nv / 2.0

        cera_detector.center -= cera_detector.u.scaled(half_width) # add half a pixel in u direction??
        cera_detector.center += cera_detector.v.scaled(half_height) # subtract half a pixel in v direction??

        # The v axis points up instead of down:
        cera_detector.rotate_around_u(angle=math.pi)

        # Construct the CERA world coordinate system:
        # --------------------------------------------------
        # z axis points in v direction of our detector CS:
        cera_z = cera_detector.v.get_copy()
        cera_z.make_unit_vector()

        z0 = cera_z.x()
        z1 = cera_z.y()
        z2 = cera_z.z()

        O0 = self.stage.center.x()
        O1 = self.stage.center.y()
        O2 = self.stage.center.z()

        S0 = self.source.center.x()
        S1 = self.source.center.y()
        S2 = self.source.center.z()

        w0 = self.stage.w.x()
        w1 = self.stage.w.y()
        w2 = self.stage.w.z()

        # x axis points from source to stage (inverted), and perpendicular to cera_z (det v):
        t = -(z0*(O0-S0) + z1*(O1-S1) + z2*(O2-S2))/(z0*w0 + z1*w1 + z2*w2)
        d = self.source.center.distance(self.stage.center)
        SOD = math.sqrt(d*d - t*t)

        if SOD > 0:
            x0 = -(O0 - S0 + t*w0)/SOD
            x1 = -(O1 - S1 + t*w1)/SOD
            x2 = -(O2 - S2 + t*w2)/SOD
        else:
            # SOD == 0
            x0 = -1
            x1 = 0
            x2 = 0

        cera_x = Vector(x0, x1, x2)
        cera_x.make_unit_vector()

        cs_CERA = CoordinateSystem()
        cs_CERA.center = self.source.center.get_copy()
        cs_CERA.set_u_w(cera_x, cera_z)

        stage_in_CERA    = self.stage.get_copy()
        detector_in_CERA = cera_detector.get_copy()
        source_in_CERA   = self.source.get_copy()

        stage_in_CERA.change_reference_frame(ctsimu_world, cs_CERA)
        detector_in_CERA.change_reference_frame(ctsimu_world, cs_CERA)
        source_in_CERA.change_reference_frame(ctsimu_world, cs_CERA)

        # Source:
        xS = source_in_CERA.center.x()
        yS = source_in_CERA.center.y()
        zS = source_in_CERA.center.z()

        # Stage:
        xO = stage_in_CERA.center.x()
        yO = stage_in_CERA.center.y()
        zO = stage_in_CERA.center.z()
        uO = stage_in_CERA.u.unit_vector()
        vO = stage_in_CERA.v.unit_vector()
        wO = stage_in_CERA.w.unit_vector()

        # Detector:
        xD = detector_in_CERA.center.x()
        yD = detector_in_CERA.center.y()
        zD = detector_in_CERA.center.z()
        uD = detector_in_CERA.u.unit_vector()
        vD = detector_in_CERA.v.unit_vector()
        wD = detector_in_CERA.w.unit_vector()

        # Detector normal:
        nx = wD.x()
        ny = wD.y()
        nz = wD.z()

        # Intersection of CERA's x axis with the stage rotation axis = ceraVolumeMidpoint (new center of stage)
        xaxis = Vector(SOD, 0, 0)
        cera_volume_midpoint = source_in_CERA.center.get_copy()
        cera_volume_midpoint.subtract(xaxis)

        if xaxis.length() != 0:
            xaxis.make_unit_vector()

        world_volume_midpoint = change_reference_frame_of_point(cera_volume_midpoint, cs_CERA, ctsimu_world)

        cera_volume_relative_midpoint = cera_volume_midpoint.to(stage_in_CERA.center)
        midpoint_x = cera_volume_relative_midpoint.x()
        midpoint_y = cera_volume_relative_midpoint.y()
        midpoint_z = cera_volume_relative_midpoint.z()

        c = uD.x()   # x component of detector u vector is c-tilt
        a = wO.x()   # x component of stage w vector is a-tilt
        b = wO.y()   # y component of stage w vector is b-tilt

        # Intersection of x axis with detector (in px):
        efoc_x = xaxis.x() # 1
        efoc_y = xaxis.y() # 0
        efoc_z = xaxis.z() # 0

        E = nx*xD + ny*yD + nz*zD
        dv = nx*efoc_x + ny*efoc_y + nz*efoc_z

        if dv > 0:
            SDD_cera = abs((E - xS*nx - yS*ny - zS*nz)/dv)
        else:
            SDD_cera = 1

        SOD_cera = source_in_CERA.center.distance(cera_volume_midpoint)

        if SDD_cera != 0:
            voxelsize_u = psu * SOD_cera / SDD_cera
            voxelsize_v = psv * SOD_cera / SDD_cera
        else:
            voxelsize_u = 1
            voxelsize_v = 1

        # Intersection point of principal ray with detector:
        detector_intersection_point = xaxis.get_copy()
        detector_intersection_point.scale(-SDD_cera)

        stage_on_detector = detector_in_CERA.center.to(detector_intersection_point)
        ufoc = stage_on_detector.dot(uD)
        vfoc = stage_on_detector.dot(vD)
        wfoc = stage_on_detector.dot(wD)

        if psu > 0:
            ufoc_px = ufoc / psu
        else:
            ufoc_px = 0

        if psv > 0:
            vfoc_px = vfoc / psv
        else:
            vfoc_px = 0

        offset_u = ufoc_px - 0.5
        offset_v = vfoc_px - 0.5

        # Detector rotation relative to stage:
        cera_x = Vector(1, 0, 0)
        cera_y = Vector(0, 1, 0)
        cera_x.scale(vO.dot(cera_x))
        cera_y.scale(vO.dot(cera_y))

        v_in_xy_plane = cera_x.get_copy()
        v_in_xy_plane.add(cera_y)

        rot = v_in_xy_plane.angle(cera_y)

        # Add this start angle to the user-defined start angle:
        start_angle += (180.0 - math.degrees(rot))

        cera_parameters = {
            "R": SOD_cera,
            "D": SDD_cera,
            "ODD": SDD_cera - SOD_cera,
            "a": a,
            "b": b,
            "c": c,
            "u0": offset_u,
            "v0": offset_v,
            "start_angle": start_angle,
            "volume_midpoint": {
                "x": midpoint_x,
                "y": midpoint_y,
                "z": midpoint_z
            },
            "voxels": {
                "x": n_voxels_x,
                "y": n_voxels_y,
                "z": n_voxels_z
            },
            "voxelsize": {
                "x": voxelsize_u,
                "y": voxelsize_u,
                "z": voxelsize_v
            }
        }

        return cera_parameters

    def projection_matrix(self,
                         volumeCS:CoordinateSystem=None,
                         imageCS:CoordinateSystem=None,
                         mode:str=None):
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
            Pre-defined modes. Either `"OpenCT"` or `"CERA"` are supported.
            They override the `volumeCS` and `imageCS`, which can be set
            to `None` when using one of the pre-defined modes.

        Returns
        -------
        P : Matrix
            Projection matrix.

        Notes
        -----
        The image coordinate system (`imageCS`) should match the location,
        scale and orientation used by the reconstruction software, and is
        expressed in terms of the detector coordinate system.

        The detector coordinate system has its origin at the detector `center`,
        the `u` unit vector points in the row vector direction, and the
        `v` unit vector points in column vector direction (they are always
        assumed to be unit vectors).

        The `center` (origin) of the `imageCS` should be where the reconstruction
        software places the origin of its own projection image coordinate
        system. For example, CERA places it at the center of the lower-left
        pixel of the projection image.

        Similarly, a volume coordinate system (`volumeCS`) can be provided
        that describes the location, scale and orientation of the reconstruction
        volume with respect to the stage coordinate system.

        If the reconstruction software expects a different unit for the image
        or volume coordinate system (e.g. mm or voxels) than the world
        coordinates (e.g. mm), you can scale the basis vectors accordingly.
        For example, if you need a pixel and voxel coordinate system instead
        of a millimeter coordinate system, scale the basis vectors by the
        respective pixel and voxel size:

        ```python
        imageCS.u.scale(pixelSize_u)
        imageCS.v.scale(pixelSize_v)
        imageCS.w.scale(1.0) # Do not scale the detector normal!

        volumeCS.u.scale(voxelSize_u)
        volumeCS.v.scale(voxelSize_v)
        volumeCS.w.scale(voxelSize_w)
        ```
        """

        validModes = ["openct", "cera"]

        if mode is not None:
            if mode.lower() in validModes:  # Override imageCS
                image = CoordinateSystem()
                volume = CoordinateSystem()

                if mode.lower() == "openct":
                    """OpenCT places the origin of the image CS at the detector
                    center. The constructor places it at (0,0,0) automatically,
                    so there is nothing to do. Comments for illustration."""
                    # image.center.set_x(0)
                    # image.center.set_y(0)
                    # image.center.set_z(0)

                    """OpenCT's image CS is in mm units. We assume that all
                    other coordinate systems are in mm as well here (at least
                    when imported from JSON file). No scaling of the basis vectors is necessary."""
                    # image.u.scale(1.0)
                    # image.v.scale(1.0)
                    # image.w.scale(1.0)

                    volume.w.invert() # mirror reconstruction volume

                elif mode.lower() == "cera":
                    if self.detector.size_is_set():
                        """CERA places the origin of the image CS in the center
                        of the lower left pixel of the projection image."""
                        image.center.set_x(-self.detector.phys_width  / 2.0 + 0.5*self.detector.pitch_u)
                        image.center.set_y( self.detector.phys_height / 2.0 - 0.5*self.detector.pitch_v)
                        # image.center.set_z(0)

                        """CERA's unit of the image CS is in px, so we need to
                        scale the image CS basis vectors by the pixel size.
                        Also, v points up instead of down."""
                        image.u.scale( self.detector.pitch_u)
                        image.v.scale(-self.detector.pitch_v)

                        volume.w.invert() # mirror reconstruction volume
                    else:
                        raise RuntimeError("Detector size not set. To calculate a projection matrix for CERA, you need to set the size of the detector. Use the set_size() function of your detector object.")
            else:
                raise RuntimeError("Unsupported mode for projection matrix: \"{}\"".format(mode))
        else:
            if imageCS is not None:
                image = copy.deepcopy(imageCS)
            else:
                # Set a standard coordinate system. Results in pure
                # detector coordinate system after transformation.
                image = CoordinateSystem()

            if volumeCS is not None:
                volume = copy.deepcopy(volumeCS)
            else:
                # Set a standard coordinate system. Results in pure
                # stage coordinate system after transformation.
                volume = CoordinateSystem()

        source = copy.deepcopy(self.source)

        # Detach the image CS from the detector CS and
        # express it in terms of the world CS:
        image.change_reference_frame(self.detector, ctsimu_world)

        # Detach the volume CS from the stage CS and
        # express it in terms of the world CS:
        volume.change_reference_frame(self.stage, ctsimu_world)

        """The volume scale factors are derived from the lengths of the basis
        vectors of the volume CS ."""
        scale_volume_u = volume.u.length()
        scale_volume_v = volume.v.length()
        scale_volume_w = volume.w.length()

        """The image scale factors are derived from the lengths of the basis
        vectors of the image CS."""
        scale_image_u = image.u.length()
        scale_image_v = image.v.length()
        scale_image_w = image.w.length()

        # Save a source CS as seen from the detector CS. This is convenient to
        # later get the SDD, ufoc and vfoc:
        source_from_image = copy.deepcopy(self.source)
        source_from_image.change_reference_frame(ctsimu_world, image)

        # Make the volume CS the new world CS:
        source.change_reference_frame(ctsimu_world, volume)
        image.change_reference_frame(ctsimu_world, volume)
        volume.change_reference_frame(ctsimu_world, volume)

        # Translation vector from volume to source:
        xfoc = source.center.x()
        yfoc = source.center.y()
        zfoc = source.center.z()

        # Focus point on detector: principal, perpendicular ray.
        # In the detector coordinate system, ufoc and vfoc are the u and v coordinates
        # of the source center; SDD (perpendicular to detector plane) is source w coordinate.
        ufoc = source_from_image.center.x() / scale_image_u
        vfoc = source_from_image.center.y() / scale_image_v
        SDD  = abs(source_from_image.center.z())

        # Scale: volume units -> world units
        A = Matrix(values=[
            [scale_volume_u, 0, 0, 0],
            [0, scale_volume_v, 0, 0],
            [0, 0, scale_volume_w, 0],
            [0, 0, 0, 1]
        ])

        # Move origin to source (the origin of the camera CS)
        F = Matrix(values=[
            [1, 0, 0, xfoc],
            [0, 1, 0, yfoc],
            [0, 0, 1, zfoc]
        ])

        # Rotations:
        R = basis_transform_matrix(volume, image)

        # Projection onto detector and scaling (world units -> volume units):
        S = Matrix(values=[
            [-SDD/scale_image_u, 0, 0],
            [0, -SDD/scale_image_v, 0],
            [0, 0, -1.0/scale_image_w]
        ])

        # Shift in detector CS: (ufoc and vfoc must be in scaled units)
        T = Matrix(values=[
            [1, 0, ufoc],
            [0, 1, vfoc],
            [0, 0, 1]
        ])

        # Multiply all together:
        P = T * (S * (R * (F * A)))

        # Renormalize:
        lower_right = P.get(col=3, row=2)
        if lower_right != 0:
            P.scale(1.0/lower_right)
            P.set(col=3, row=2, value=1.0) # avoids rounding issues

        return P

    def create_detector_flat_field_rays(self):
        """ Calculate an analytical free beam intensity distribution
            picture for the given detector, to be used for an
            ideal flat field correction. """
        width      = self.detector.cols()
        height     = self.detector.rows()
        pixelSizeU = self.detector.pitch_u
        pixelSizeV = self.detector.pitch_v

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
        dx = self.detector.center.x()
        dy = self.detector.center.y()
        dz = self.detector.center.z()

        sx = self.source.center.x()
        sy = self.source.center.y()
        sz = self.source.center.z()

        # Vectors of the detector coordinate system:
        ux = self.detector.u.x()
        uy = self.detector.u.y()
        uz = self.detector.u.z()

        vx = self.detector.v.x()
        vy = self.detector.v.y()
        vz = self.detector.v.z()

        wx = self.detector.w.x()
        wy = self.detector.w.y()
        wz = self.detector.w.z()


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
   Detector Vector W: {}, {}, {}".format(alpha, dist, SDD, width, height, pixelSizeU, pixelSizeV, sx, sy, sz, dx, dy, dz, connectionLine.x(), connectionLine.y(), connectionLine.z(), ux, uy, uz, vx, vy, vz, wx, wy, wz))

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
                        pixel = self.detector.pixel_vector(x+(gx+1)*stepSize, y+(gy+1)*stepSize)

                        # Grid with no margin:
                        #if gridSize > 1:
                        #    stepSize = 1.0 / (gridSize-1)
                        #    pixel = self.detector.pixel_vector(x+gx*stepSize, y+gy*stepSize)
                        #else:
                        #    pixel = self.detector.pixel_vector_center(x, y)

                        distToSource = self.source.center.distance(pixel)

                        # Angle of incident rays:
                        vecSourceToPixel = Vector(pixel.x()-sx, pixel.y()-sy, pixel.z()-sz)
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

    def pixel_area_on_unit_sphere(self, A, B, C, D):
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

    def triangle_area_on_unit_sphere(self, A, B, C):
        # Source must be at (0, 0, 0) relative to given detector,
        # and A, B, C must be vectors pointing to triangle corners in
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

    def polygon_area_on_unit_sphere(self, polygon):
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

                area += self.triangle_area_on_unit_sphere(p1, p2, p3)

            return area
        else:
            return 0

    def create_detector_flat_field_sphere(self, *coverPolygons):
        """ Calculate an analytical free beam intensity distribution
            picture for the given detector, to be used for an
            ideal flat field correction.

            Geometrical approach using spherical geometry. """

        # Change to the detector coordinate system:
        D = copy.deepcopy(self.detector)
        S = copy.deepcopy(self.source)
        world = CoordinateSystem()  # will be initialized as world

        S.change_reference_frame(world, D)
        D.change_reference_frame(world, D)
        D.compute_geometry_parameters()

        # Source - Detector distance (SDD) defined by shortest distance between source and detector,
        # or distance between source and spot of highest intensity on detector.
        SDD = abs(S.center.z())

        # Calculate the area of the theoretical "brightest pixel" on the unit sphere:
        pu = D.pitch_u
        pv = D.pitch_v
        nRows = D.rows()
        nCols = D.cols()

        hpu = 0.5*pu
        hpv = 0.5*pv
        pA = Vector(SDD,  hpu,  hpv)
        pB = Vector(SDD, -hpu,  hpv)
        pC = Vector(SDD, -hpu, -hpv)
        pD = Vector(SDD,  hpu, -hpv)
        areaOfBrightestPixel = self.pixel_area_on_unit_sphere(pA, pB, pC, pD)

        # Full flat field image (without any clipping bodies):
        flatField = Image()
        flatField.shape(D.cols(), D.rows(), 0, flatField.getInternalDataType())

        # A second image with a clipping body under consideration: (both will be returned)
        clippedFlatField = None
        if len(coverPolygons) > 0:
            clippedFlatField = Image()
            clippedFlatField.shape(D.cols(), D.rows(), 0, flatField.getInternalDataType())

        # Upper left detector corner in world coordinates (remember: world is now the detector CS)
        p00 = D.pixel_vector(0, 0)

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
                pxSphericalArea  = self.polygon_area_on_unit_sphere(pixelPolygon)

                flatField.setPixel(x, y, pxSphericalArea)

                if len(coverPolygons) > 0:
                    for coverPolygon in coverPolygons:
                        pixelPolygon = pixelPolygon.clip(coverPolygon)

                    # Remove the intensity covered by the clipping polygon:
                    pixelPolygon.make_3D(z_component=SDD)
                    subarea = self.polygon_area_on_unit_sphere(pixelPolygon)
                    pxSphericalArea -= subarea

                    clippedFlatField.setPixel(x, y, pxSphericalArea)

            progress = 100*(float(x+1)/float(D.cols()))
            print("\rCalculating analytical intensity profile... {:0.1f}%    ".format(progress), end='')

        # Method #1: renormalize to area of theoretically brightest pixel:
        flatField.divide(areaOfBrightestPixel)
        if clippedFlatField is not None:
            clippedFlatField.divide(areaOfBrightestPixel)

        # Method #2: rescale maximum of actual brightest pixel to 1.0:
        #flatField.renormalize(newMin=0, newMax=1.0, currentMin=0)

        #flatField.save("ff.tif", numpy.dtype('float32'))

        print("\rCalculating analytical intensity profile... 100%  ")

        return flatField, clippedFlatField

    def solid_angle(self, l, m):
        """ Solid angle helper function for intensity profile. Approach by Florian Wohlgemuth. """
        if l != 0:
            return (l/abs(l)) * math.atan(abs(l)*m/math.sqrt(1.0+l**2+m**2))
        else:
            return 0

    def create_detector_flat_field_analytical(self):
        """ Calculate an analytical free beam intensity distribution
            picture for the given detector, to be used for an
            ideal flat field correction.

            Analytical approach by Florian Wohlgemuth. """

        width  = self.detector.cols()
        height = self.detector.rows()
        pitch_u = self.detector.pitch_u
        pitch_v = self.detector.pitch_v

        if(width is None):
            raise Exception("The detector width (in pixels) must be provided through a valid CTSimU JSON file.")
        if(height is None):
            raise Exception("The detector height (in pixels) must be provided through a valid CTSimU JSON file.")
        if(pitch_u is None):
            raise Exception("The pixel size (in mm) in u direction must be provided through a valid CTSimU JSON file.")
        if(pitch_v is None):
            raise Exception("The pixel size (in mm) in v direction must be provided through a valid CTSimU JSON file.")

        flatField = Image()
        flatField.shape(width, height, 0, flatField.getInternalDataType())

        # Positions of detector and source center:
        dx = self.detector.center.x()
        dy = self.detector.center.y()
        dz = self.detector.center.z()

        sx = self.source.center.x()
        sy = self.source.center.y()
        sz = self.source.center.z()

        # Vectors of the detector coordinate system:
        ux = self.detector.u.x()
        uy = self.detector.u.y()
        uz = self.detector.u.z()

        vx = self.detector.v.x()
        vy = self.detector.v.y()
        vz = self.detector.v.z()

        wx = self.detector.w.x()
        wy = self.detector.w.y()
        wz = self.detector.w.z()


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
        translation_vector = Vector(-sx, -sy, -sz)
        det.translate(translation_vector)
        det.compute_geometry_parameters()

        upperLeft_u = det.pixel_vector(0, 0).dot(self.detector.u)
        upperLeft_v = det.pixel_vector(0, 0).dot(self.detector.v)
        upperLeft_w = det.pixel_vector(0, 0).dot(self.detector.w)

        if upperLeft_w != 0:   # check if detector is not facing its edge towards the source
            for x in range(width):
                for y in range(height):
                    nu = x
                    nv = y
                    lambda0 = (upperLeft_u + nu*pitch_u) / upperLeft_w
                    lambda1 = (upperLeft_u + (nu+1)*pitch_u) / upperLeft_w
                    mu0     = (upperLeft_v + nv*pitch_v) / upperLeft_w
                    mu1     = (upperLeft_v + (nv+1)*pitch_v) / upperLeft_w

                    omega = self.solid_angle(lambda0, mu0) + self.solid_angle(lambda1, mu1) - self.solid_angle(lambda0, mu1) - self.solid_angle(lambda1, mu0)

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
        hpu = 0.5*det.pitch_u
        hpv = 0.5*det.pitch_v
        A = Vector(SDD,  hpu,  hpv)
        B = Vector(SDD, -hpu,  hpv)
        C = Vector(SDD, -hpu, -hpv)
        D = Vector(SDD,  hpu, -hpv)
        areaOfBrightestPixel = self.pixel_area_on_unit_sphere(A, B, C, D)
        flatField.divide(areaOfBrightestPixel)

        # Method #2: rescale actual maximum to 1.
        #flatField.renormalize(newMin=0, newMax=1.0, currentMin=0)

        #flatField.save("ff.tif", dataType="float32")
        return flatField


    def create_detector_flat_field(self):
        return create_detector_flat_field_analytical()

def _cera_bool(truth:bool) -> str:
    """Convert a Pythonian boolean into a CERA boolean string.

    Parameters
    ----------
    truth : bool
        Python boolean.

    Returns
    -------
    _cera_boolean : str
        Either `"true"` or `"false"`.
    """
    if isinstance(truth, str):
        # Catch if truth is already a string
        if truth == "true" or truth == "false":
            return truth
        else:
            raise Exception(f"Cannot convert '{truth}' into the CERA boolean 'true' or 'false'.")

    if truth:
        return "true"

    return "false"

def create_CERA_config(geo:'Geometry', projection_file_pattern:str, basename:str, save_dir:str=None, n_projections:int=None, flip_u:bool=False, flip_v:bool=True, projection_datatype:str="float32", projection_filetype:str="tiff", projection_byteorder:str="little", projection_headersize:int=0, start_angle:float=0, total_angle:float=360, scan_direction="CCW", voxels_x:int=None, voxels_y:int=None, voxels_z:int=None, voxelsize_x:float=None, voxelsize_y:float=None, voxelsize_z:float=None, i0max:float=60000, output_datatype:str="float32", matrices:list=None):
    """Write a CERA config file for the given geometry.

    A circular trajectory for the given angular range is assumed, all parameters
    of the output config file will reflect a circular behaviour (obeying static tilts).
    For non-circular trajectories, provide a list of projection matrices.

    Parameters
    ----------
    geo : ctsimu.geometry.Geometry
        A geometry that should represent the CT setup at frame zero, as seen
        by the reconstruction software.

    projection_file_pattern : str
        Pattern for the sequential projection files.

        Example: `"../projections/corrected/example_%04d.tif"`

    basename : str
        Base name for the configuration files and table of projection matrices.

    save_dir : str
        Directory where the configuration files will be stored.

        Standard value: `None` (local script directory)

    n_projections : int
        Number of projections. Set to `None` if number of projections
        should be inferred from the number of provided projection matrices.

        Standard value: `None`

    flip_u : bool
        Flip projection images horizontally before reconstruction?

        Standard value: `False`

    flip_v : bool
        Flip projection images vertically before reconstruction?

        Standard value: `True`

    projection_datatype : str
        Data type of the projection images, as well as possible bright and dark
        images and the bad pixel map.

        Allowed values: `"uint16"`, `"float32"`

        Standard value: `"float32"`

    projection_filetype : str
        File type of the projection images, as well as possible bright and dark
        images and the bad pixel map.

        Allowed values: `"raw"` or `"tiff"`

        Standard value: `"tiff"`

    projection_headersize : int
        For RAW projection images: header size to skip (in bytes).

        Standard value: `0`

    projection_byteorder : str
        For RAW projection images: endianness of the files.

        Allowed values: `"little"` or `"big"`

        Standard value: `"little"`

    start_angle : float
        Reconstruction start angle (in degrees). The start angle can be
        tuned to change the in-plane rotation of the reconstruction images.
        Depending on the current stage rotation in this geometry, the
        start angle will be adjusted. Consider this parameter more like
        an offset to the start angle of stage rotation.

        Standard value: `0`

    total_angle : float
        Total angular range of the CT scan (in degrees).

        Standard value: `360`

    scan_direction : str
        Direction of stage rotation, either `"CCW"` for counter-clockwise
        or `"CW"` for clockwise rotation.

        Standard value: `"CCW"`

    voxels_x : int
        Number of voxels in x direction of the reconstruction volume.
        Set to `None` for a default value based on the detector pixels.

        Standard value: `None`

    voxels_y : int
        Number of voxels in y direction of the reconstruction volume.
        Set to `None` for a default value based on the detector pixels.

        Standard value: `None`

    voxels_z : int
        Number of voxels in z direction of the reconstruction volume.
        Set to `None` for a default value based on the detector pixels.

        For helix scans, this parameter should be increased because the
        default voxel size in z direction corresponds to the detector height,
        whereas a helix scan usually covers much more than the detector height.

        Standard value: `None`

    voxelsize_x : float
        Voxel size in x direction of the reconstruction volume.
        Set to `None` for a default value based on the detector pixel size
        and magnification.

        Standard value: `None`

    voxelsize_y : float
        Voxel size in y direction of the reconstruction volume.
        Set to `None` for a default value based on the detector pixel size
        and magnification.

        Standard value: `None`

    voxelsize_z : float
        Voxel size in z direction of the reconstruction volume.
        Set to `None` for a default value based on the detector pixel size
        and magnification.

        Standard value: `None`

    i0max : float
        Grey value for the maximum free-beam intensity in the
        projection images.

        Standard value: `60000`

    output_datatype : str
        Data type for the reconstruction volume output file.
        Either `"float32"` or `"uint16"`.

        Standard value: `"float32"`

    matrices : list
        List of projection matrices of type `ctsimu.primitives.Matrix`.
        One matrix for each projection image is required. If this parameter
        is set to `None`, the configuration will be set up for a circular scan
        trajectory (obeying static tilts).
    """

    now = datetime.now()

    cera_parameters = geo.get_CERA_standard_circular_parameters(start_angle=start_angle)

    cera_config_filename = join_dir_and_filename(save_dir, f"{basename}.config")
    touch_directory(cera_config_filename)

    big_endian = False
    if projection_byteorder == "big":
        big_endian = True

    cera_projection_filetype = "tiff"
    if projection_filetype == "raw":
        if projection_datatype == "uint16":
            cera_projection_filetype = "raw_uint16"
        elif projection_datatype == "float32":
            cera_projection_filetype = "raw_float"
        else:
            raise Exception(f"Projection datatype not supported by CERA: {projection_datatype}. Supported datatypes: 'uint16' and 'float32'.")

    # Use default number of voxels if none is given:
    if voxels_x is None:
        voxels_x = cera_parameters["voxels"]["x"]
    if voxels_y is None:
        voxels_y = cera_parameters["voxels"]["y"]
    if voxels_z is None:
        voxels_z = cera_parameters["voxels"]["z"]

    # Use default voxel size if none is given:
    if voxelsize_x is None:
        voxelsize_x = cera_parameters["voxelsize"]["x"]
    if voxelsize_y is None:
        voxelsize_y = cera_parameters["voxelsize"]["y"]
    if voxelsize_z is None:
        voxelsize_z = cera_parameters["voxelsize"]["z"]

    # Flip scan direction:
    # we assume object scan direction, CERA assumes gantry scan direction.
    if scan_direction == "CW":
        cera_scan_direction = "CCW"
    else:
        cera_scan_direction = "CW"

    midpoint_comment = ""

    # Projection matrix file
    # -------------------------------
    if isinstance(matrices, list):
        nMatrices = len(matrices)
        if n_projections is None:
            n_projections = nMatrices

        if nMatrices != n_projections:
            raise Exception("Error in write_CERA_config: given number of projections does not match number of projection matrices.")

        # Volume midpoints will be set to zero when using projection matrices,
        # their values for circular trajectory reconstructions is commented
        # so they can be easily activated again if the user decides not to
        # use projection matrices.
        midpoint_comment = "0 # "

        projTableString = """projtable.txt version 3
{timestring}

# format: angle / entries of projection matrices
{nMatrices}
""".format(
        nMatrices=nMatrices,
        timestring=now.strftime("%a %b %d %H:%M:%S %Y")
        )

        for i in range(nMatrices):
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

        projtable_filename = f"{basename}_projtable.txt"
        cera_projtable_filename = join_dir_and_filename(save_dir, projtable_filename)
        touch_directory(cera_projtable_filename)

        with open(cera_projtable_filename, 'w') as f:
            f.write(projTableString)
            f.close()

    if n_projections is None:
        raise Exception("write_CERA_config: Please provide the number of projection images in the parameter `n_projections`.")

    # CERA config file
    # ---------------------
    configFileString = """#CERACONFIG

[Projections]
NumChannelsPerRow = {nCols}
NumRows = {nRows}
PixelSizeU = {psu}
PixelSizeV = {psv}
Rotation = None
FlipU = {flip_u}
FlipV = {flip_v}
Padding = 0
BigEndian = {big_endian}
CropBorderRight = 0
CropBorderLeft = 0
CropBorderTop = 0
CropBorderBottom = 0
BinningFactor = None
SkipProjectionInterval = 1
ProjectionDataDomain = Intensity
RawHeaderSize = {projection_headersize}

[Volume]
SizeX = {volx}
SizeY = {voly}
SizeZ = {volz}
# Midpoints are only necessary for reconstructions
# without projection matrices.
MidpointX = {midpoint_comment}{midpointx}
MidpointY = {midpoint_comment}{midpointy}
MidpointZ = {midpoint_comment}{midpointz}
VoxelSizeX = {vsx}
VoxelSizeY = {vsy}
VoxelSizeZ = {vsz}
OutputDatatype = {output_datatype}

[CustomKeys]
NumProjections = {nProjections}
ProjectionFileType = {projection_filetype}
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
AquisitionDirection = {scanDirection}
a = {a}
b = {b}
c = {c}
ProjectionMatrixFilename = {projtable_filename}

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
    psu=geo.detector.pitch_u,
    psv=geo.detector.pitch_v,
    flip_u=_cera_bool(flip_u),
    flip_v=_cera_bool(flip_v),
    big_endian=_cera_bool(big_endian),
    projection_headersize=projection_headersize,
    volx=int(voxels_x),
    voly=int(voxels_y),
    volz=int(voxels_z),
    midpoint_comment=midpoint_comment,
    midpointx=cera_parameters["volume_midpoint"]["x"],
    midpointy=cera_parameters["volume_midpoint"]["y"],
    midpointz=cera_parameters["volume_midpoint"]["z"],
    vsx=voxelsize_x,
    vsy=voxelsize_y,
    vsz=voxelsize_z,
    output_datatype=output_datatype,
    nProjections=int(n_projections),
    projection_filetype=cera_projection_filetype,
    projFilePattern=projection_file_pattern,
    SOD=cera_parameters["R"],
    SDD=cera_parameters["D"],
    offu=cera_parameters["u0"],
    offv=cera_parameters["v0"],
    startAngle=cera_parameters["start_angle"],
    scanAngle=total_angle,
    scanDirection=cera_scan_direction,
    a=cera_parameters["a"],
    b=cera_parameters["b"],
    c=cera_parameters["c"],
    projtable_filename=projtable_filename,
    i0max=i0max
    )

    with open(cera_config_filename, 'w') as f:
        f.write(configFileString)
        f.close()

def create_OpenCT_config(geo:'Geometry', filename:str=None, variant:str="free", projection_files:list=None, projection_dir:str=None, flip_u:bool=False, flip_v:bool=False, projection_datatype:str="float32", projection_filetype:str="tiff", projection_headersize:int=0, projection_byteorder:str="little", detector_coordinate_frame="OriginAtDetectorCenter.VerticalAxisRunningDownwards", detector_coordinate_dimension="Length", total_angle:float=None, scan_direction:str="CCW", matrices:list=None, volumename:str=None, bb_center_x:float=0, bb_center_y:float=0, bb_center_z:float=0, voxels_x:int=None, voxels_y:int=None, voxels_z:int=None, voxelsize_x:float=None, voxelsize_y:float=None, voxelsize_z:float=None, bright_image_dir:str=None, bright_images:list=None, dark_image:str=None, bad_pixel_mask:str=None) -> dict:
    """Create an OpenCT free trajectory CBCT configuration and optionally write to file.

    Parameters
    ----------
    geo : ctsimu.geometry.Geometry
        A geometry that should represent the CT setup at frame zero, as seen
        by the reconstruction software.

    variant : str
        Which variant of the OpenCT file format is created.
        Options: `"free"` and `"circular"`.

        Standard value: `"free"`

    filename : str
        Path and filename for the OpenCT configuration file to be written.
        If no file should be written, set this to `None`.

        Standard value: `None`

    projection_files : list
        List of projection file names.

    projection_dir : str
        Path where the projection images are stored.

        Standard value: `None`

    flip_u : bool
        Flip projection images horizontally before reconstruction?

        Standard value: `False`

    flip_v : bool
        Flip projection images vertically before reconstruction?

        Standard value: `False`

    projection_datatype : str
        Data type of the projection images, as well as possible bright and dark
        images and the bad pixel map.

        Allowed values: `"uint8"`, `"uint16"`, `"uint32"`, `"int8"`, `"int16"`, `"int32"`, `"float32"`

        Standard value: `"float32"`

    projection_filetype : str
        File type of the projection images, as well as possible bright and dark
        images and the bad pixel map.

        Allowed values: `"raw"` or `"tiff"`

        Standard value: `"tiff"`

    projection_headersize : int
        For RAW projection images: header size to skip (in bytes).

        Standard value: `0`

    projection_byteorder : str
        For RAW projection images: endianness of the files.

        Allowed values: `"little"` or `"big"`

        Standard value: `"little"`

    detector_coordinate_frame : str
        A string that defines the orientation of the detector coordinate system
        with respect to the projection images.

        Possible values:

        + `"OriginAtDetectorCenter.VerticalAxisRunningUpwards"`
        + `"OriginAtDetectorCenter.VerticalAxisRunningDownwards"`
        + `"OriginAtDetectorTopLeftCorner"`
        + `"OriginAtDetectorBottomLeftCorner"`

        Standard value: `"OriginAtDetectorCenter.VerticalAxisRunningDownwards"`

    detector_coordinate_dimension : str
        A string that defines the unit of the detector coordinate system.

        Possible values:

        + `"Length"` for physical length units (usually mm).
        + `"PixelCount"` if projection matrices refer to a pixel coordinate system.

        Standard value: `"Length"`

    total_angle : float
        Total angular range (in deg) of the CT scan. This parameter is only
        really needed for the strict circular trajectory variant of the file
        format, but not for free trajectories.

        Standard value: `None`

    scan_direction : str
        Direction of stage rotation, either `"CCW"` for counter-clockwise
        or `"CW"` for clockwise rotation. Only relevant for circular
        trajectory variant.

        Standard value: `"CCW"`

    matrices : list
        List of projection matrices of type `ctsimu.primitives.Matrix`.
        One matrix for each projection image is required for the free trajectory
        variant of the OpenCT file format. For the circular trajectory variant,
        the matrices are not required.

        Standard value: `None`

    volumename : str
        Optional name for the reconstruction volume.

        Standard value: `None`

    bb_center_x : float
        Position (in x direction) of the reconstruction volume bounding box.
        Center position in respect to the stage coordinate system.

        Standard value: `0`

    bb_center_y : float
        Position (in y direction) of the reconstruction volume bounding box.
        Center position in respect to the stage coordinate system.

        Standard value: `0`

    bb_center_z : float
        Position (in z direction) of the reconstruction volume bounding box.
        Center position in respect to the stage coordinate system.

        Standard value: `0`

    voxels_x : int
        Number of voxels in x direction of the reconstruction volume.
        Set to `None` for a default value based on the detector pixels.

        Standard value: `None`

    voxels_y : int
        Number of voxels in y direction of the reconstruction volume.
        Set to `None` for a default value based on the detector pixels.

        Standard value: `None`

    voxels_z : int
        Number of voxels in z direction of the reconstruction volume.
        Set to `None` for a default value based on the detector pixels.

        For helix scans, this parameter should be increased because the
        default voxel size in z direction corresponds to the detector height,
        whereas a helix scan usually covers much more than the detector height.

        Standard value: `None`

    voxelsize_x : float
        Voxel size in x direction of the reconstruction volume.
        Set to `None` for a default value based on the detector pixel size
        and magnification.

        Standard value: `None`

    voxelsize_y : float
        Voxel size in y direction of the reconstruction volume.
        Set to `None` for a default value based on the detector pixel size
        and magnification.

        Standard value: `None`

    voxelsize_z : float
        Voxel size in z direction of the reconstruction volume.
        Set to `None` for a default value based on the detector pixel size
        and magnification.

        Standard value: `None`

    bright_image_dir : str
        Optional directory where bright correction images are stored.

        Standard value: `None`

    bright_images : list
        List of file names of bright correction images.

        Standard value: `None`

    dark_image : str
        Path to a dark correction image.

        Standard value: `None`

    bad_pixel_mask : str
        Path to a bad pixel mask for correction.

        Standard value: `None`

    Returns
    -------
    openct_config : dict
        Dictionary that represents the JSON structure of the OpenCT file.
    """

    # Convert some settings to OpenCT keywords:
    projection_datatype = convert(openct_converter["datatype"], projection_datatype)
    projection_byteorder = convert(openct_converter["endian"], projection_byteorder)

    n_projections = 0
    if projection_files is not None:
        n_projections = len(projection_files)
    else:
        # Initialize as an empty list:
        projection_files = list()

    n_matrices = 0
    openct_matrices = list()
    if matrices is not None:
        if isinstance(matrices, list):
            n_matrices = len(matrices)
            if n_matrices != n_projections:
                warnings.warn(f"In write_openCT_config: Number of projection matrices ({n_matrices}) does not match number of projection images ({n_projections}).")

            for m in matrices:
                openct_matrices.append(m.value.tolist())
        else:
            raise Exception("Error in write_openCT_config: 'matrices' must be a list.")

    mirror_detector_axis = ""
    if (flip_u is True) and (flip_v is False):
        mirror_detector_axis = "U"
    elif (flip_u is False) and (flip_v is True):
        mirror_detector_axis = "V"
    elif (flip_u is True) and (flip_v is True):
        mirror_detector_axis = "UV"

    # ### Geometry parameters
    # Take some standard parameters from CERA:
    cera_parameters = geo.get_CERA_standard_circular_parameters()

    # Use default number of voxels if none is given:
    if voxels_x is None:
        voxels_x = cera_parameters["voxels"]["x"]
    if voxels_y is None:
        voxels_y = cera_parameters["voxels"]["y"]
    if voxels_z is None:
        voxels_z = cera_parameters["voxels"]["z"]

    # Use default voxel size if none is given:
    if voxelsize_x is None:
        voxelsize_x = cera_parameters["voxelsize"]["x"]
    if voxelsize_y is None:
        voxelsize_y = cera_parameters["voxelsize"]["y"]
    if voxelsize_z is None:
        voxelsize_z = cera_parameters["voxelsize"]["z"]

    # Calculate bounding box size:
    bb_size_x = voxels_x * voxelsize_x
    bb_size_y = voxels_y * voxelsize_y
    bb_size_z = voxels_z * voxelsize_z

    SOD = cera_parameters["R"]
    SDD = cera_parameters["D"]
    ODD = cera_parameters["ODD"]

    openct_variant = "FreeTrajectoryCBCTScan"
    if variant == "circular":
        openct_variant = "CircularTrajectoryCBCTScan"

        # Remove matrices to avoid license issues,
        # though technically matrices are allowed in the circular
        # trajectory format.
        openct_matrices = None

        # For clockwise scan direction, the projection filenames
        # need to be in reverse order:
        if scan_direction == "CW":
            projection_files.reverse()

    openct_config = {
        "version": {"major": 1, "minor": 0},
        "OpenCTJSON": {
            "versionMajor": 1,
            "versionMinor": 0,
            "revisionNumber": 0,
            "variant": openct_variant
        },
        "hints": None,
        "units": {
            "length": "Millimeter",
            "angle": "Degree"
        },
        "volumeName":  volumename,
        "projections": {
            "numProjections":  n_projections,
            "intensityDomain": True,
            "images": {
                "dataType":  projection_datatype,
                "fileType":  projection_filetype.upper(),
                "skipBytes": projection_headersize,
                "endianness": projection_byteorder,
                "directory": projection_dir,
                "files": projection_files
            },
            "detectorCoordinateFrame": detector_coordinate_frame,
            "detectorCoordinateDimension": detector_coordinate_dimension,
            "matrices": openct_matrices
        },
        "geometry": {
            "detectorPixel": [
                int(geo.detector.cols()),
                int(geo.detector.rows())
            ],
            "detectorSize": [
                geo.detector.phys_width,
                geo.detector.phys_height
            ],
            "distanceSourceObject": SOD,
            "distanceObjectDetector": ODD,
            "mirrorDetectorAxis": mirror_detector_axis,
            "skipAngle": 0,
            "totalAngle": total_angle,
            "objectBoundingBox": {
                "centerXYZ": [bb_center_x, bb_center_y, bb_center_z],
                "sizeXYZ": [bb_size_x, bb_size_y, bb_size_z]
            }
        },
        "corrections": {
            "brightImages": {
                "dataType":  projection_datatype,
                "fileType":  projection_filetype.upper(),
                "skipBytes": projection_headersize,
                "endianness": projection_byteorder,
                "directory": bright_image_dir,
                "files":     bright_images
            },
            "darkImage": {
                "file":     dark_image,
                "dataType": projection_datatype,
                "fileType": projection_filetype.upper(),
                "skipBytes": projection_headersize,
                "endianness": projection_byteorder
            },
            "badPixelMask": {
                "file":     bad_pixel_mask,
                "dataType": projection_datatype,
                "fileType": projection_filetype.upper(),
                "skipBytes": projection_headersize,
                "endianness": projection_byteorder
            }
        }
    }

    # Clean up corrections if they are not needed:
    if (bright_images is None) and (dark_image is None) and (bad_pixel_mask is None):
        # Set all corrections to null if none are required.
        openct_config["corrections"] = None
    else:
        if bright_images is None:
            openct_config["corrections"]["brightImages"] = None
        if dark_image is None:
            openct_config["corrections"]["darkImage"] = None
        if bad_pixel_mask is None:
            openct_config["corrections"]["badPixelMask"] = None

    if filename is not None:
        write_json_file(filename=filename, dictionary=openct_config)

    return openct_config