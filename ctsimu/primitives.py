# -*- coding: UTF-8 -*-

import math
import numpy
import copy

from .helpers import *

class Matrix:
    def __init__(self, rows=None, cols=None, values=None):
        self.rows  = None
        self.cols  = None
        self.value = None

        if values is None:
            if not((rows is None) or (cols is None)):
                self.rows  = rows
                self.cols  = cols

                clist = cols*[0]
                self.value = []
                for i in range(rows):
                    self.value.append(clist.copy())
            else:
                raise Exception("Cannot initialize matrix.")
        else:
            self.value = values
            if len(values) > 0:
                self.rows = len(self.value)
                if len(self.value[0]) > 0:
                    self.cols = len(self.value[0])

    def __str__(self):
        string = ""
        for r in range(self.rows):
            string += "["
            for c in range(self.cols):
                string += "{space}{x:.7f}  ".format(x=self.value[r][c], space=" "*(self.value[r][c]>=0))

            string += "]\n"

        return string

    def __mul__(self, other):
        if isinstance(other, Vector):
            if self.cols >= 3 and self.rows >= 3:
                result = Vector()
                for r in range(3):
                    s = 0
                    for c in range(3):
                        s += self.value[r][c] * other.get(c)

                    result.setIdx(i=r, value=s)

                return result
            else:
                raise Exception("Error: Cannot multiply matrix of size ({rows} rows x {cols} columns) with a 3-vector.".format(rows=self.rows, cols=self.cols))
        elif isinstance(other, Matrix):
            result_rows = self.rows
            result_cols = other.cols
            result = Matrix(rows=result_rows, cols=result_cols)

            for row in range(result_rows):
                for col in range(result_cols):
                    s = 0
                    for idx in range(min(self.cols, other.rows)):
                        s += self.get(col=idx, row=row) * other.get(col=col, row=idx)

                    result.set(row=row, col=col, value=s)

            return result
        else:
            raise Exception("Matrix multiplication failed.")

    def set(self, col, row, value):
        self.value[row][col] = value

    def get(self, col, row):
        if len(self.value) > row:
            if len(self.value[row]) > col:
                return self.value[row][col]

        raise Exception("Matrix.get(col={}, row={}): requested index does not exist.".format(col, row))

    def scale(self, factor):
        for row in range(self.rows):
            for col in range(self.cols):
                self.set(col=col, row=row, value=self.get(col=col, row=row)*factor)

class Vector:
    """ A 3D vector in space. """

    def __init__(self, x=0, y=0, z=0):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        
        # The following member variables are private and should not
        # be accessed from outside. They are invalidated whenever the
        # vector changes.
        self._unitVector = None
        self._length = None

        self.update()

    def __str__(self):
        return "({spaceX}{x:.7f}, {spaceY}{y:.7f}, {spaceZ}{z:.7f})".format(x=self.x, y=self.y, z=self.z, spaceX=" "*(self.x>=0), spaceY=" "*(self.y>=0), spaceZ=" "*(self.z>=0))

    def __add__(self, other):
        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z
        return Vector(x, y, z)

    def __sub__(self, other):
        x = self.x - other.x
        y = self.y - other.y
        z = self.z - other.z
        return Vector(x, y, z)

    def __mul__(self, other):
        x = self.x * other.x
        y = self.y * other.y
        z = self.z * other.z
        return Vector(x, y, z)

    def __truediv__(self, other):
        x = self.x / other.x
        y = self.y / other.y
        z = self.z / other.z
        return Vector(x, y, z)

    def __floordiv__(self, other):
        x = self.x // other.x
        y = self.y // other.y
        z = self.z // other.z
        return Vector(x, y, z)

    def update(self):
        """Called when vector is changed.

        Invalidates private member variables, so they will be re-calculated the
        next time their public getter functions are called.
        """
        self._unitVector = None
        self._length = None

    def get(self, i):
        if i == 0:
            return self.x
        elif i == 1:
            return self.y
        elif i == 2:
            return self.z

        return 0

    def setx(self, value):
        self.x = float(value)
        self.update()

    def sety(self, value):
        self.y = float(value)
        self.update()

    def setz(self, value):
        self.z = float(value)
        self.update()

    def setxy(self, x=0, y=0):
        """ Set x and y component (relevant for 2D computations) """
        self.x = float(x)
        self.y = float(y)
        self.update()

    def set(self, x=0, y=0, z=0):
        """ Set all vector components. """
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.update()

    def setIdx(self, i, value):
        if i == 0:
            self.setx(value)
        elif i == 1:
            self.sety(value)
        elif i == 2:
            self.setz(value)

    def length(self):
        """Get the length of the vector."""
        if self._length is None:
            self._length = math.sqrt(self.x**2 + self.y**2 + self.z**2)

        return self._length

    def angle(self, other):
        """ Calculate angle between this vector and another vector. """
        dotProd = self.dot(other)
        l1 = self.length()
        l2 = other.length()
        lp = l1 * l2

        if lp > 0:
            cs = dotProd / lp
            alpha = 0

            # Avoid out-of-domain due to rounding errors:
            if cs >= 1.0:
                alpha = 0
            elif cs <= -1.0:
                alpha = math.pi
            else:
                alpha = math.acos(cs)

            return alpha
        else:
            return 0

    def makeUnitVector(self):
        """ Normalize vector length to 1. """
        vectorLength = self.length()
        if vectorLength != 0:
            if vectorLength != 1.0:
                self.x /= float(vectorLength)
                self.y /= float(vectorLength)
                self.z /= float(vectorLength)
                self.update()
        else:
            raise Exception("Unit vector: a zero length vector cannot be converted into a unit vector.")

    def unitVector(self):
        if self._unitVector is None:
            self._unitVector = Vector(self.x, self.y, self.z)
            self._unitVector.makeUnitVector()

        return self._unitVector

    def scale(self, factor):
        """ Scale vector length by a certain factor. """
        self.x *= factor
        self.y *= factor
        self.z *= factor
        self.update()

    def scaled(self, factor):
        """ Return a copy of this vector, scaled by a certain factor. """
        return Vector(self.x*factor, self.y*factor, self.z*factor)

    def add(self, vec):
        """ Add another vector to this vector. """
        self.x += vec.x
        self.y += vec.y
        self.z += vec.z
        self.update()

    def subtract(self, vec):
        """ Subtract another vector from this vector. """
        self.x -= vec.x
        self.y -= vec.y
        self.z -= vec.z
        self.update()

    def multiply(self, vec):
        """ Multiply another vector element-wise to this vector. """
        self.x *= vec.x
        self.y *= vec.y
        self.z *= vec.z
        self.update()

    def divide(self, vec):
        """ Element-wise divide this vector by another vector. """
        if vec.x != 0 and vec.y != 0 and vec.z != 0:
            self.x /= vec.x
            self.y /= vec.y
            self.z /= vec.z
            self.update()
        else:
            raise ZeroDivisionError("Vector division by zero: one component of divisor is zero: {}".format(vec))

    def square(self):
        self.x *= self.x
        self.y *= self.y
        self.z *= self.z
        self.update()

    def sqrt(self):
        self.x = math.sqrt(self.x)
        self.y = math.sqrt(self.y)
        self.z = math.sqrt(self.z)
        self.update()

    def distance(self, vec):
        """ Distance between target points of this and another vector. """
        return math.sqrt(math.pow(self.x-vec.x, 2) + math.pow(self.y-vec.y, 2) + math.pow(self.z-vec.z, 2))

    def dot(self, other):
        """ Calculate vector dot product. """
        return (self.x*other.x + self.y*other.y + self.z*other.z)

    def cross_z(self, other):
        """ Calculate the z component of the cross product. """
        return self.x*other.y - self.y*other.x

    def cross(self, other):
        """ Calculate vector cross product. """
        cx = self.y*other.z - self.z*other.y
        cy = self.z*other.x - self.x*other.z
        cz = self.x*other.y - self.y*other.x

        c = Vector(cx, cy, cz)
        return c

    def inverse(self):
        return Vector(-self.x, -self.y, -self.z)

    def invert(self):
        self.x = -self.x
        self.y = -self.y
        self.z = -self.z
        self.update()

    def rotate(self, axis, angle):
        """ Rotate vector around given axis by given angle [rad]. """

        # Implementing a general rotation matrix.
        cs = math.cos(angle)
        sn = math.sin(angle)

        vx = self.x
        vy = self.y
        vz = self.z

        nx = axis.unitVector().x
        ny = axis.unitVector().y
        nz = axis.unitVector().z

        rx = vx*(nx*nx*(1.0-cs) + cs)    + vy*(nx*ny*(1.0-cs) - nz*sn) + vz*(nx*nz*(1.0-cs) + ny*sn)
        ry = vx*(ny*nx*(1.0-cs) + nz*sn) + vy*(ny*ny*(1.0-cs) + cs)    + vz*(ny*nz*(1.0-cs) - nx*sn)
        rz = vx*(nz*nx*(1.0-cs) - ny*sn) + vy*(nz*ny*(1.0-cs) + nx*sn) + vz*(nz*nz*(1.0-cs) + cs)

        self.set(rx, ry, rz)

    @staticmethod
    def connection(vecFrom, vecTo):
        return Vector(vecTo.x-vecFrom.x, vecTo.y-vecFrom.y, vecTo.z-vecFrom.z)

class Vector2D:
    """A 2D vector in a plane.

    Separate implementation from 3D vectors to gain speed.
    """

    def __init__(self, x=0, y=0):
        self.x = float(x)
        self.y = float(y)

        # The following member variables are private and should not
        # be accessed from outside. They are invalidated whenever the
        # vector changes.
        self._unitVector = None
        self._length = None

        self.update()

    def __str__(self):
        return "({spaceX}{x:.7f}, {spaceY}{y:.7f})".format(x=self.x, y=self.y, spaceX=" "*(self.x>=0), spaceY=" "*(self.y>=0))

    def __add__(self, other):
        x = self.x + other.x
        y = self.y + other.y
        return Vector2D(x, y)

    def __sub__(self, other):
        x = self.x - other.x
        y = self.y - other.y
        return Vector2D(x, y)

    def __mul__(self, other):
        x = self.x * other.x
        y = self.y * other.y
        return Vector2D(x, y)

    def __truediv__(self, other):
        x = self.x / other.x
        y = self.y / other.y
        return Vector2D(x, y)

    def __floordiv__(self, other):
        x = self.x // other.x
        y = self.y // other.y
        return Vector2D(x, y)

    def update(self):
        """Called when vector is changed.

        Invalidates private member variables, so they will be re-calculated the
        next time their public getter functions are called.
        """
        self._unitVector = None
        self._length = None

    def setx(self, value):
        self.x = value
        self.update()

    def sety(self, value):
        self.y = value
        self.update()

    def set(self, x=0, y=0):
        """ Set all vector components. """
        self.x = x
        self.y = y
        self.update()

    def length(self):
        """Get the length of the vector."""
        if self._length is None:
            self._length = math.sqrt(self.x**2 + self.y**2)

        return self._length

    def angle(self, other):
        """ Calculate angle between this vector and another vector. """
        dotProd = self.dot(other)
        l1 = self.length()
        l2 = other.length()
        lp = l1 * l2

        if lp > 0:
            cs = dotProd / lp
            alpha = 0

            # Avoid out-of-domain due to rounding errors:
            if cs >= 1.0:
                alpha = 0
            elif cs <= -1.0:
                alpha = math.pi
            else:
                alpha = math.acos(cs)

            return alpha
        else:
            return 0

    def makeUnitVector(self):
        """ Normalize vector length to 1. """
        vectorLength = self.length()
        if vectorLength != 0:
            if vectorLength != 1.0:
                self.x /= float(vectorLength)
                self.y /= float(vectorLength)
                self.update()
        else:
            raise Exception("Unit vector: a zero length vector cannot be converted into a unit vector.")

    def unitVector(self):
        if self._unitVector is None:
            self._unitVector = Vector2D(self.x, self.y)
            self._unitVector.makeUnitVector()

        return self._unitVector

    def scale(self, factor):
        """ Scale vector length by a certain factor. """
        self.x *= factor
        self.y *= factor
        self.update()

    def add(self, vec):
        """ Add another vector to this vector. """
        self.x += vec.x
        self.y += vec.y
        self.update()

    def subtract(self, vec):
        """ Subtract another vector from this vector. """
        self.x -= vec.x
        self.y -= vec.y
        self.update()

    def multiply(self, vec):
        """ Multiply another vector element-wise to this vector. """
        self.x *= vec.x
        self.y *= vec.y
        self.update()

    def divide(self, vec):
        """ Element-wise divide this vector by another vector. """
        if vec.x != 0 and vec.y != 0:
            self.x /= vec.x
            self.y /= vec.y
            self.update()
        else:
            raise ZeroDivisionError("Vector division by zero: one component of divisor is zero: {}".format(vec))

    def square(self):
        self.x *= self.x
        self.y *= self.y
        self.update()

    def sqrt(self):
        self.x = math.sqrt(self.x)
        self.y = math.sqrt(self.y)
        self.update()

    def distance(self, vec):
        """ Distance between target points of this and another vector. """
        return math.sqrt(math.pow(self.x-vec.x, 2) + math.pow(self.y-vec.y, 2))

    def dot(self, other):
        """ Calculate vector dot product. """
        return (self.x*other.x + self.y*other.y)

    def cross_z(self, other):
        """ Calculate the z component of the cross product. """
        return self.x*other.y - self.y*other.x

    def cross(self, other):
        """ Calculate vector cross product. """
        c = Vector(0, 0, self.cross_z(other))
        return c

    def inverse(self):
        return Vector(-self.x, -self.y)

    def invert(self):
        self.x = -self.x
        self.y = -self.y
        self.update()

    def rotate(self, angle):
        """ Rotate vector in plane by given angle [rad]. """

        # Implementing a general rotation matrix.
        cs = math.cos(angle)
        sn = math.sin(angle)

        self.set(self.x*cs - self.y*sn, self.x*sn + self.y*cs)

    @staticmethod
    def connection(vecFrom, vecTo):
        return Vector2D(vecTo.x-vecFrom.x, vecTo.y-vecFrom.y)

class Line2D:
    def __init__(self):
        # y = mx + n
        self.m = None
        self.n = None

    def __str__(self):
        return "m={}, n={}".format(self.m, self.n)

    def set(self, m, n):
        self.m = m
        self.n = n

    def setFromPoints(self, point0, point1):
        # points are defined by 2D vectors
        x0 = point0.x
        y0 = point0.y
        x1 = point1.x
        y1 = point1.y

        if x0 != x1:
            self.m = (y1-y0) / (x1-x0)
            self.n = y0 - self.m*x0
        else:
            # Vertical line
            self.m = math.inf
            self.n = x0  # Store x intersection in n if line is vertical

    def intersection(self, otherLine):
        m0 = self.m
        n0 = self.n

        m1 = otherLine.m
        n1 = otherLine.n

        if m0 == m1:
            # Lines are parallel
            raise Exception("Lines are parallel.")

        if m0 != math.inf and m1 != math.inf:
            xs = (n1-n0)/(m0-m1)
            ys = m0*xs + n0
            return Vector2D(xs, ys)
        elif m0 == math.inf and m1 != math.inf:
            xs = n0
            ys = m1*xs + n1
            return Vector2D(xs, ys)
        elif m0 != math.inf and m1 == math.inf:
            xs = n1
            ys = m0*xs + n0
            return Vector2D(xs, ys)


class Polygon:
    """ A general polygon with N points in space. """

    def __init__(self, *points):
        """ Points should be defined in counter-clockwise direction. They are of class Vector. """
        self.points = []
        self._area = None
        self.vertexOrderCCW = True  # Vertices defined in counter-clockwise order. Should be False for clockwise.
        self.set(*points)
    
    def __str__(self):
        s = ""
        for i in range(len(self.points)):
            s += "P{i}: ({x}, {y})\n".format(i=i, x=self.points[i].x, y=self.points[i].y)

        return s

    def make3D(self, zComponent):
        """ Convert all points in xy-plane from 2D vectors to 3D vectors,
            using the provided zComponent. """

        for i in range(len(self.points)):
            p = self.points[i]
            newPoint = Vector(p.x, p.y, zComponent)
            self.points[i] = newPoint

    def set(self, *points):
        """ Points should be defined in counter-clockwise direction. They are of class Vector. """
        #print("Setting: {}".format(points))
        #print("points has a length of {}".format(len(points)))

        self.points = []
        self.points.extend(points)

        self._area = None

    def area(self):
        if self._area is None:
            self.calculateArea()

        return self._area

    def calculateArea(self):
        self._area = 0

        # Split polygon into triangles and calculate area of each
        # triangle using the trapezoid method.

        if len(self.points) >= 3:
            # Start at first point
            p1 = self.points[0]
            x1 = p1.x
            y1 = p1.y

            for i in range(1, len(self.points)-1):
                p2 = self.points[i]
                p3 = self.points[i+1]                
                x2 = p2.x
                y2 = p2.y
                x3 = p3.x
                y3 = p3.y
                self._area += 0.5 * ( (y1+y3)*(x3-x1) + (y2+y3)*(x2-x3) - (y1+y2)*(x2-x1) )

    def getBoundingBox(self):
        leftMost  = self.points[0].x
        rightMost = -1
        upMost    = self.points[0].y
        downMost  = -1

        for p in self.points:
            if p.x < leftMost:
                leftMost = math.floor(p.x)
            if p.x > rightMost:
                rightMost = math.ceil(p.x)

            if p.y < upMost:
                upMost = math.floor(p.y)
            if p.y > downMost:
                downMost = math.ceil(p.y)

        return int(leftMost), int(upMost), int(rightMost), int(downMost)

    def isInside2D(self, point):
        """ Check if the given point is inside the polygon or on an edge. """
        x = point.x
        y = point.y

        if len(self.points) >= 3:
            p1 = self.points[0]
            x1 = p1.x
            y1 = p1.y

            # Set up sub-triangles and check if point is in any of those:
            for i in range(1, len(self.points)-1):
                p2 = self.points[i]
                p3 = self.points[i+1]                
                x2 = p2.x
                y2 = p2.y
                x3 = p3.x
                y3 = p3.y

                # Calculate the barycentric coordinates of the point with respect to the triangle:
                D = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3)

                lambda1 = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3)) / D
                lambda2 = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3)) / D
                lambda3 = 1 - lambda1 - lambda2

                #print("Is {} inside triangle?   D: {}, l1: {}, l2: {}, l3: {}".format(point, D, lambda1, lambda2, lambda3))

                if (lambda1>=0 and lambda2>=0 and lambda3>=0):
                    return True

        return False

    def insideEdge(self, edgePoint0, edgePoint1, vertexToTest):
        """ Helper function for clip():
            decide if vertex point is on the "inside" of the clipping edge.
            Inside means "to the left" if vertices are in counter-clockwise direction,
            otherwise "to the right". """

        edge  = edgePoint1 - edgePoint0
        point = vertexToTest - edgePoint0

        cpz = edge.cross_z(point)

        if cpz >= 0:  # on the edge
            return True
        
        return False

    def clip(self, clipPolygon):
        """ Implementation of the Sutherland-Hodgman algorithm.

            clipPolygon must be convex. """

        # Make a list of edges (lines) of the clip polygon:
        outputVertices = self.points # copy.deepcopy(self.points)

        nPoints = len(clipPolygon.points)
        for i in range(nPoints):
            edgePoint0 = clipPolygon.points[i]
            edgePoint1 = clipPolygon.points[((i+1)%nPoints)]

            edgeLine = Line2D()
            edgeLine.setFromPoints(point0=edgePoint0, point1=edgePoint1)

            inputVertices = outputVertices
            outputVertices = []

            #print("## Clip line: {}".format(edgeLine))

            for i in range(len(inputVertices)):
                currentPoint  = inputVertices[i]
                previousPoint = inputVertices[(i+len(inputVertices)-1)%len(inputVertices)]

                #print("  Current Point:  {}".format(currentPoint))
                #print("  Previous Point: {}".format(previousPoint))

                currentLine = Line2D()
                currentLine.setFromPoints(point0=currentPoint, point1=previousPoint)

                #print("  -> current Line: {}".format(currentLine))

                if self.insideEdge(edgePoint0, edgePoint1, currentPoint):
                    #print("  Current point is inside clip polygon.")
                    if not self.insideEdge(edgePoint0, edgePoint1, previousPoint):
                        try:
                            intersectionPoint = currentLine.intersection(edgeLine)
                            #print("  Added intersection to output list.")
                            #print("  -> intersection: {}".format(intersectionPoint))
                            outputVertices.append(intersectionPoint)
                        except:  # Parallel lines -> no intersection
                            pass

                    #print("  Added currentPoint to output list.")
                    outputVertices.append(currentPoint)
                elif self.insideEdge(edgePoint0, edgePoint1, previousPoint):
                    #print("  Current point is not inside clip polygon, but previous point is.")
                    try:
                        intersectionPoint = currentLine.intersection(edgeLine)
                        #print("  Only added intersection point to output list.")
                        #print("  -> intersection: {}".format(intersectionPoint))
                        outputVertices.append(intersectionPoint)
                    except:
                        pass
                #else:
                    #print("Neither current nor previous point is inside clip polygon.")


                #print("\n")

        result = Polygon(*outputVertices)
        return result
