# -*- coding: UTF-8 -*-

import math
import numpy
import copy

from .general import *

# Unit conversions
def inMM(jsonVal):
    """ Convert value/unit pair to mm. """
    if ("value" in jsonVal) and ("unit" in jsonVal):
        value = jsonVal["value"]
        unit  = jsonVal["unit"]

        if(unit == "mm"):
            return value
        elif(unit == "nm"):
            return (value * 1e-6)
        elif(unit == "um"):
            return (value * 1e-3)
        elif(unit == "cm"):
            return (value * 10)
        elif(unit == "dm"):
            return (value * 100)
        elif(unit == "m"):
            return (value * 1000)

        raise Exception(unit + " is not a valid unit of length.")
    else:
        raise Exception("\"value\" or \"unit\" missing.")

def inRad(jsonVal):
    """ Convert value/unit pair to radians. """
    if ("value" in jsonVal) and ("unit" in jsonVal):
        value = jsonVal["value"]
        unit  = jsonVal["unit"]

        if(unit == "rad"):
            return value
        elif(unit == "deg"):
            return ((value * math.pi) / 180.0)

        raise Exception(unit + " is not a valid angular unit.")
    else:
        raise Exception("\"value\" or \"unit\" missing.")

class Matrix:
    def __init__(self, rows, cols):
        self._rows  = rows
        self._cols  = cols

        clist = cols*[0]
        self._value = []
        for i in range(rows):
            self._value.append(clist.copy())

    def __str__(self):
        string = ""
        for r in range(self._rows):
            string += "["
            for c in range(self._cols):
                string += "{space}{x:.7f}  ".format(x=self._value[r][c], space=" "*(self._value[r][c]>=0))

            string += "]\n"

        return string

    def __mul__(self, other):
        if isinstance(other, Vector):
            if self._cols >= 3 and self._rows >= 3:
                result = Vector()
                for r in range(3):
                    s = 0
                    for c in range(3):
                        s += self._value[r][c] * other.get(c)

                    result.setIdx(i=r, value=s)

                return result
            else:
                raise Exception("Error: Cannot multiply matrix of size ({rows} rows x {cols} columns) with a 3-vector.".format(rows=self._rows, cols=self._cols))

class Vector:
    """ A 3D vector in space. """

    def __init__(self, x=0, y=0, z=0):
        self._x = x
        self._y = y
        self._z = z

    def __str__(self):
        return "({spaceX}{x:.7f}, {spaceY}{y:.7f}, {spaceZ}{z:.7f})".format(x=self._x, y=self._y, z=self._z, spaceX=" "*(self._x>=0), spaceY=" "*(self._y>=0), spaceZ=" "*(self._z>=0))

    def __add__(self, other):
        x = self._x + other._x
        y = self._y + other._y
        z = self._z + other._z
        return Vector(x, y, z)

    def __sub__(self, other):
        x = self._x - other._x
        y = self._y - other._y
        z = self._z - other._z
        return Vector(x, y, z)

    def __mul__(self, other):
        x = self._x * other._x
        y = self._y * other._y
        z = self._z * other._z
        return Vector(x, y, z)

    def __truediv__(self, other):
        x = self._x / other._x
        y = self._y / other._y
        z = self._z / other._z
        return Vector(x, y, z)

    def __floordiv__(self, other):
        x = self._x // other._x
        y = self._y // other._y
        z = self._z // other._z
        return Vector(x, y, z)

    def get(self, i):
        if i == 0:
            return self._x
        elif i == 1:
            return self._y
        elif i == 2:
            return self._z

        return 0

    def x(self):
        return self._x

    def setx(self, value):
        self._x = value

    def y(self):
        return self._y

    def sety(self, value):
        self._y = value

    def z(self):
        return self._z

    def setz(self, value):
        self._z = value

    def setxy(self, x=0, y=0):
        """ Set x and y component (relevant for 2D computations) """
        self._x = x
        self._y = y

    def set(self, x=0, y=0, z=0):
        """ Set all vector components. """
        self.setx(x)
        self.sety(y)
        self.setz(z)

    def setIdx(self, i, value):
        if i == 0:
            self._x = value
        elif i == 1:
            self._y = value
        elif i == 2:
            self._z = value

    def length(self):
        """ Calculate length of vector. """
        return math.sqrt(self._x**2 + self._y**2 + self._z**2)

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
        """ Normalize vector lenth to 1. """
        vectorLength = self.length()
        if vectorLength != 0:
            self._x /= float(vectorLength)
            self._y /= float(vectorLength)
            self._z /= float(vectorLength)
        else:
            raise Exception("Unit vector: a zero length vector cannot be converted into a unit vector.")

    def scale(self, factor):
        """ Scale vector length by a certain factor. """
        self._x *= factor
        self._y *= factor
        self._z *= factor

    def add(self, vec):
        """ Add another vector to this vector. """
        self._x += vec._x
        self._y += vec._y
        self._z += vec._z

    def subtract(self, vec):
        """ Subtract another vector from this vector. """
        self._x -= vec._x
        self._y -= vec._y
        self._z -= vec._z

    def square(self):
        self._x *= self._x
        self._y *= self._y
        self._z *= self._z

    def sqrt(self):
        self._x = math.sqrt(self._x)
        self._y = math.sqrt(self._y)
        self._z = math.sqrt(self._z)

    def distance(self, vec):
        """ Distance between target points of this and another vector. """
        return math.sqrt(math.pow(self._x-vec._x, 2) + math.pow(self._y-vec._y, 2) + math.pow(self._z-vec._z, 2))

    def dot(self, other):
        """ Calculate vector dot product. """
        return (self._x*other._x + self._y*other._y + self._z*other._z)

    def cross_z(self, other):
        """ Calculate the z component of the cross product. """
        return self._x*other._y - self._y*other._x

    def cross(self, other):
        """ Calculate vector cross product. """
        cx = self._y*other._z - self._z*other._y
        cy = self._z*other._x - self._x*other._z
        cz = self._x*other._y - self._y*other._x

        c = Vector(cx, cy, cz)
        return c

    def inverse(self):
        return Vector(-self._x, -self._y, -self._z)

    def rotate(self, axis, angleInRad):
        """ Rotate vector around given axis by given angle [rad]. """
        axis.makeUnitVector()

        # Implementing a general rotation matrix.
        cs = math.cos(angleInRad)
        sn = math.sin(angleInRad)

        vx = self._x
        vy = self._y
        vz = self._z

        nx = axis._x
        ny = axis._y
        nz = axis._z

        rx = vx*(nx*nx*(1-cs) + cs)    + vy*(nx*ny*(1-cs) - nz*sn) + vz*(nx*nz*(1-cs) + ny*sn)
        ry = vx*(ny*nx*(1-cs) + nz*sn) + vy*(ny*ny*(1-cs) + cs)    + vz*(ny*nz*(1-cs) - nx*sn)
        rz = vx*(nz*nx*(1-cs) - ny*sn) + vy*(nz*ny*(1-cs) + nx*sn) + vz*(nz*nz*(1-cs) + cs)

        self.set(rx, ry, rz)

    @staticmethod
    def connection(vecFrom, vecTo):
        return Vector(vecTo._x-vecFrom._x, vecTo._y-vecFrom._y, vecTo._z-vecFrom._z)

class Vector2D:
    """ A 2D vector in a plane. """

    def __init__(self, x=0, y=0):
        self._x = x
        self._y = y

    def __str__(self):
        return "({spaceX}{x:.7f}, {spaceY}{y:.7f})".format(x=self._x, y=self._y, spaceX=" "*(self._x>=0), spaceY=" "*(self._y>=0))

    def __add__(self, other):
        x = self._x + other._x
        y = self._y + other._y
        return Vector2D(x, y)

    def __sub__(self, other):
        x = self._x - other._x
        y = self._y - other._y
        return Vector2D(x, y)

    def __mul__(self, other):
        x = self._x * other._x
        y = self._y * other._y
        return Vector2D(x, y)

    def __truediv__(self, other):
        x = self._x / other._x
        y = self._y / other._y
        return Vector2D(x, y)

    def __floordiv__(self, other):
        x = self._x // other._x
        y = self._y // other._y
        return Vector2D(x, y)

    def x(self):
        return self._x

    def setx(self, value):
        self._x = value

    def y(self):
        return self._y

    def sety(self, value):
        self._y = value

    def set(self, x=0, y=0):
        """ Set all vector components. """
        self.setx(x)
        self.sety(y)

    def length(self):
        """ Calculate length of vector. """
        return math.sqrt(self._x**2 + self._y**2)

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
        """ Normalize vector lenth to 1. """
        vectorLength = self.length()
        if vectorLength != 0:
            self._x /= float(vectorLength)
            self._y /= float(vectorLength)
        else:
            raise Exception("Unit vector: a zero length vector cannot be converted into a unit vector.")

    def scale(self, factor):
        """ Scale vector length by a certain factor. """
        self._x *= factor
        self._y *= factor

    def add(self, vec):
        """ Add another vector to this vector. """
        self._x += vec._x
        self._y += vec._y

    def subtract(self, vec):
        """ Subtract another vector from this vector. """
        self._x -= vec._x
        self._y -= vec._y

    def square(self):
        self._x *= self._x
        self._y *= self._y

    def sqrt(self):
        self._x = math.sqrt(self._x)
        self._y = math.sqrt(self._y)

    def distance(self, vec):
        """ Distance between target points of this and another vector. """
        return math.sqrt(math.pow(self._x-vec._x, 2) + math.pow(self._y-vec._y, 2))

    def dot(self, other):
        """ Calculate vector dot product. """
        return (self._x*other._x + self._y*other._y)

    def cross_z(self, other):
        """ Calculate the z component of the cross product. """
        return self._x*other._y - self._y*other._x

    def cross(self, other):
        """ Calculate vector cross product. """
        c = Vector(0, 0, self.cross_z(other))
        return c

    def inverse(self):
        return Vector(-self._x, -self._y)

    def rotate(self, angleInRad):
        """ Rotate vector in plane by given angle [rad]. """

        # Implementing a general rotation matrix.
        cs = math.cos(angleInRad)
        sn = math.sin(angleInRad)

        self.set(self._x*cs - self._y*sn, self._x*sn + self._y*cs)

    @staticmethod
    def connection(vecFrom, vecTo):
        return Vector2D(vecTo._x-vecFrom._x, vecTo._y-vecFrom._y)

class Line2D:
    def __init__(self):
        # y = mx + n
        self._m = None
        self._n = None

    def __str__(self):
        return "m={}, n={}".format(self._m, self._n)

    def set(self, m, n):
        self._m = m
        self._n = n

    def setFromPoints(self, point0, point1):
        # points are defined by 2D vectors
        x0 = point0._x
        y0 = point0._y
        x1 = point1._x
        y1 = point1._y

        if x0 != x1:
            self._m = (y1-y0) / (x1-x0)
            self._n = y0 - self._m*x0
        else:
            # Vertical line
            self._m = math.inf
            self._n = x0  # Store x intersection in n if line is vertical

    def intersection(self, otherLine):
        m0 = self._m
        n0 = self._n

        m1 = otherLine._m
        n1 = otherLine._n

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
        self._points = []
        self._area = None
        self._vertexOrderCCW = True  # Vertices defined in counter-clockwise order. Should be False for clockwise.
        self.set(*points)
    
    def __str__(self):
        s = ""
        for i in range(len(self._points)):
            s += "P{i}: ({x}, {y})\n".format(i=i, x=self._points[i]._x, y=self._points[i]._y)

        return s

    def make3D(self, zComponent):
        """ Convert all points in xy-plane from 2D vectors to 3D vectors,
            using the provided zComponent. """

        for i in range(len(self._points)):
            p = self._points[i]
            newPoint = Vector(p._x, p._y, zComponent)
            self._points[i] = newPoint

    def set(self, *points):
        """ Points should be defined in counter-clockwise direction. They are of class Vector. """
        #print("Setting: {}".format(points))
        #print("points has a length of {}".format(len(points)))

        self._points = []
        self._points.extend(points)

        self._area = None

    def area(self):
        if self._area == None:
            self.calculateArea()

        return self._area

    def calculateArea(self):
        self._area = 0

        # Split polygon into triangles and calculate area of each
        # triangle using the trapezoid method.

        if len(self._points) >= 3:
            # Start at first point
            p1 = self._points[0]
            x1 = p1._x
            y1 = p1._y

            for i in range(1, len(self._points)-1):
                p2 = self._points[i]
                p3 = self._points[i+1]                
                x2 = p2._x
                y2 = p2._y
                x3 = p3._x
                y3 = p3._y
                self._area += 0.5 * ( (y1+y3)*(x3-x1) + (y2+y3)*(x2-x3) - (y1+y2)*(x2-x1) )

    def getBoundingBox(self):
        leftMost  = self._points[0]._x
        rightMost = -1
        upMost    = self._points[0]._y
        downMost  = -1

        for p in self._points:
            if p._x < leftMost:
                leftMost = math.floor(p._x)
            if p._x > rightMost:
                rightMost = math.ceil(p._x)

            if p._y < upMost:
                upMost = math.floor(p._y)
            if p._y > downMost:
                downMost = math.ceil(p._y)

        return int(leftMost), int(upMost), int(rightMost), int(downMost)

    def isInside2D(self, point):
        """ Check if the given point is inside the polygon or on an edge. """
        x = point._x
        y = point._y

        if len(self._points) >= 3:
            p1 = self._points[0]
            x1 = p1._x
            y1 = p1._y

            # Set up sub-triangles and check if point is in any of those:
            for i in range(1, len(self._points)-1):
                p2 = self._points[i]
                p3 = self._points[i+1]                
                x2 = p2._x
                y2 = p2._y
                x3 = p3._x
                y3 = p3._y

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
        outputVertices = self._points # copy.deepcopy(self._points)

        nPoints = len(clipPolygon._points)
        for i in range(nPoints):
            edgePoint0 = clipPolygon._points[i]
            edgePoint1 = clipPolygon._points[((i+1)%nPoints)]

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
