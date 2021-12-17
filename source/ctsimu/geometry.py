# -*- coding: UTF-8 -*-

import numpy
import os    # File and path handling
import json
import math
import copy
import pkgutil

from .geoprimitives import *
from .image import Image  # To create detector flat field

def basisTransformMatrix(fromCS, toCS):
    T = Matrix(3, 3)

    # Row 1:
    T._value[0][0] = toCS._vectorU.dot(fromCS._vectorU)
    T._value[0][1] = toCS._vectorU.dot(fromCS._vectorV)
    T._value[0][2] = toCS._vectorU.dot(fromCS._vectorW)

    # Row 2:
    T._value[1][0] = toCS._vectorV.dot(fromCS._vectorU)
    T._value[1][1] = toCS._vectorV.dot(fromCS._vectorV)
    T._value[1][2] = toCS._vectorV.dot(fromCS._vectorW)

    # Row 3:
    T._value[2][0] = toCS._vectorW.dot(fromCS._vectorU)
    T._value[2][1] = toCS._vectorW.dot(fromCS._vectorV)
    T._value[2][2] = toCS._vectorW.dot(fromCS._vectorW)

    return T

class GeometryObject:
    """ An object according to the CTSimU scenario specification,
        containing a centre coordinate and an orientation in 3D space. 
        So far, geometrical objects can be source, stage or detector.
        Samples would need additional attention due to
        possible stage coordinate system (instead of world)."""

    def __init__(self):
        """ Initialize all vectors as world CS: """
        self._centre  = Vector(0, 0, 0)
        self._vectorU = Vector(1, 0, 0)
        self._vectorV = Vector(0, 1, 0)
        self._vectorW = Vector(0, 0, 1)

    def setupFromGeometryDefinition(self, geometry):
        """ Set up this geometrical object from a JSON dictionary. """
        if "centre" in geometry:
            if "x" in geometry["centre"]:
                cx = inMM(geometry["centre"]["x"])
            else:
                raise Exception("No \"x\" coordinate found for the center.")

            if "y" in geometry["centre"]:
                cy = inMM(geometry["centre"]["y"])
            else:
                raise Exception("No \"y\" coordinate found for the center.")
            
            if "z" in geometry["centre"]:
                cz = inMM(geometry["centre"]["z"])
            else:
                raise Exception("No \"z\" coordinate found for the center.")

            self._centre.set(cx, cy, cz)
        else:
            raise Exception("JSON file is missing a geometry \"centre\" section.")


        if "vector_u" in geometry:
            if "x" in geometry["vector_u"]:
                ux = geometry["vector_u"]["x"]
            else:
                raise Exception("No \"x\" component found for vector u.")

            if "y" in geometry["vector_u"]:
                uy = geometry["vector_u"]["y"]
            else:
                raise Exception("No \"y\" component found for vector u.")

            if "z" in geometry["vector_u"]:
                uz = geometry["vector_u"]["z"]
            else:
                raise Exception("No \"z\" component found for vector u.")
        else:
            raise Exception("JSON file is missing a geometry \"vector_u\" section.")

        
        if "vector_w" in geometry:
            if "x" in geometry["vector_w"]:
                wx = geometry["vector_w"]["x"]
            else:
                raise Exception("No \"x\" component found for vector w.")

            if "y" in geometry["vector_w"]:
                wy = geometry["vector_w"]["y"]
            else:
                raise Exception("No \"y\" component found for vector w.")

            if "z" in geometry["vector_w"]:
                wz = geometry["vector_w"]["z"]
            else:
                raise Exception("No \"z\" component found for vector w.")
        else:
            raise Exception("JSON file is missing a geometry \"vector_w\" section.")


        self.setup(cx, cy, cz, ux, uy, uz, wx, wy, wz)

        if "deviation" in geometry:
            if "position" in geometry["deviation"]:
                if "x" in geometry["deviation"]["position"]:
                    if geometry["deviation"]["position"]["x"] != None:
                        translationX = inMM(geometry["deviation"]["position"]["x"])
                        self.translateX(translationX)
    
                if "y" in geometry["deviation"]["position"]:
                    if geometry["deviation"]["position"]["y"] != None:
                        translationY = inMM(geometry["deviation"]["position"]["y"])
                        self.translateY(translationY)

                if "z" in geometry["deviation"]["position"]:
                    if geometry["deviation"]["position"]["z"] != None:
                        translationZ = inMM(geometry["deviation"]["position"]["z"])
                        self.translateZ(translationZ)

            # Rotations according to w''v'u convention:
            if "rotation" in geometry["deviation"]:
                if "w" in geometry["deviation"]["rotation"]:
                    if geometry["deviation"]["rotation"]["w"] != None:
                        angleAroundW = inRad(geometry["deviation"]["rotation"]["w"])
                        self.rotateAroundW(angleAroundW)

                if "v" in geometry["deviation"]["rotation"]:
                    if geometry["deviation"]["rotation"]["v"] != None:
                        angleAroundV = inRad(geometry["deviation"]["rotation"]["v"])
                        self.rotateAroundV(angleAroundV)

                if "u" in geometry["deviation"]["rotation"]:
                    if geometry["deviation"]["rotation"]["u"] != None:
                        angleAroundU = inRad(geometry["deviation"]["rotation"]["u"])
                        self.rotateAroundU(angleAroundU)

    def setup(self, centerX=0, centerY=0, centerZ=0, uX=0, uY=0, uZ=0, wX=0, wY=0, wZ=0):
        """ Set up centre and orientation from JSON geometry components. """
        self._centre  = Vector(centerX, centerY, centerZ)
        self._vectorU = Vector(uX, uY, uZ)
        self._vectorW = Vector(wX, wY, wZ)
        
        self._vectorV = self._vectorW.cross(self._vectorU)

        self.makeUnitCS()

    def makeUnitCS(self):
        """ Make coordinate system base unit vectors. """
        self._vectorU.makeUnitVector()
        self._vectorV.makeUnitVector()
        self._vectorW.makeUnitVector()

    def translate(self, translationVector):
        """ Move object in space. """
        self._centre.add(translationVector)

    def translateX(self, dx):
        """ Move object in x direction. """
        self._centre.setx(self._centre.x() + dx)

    def translateY(self, dy):
        """ Move object in y direction. """
        self._centre.sety(self._centre.y() + dy)

    def translateZ(self, dz):
        """ Move object in z direction. """
        self._centre.setz(self._centre.z() + dz)

    def rotateAroundU(self, angleInRad):
        """ Rotate object around its u axis by given angle [rad]. """
        self._vectorV.rotate(self._vectorU, angleInRad)
        self._vectorW.rotate(self._vectorU, angleInRad)

    def rotateAroundV(self, angleInRad):
        """ Rotate object around its v axis by given angle [rad]. """
        self._vectorU.rotate(self._vectorV, angleInRad)
        self._vectorW.rotate(self._vectorV, angleInRad)

    def rotateAroundW(self, angleInRad):
        """ Rotate object around its w axis by given angle [rad]. """
        self._vectorU.rotate(self._vectorW, angleInRad)
        self._vectorV.rotate(self._vectorW, angleInRad)

    def rotate(self, axis, angleInRad):
        """ Rotate coordinate system around a given axis by angle [rad]. """
        self._vectorU.rotate(axis, angleInRad)
        self._vectorV.rotate(axis, angleInRad)
        self._vectorW.rotate(axis, angleInRad)

    def changeReferenceFrame(self, fromCS, toCS):
        # Rotate basis vectors into toCS:
        T = basisTransformMatrix(fromCS, toCS)
        self._vectorU = T * self._vectorU
        self._vectorV = T * self._vectorV
        self._vectorW = T * self._vectorW

        # Move center to toCS:
        centreDiff = fromCS._centre - toCS._centre
        newRelativeCentreInFrom = self._centre + centreDiff
        self._centre = T * newRelativeCentreInFrom


class Detector(GeometryObject):
    """ Detector class to get pixel coordinates, etc. """

    def __init__(self):
        GeometryObject.__init__(self)

        self._pixelsU     = None  # Detector pixels in u direction
        self._pixelsV     = None  # Detector pixels in v direction
        self._pixelPitchU = None
        self._pixelPitchV = None
        self._physWidth    = 0    # Physical width in units of pitch U
        self._physHeight   = 0    # Physical height in units of pitch V

        self._upperLeftX = 0
        self._upperLeftY = 0
        self._upperLeftZ = 0

    def setSize(self, nPixelsU=None, nPixelsV=None, pitchU=None, pitchV=None):
        self._pixelsU = nPixelsU
        self._pixelsV = nPixelsV
        self._pixelPitchU = pitchU
        self._pixelPitchV = pitchV

        self.computeGeometryParameters()

    def computeGeometryParameters(self):
        if not(self._pixelsU==None or self._pixelsV==None or self._pixelPitchU==None or self._pixelPitchV==None):
            # Physical width and height:
            self._physWidth  = self._pixelsU * self._pixelPitchU
            self._physHeight = self._pixelsV * self._pixelPitchV

            # Vectors of the detector coordinate system:
            ux = self._vectorU.x()
            uy = self._vectorU.y()
            uz = self._vectorU.z()
            vx = self._vectorV.x()
            vy = self._vectorV.y()
            vz = self._vectorV.z()

            # World coordinates of corner (0, 0) of detector CS:
            self._upperLeftX = self._centre.x() - 0.5*(ux*self._physWidth + vx*self._physHeight)
            self._upperLeftY = self._centre.y() - 0.5*(uy*self._physWidth + vy*self._physHeight)
            self._upperLeftZ = self._centre.z() - 0.5*(uz*self._physWidth + vz*self._physHeight)

    def cols(self):
        return self._pixelsU

    def rows(self):
        return self._pixelsV

    def physWidth(self):
        return self._physWidth

    def physHeight(self):
        return self._physHeight

    def pitchU(self):
        return self._pixelPitchU

    def pitchV(self):
        return self._pixelPitchV

    def pixelVectorUpperLeft(self, x, y):
        # x, y are coordinates in pixel coordinates system
        px = self._upperLeftX + self._vectorU.x()*x*self._pixelPitchU + self._vectorV.x()*y*self._pixelPitchV
        py = self._upperLeftY + self._vectorU.y()*x*self._pixelPitchU + self._vectorV.y()*y*self._pixelPitchV
        pz = self._upperLeftZ + self._vectorU.z()*x*self._pixelPitchU + self._vectorV.z()*y*self._pixelPitchV
        pixelVector = Vector(px, py, pz)
        return pixelVector

    def pixelVectorCenter(self, x, y):
        return self.pixelVectorUpperLeft(x+0.5, y+0.5)


class Geometry:
    """ Gathers information about the CT setup (source and detector)
        that is necessary for an analytical flat field correction.
        Calculates the analytical intensity distribution (flat field). """

    def __init__(self, jsonFile=None, jsonFileFromPkg=None):
        """ Initialize using the provided JSON geometry specification. """
        self._detector    = Detector()
        self._source      = GeometryObject()

        self._detectorTotalTilt = None
        self._SDD = None
        self._brightestSpotWorld = None
        self._brightestSpotDetector = None

        jsonText = None
        if jsonFileFromPkg != None:  # from package
            jsonFile = jsonFileFromPkg
            jsonText = pkgutil.get_data(__name__, jsonFileFromPkg).decode()
        elif jsonFile != None:  # from file
            if os.path.isfile(jsonFile):
                log("JSON File: {}".format(jsonFile))
                jsonFilePtr = open(jsonFile, "r")
                jsonText = jsonFilePtr.read()
                jsonFilePtr.close()
            else:
                raise Exception("Can't find " + jsonFile)
        else:
            raise Exception("Please provide a valid JSON file.")

        if(jsonText != None):
            try:
                jsonDict = json.loads(jsonText)
            except:
                raise Exception("Error parsing JSON file: {}".format(jsonFile))

            # Detector size and pixel pitch:
            pixelsU = getFieldOrNone(jsonDict, "detector", "columns", "value")
            pixelsV = getFieldOrNone(jsonDict, "detector", "rows", "value")
            pixelPitchU = getFieldOrNone(jsonDict, "detector", "pixel_pitch", "u", "value")
            pixelPitchV = getFieldOrNone(jsonDict, "detector", "pixel_pitch", "v", "value")

            try:
                detectorGeometry = getFieldOrNone(jsonDict, "geometry", "detector")
                if detectorGeometry != None:
                    self._detector.setupFromGeometryDefinition(detectorGeometry)
                    self._detector.setSize(pixelsU, pixelsV, pixelPitchU, pixelPitchV)
                else:
                    raise Exception("JSON file does not contain a valid \"detector\" section in \"geometry\".")
            except Exception as e:
                log("Something went wrong when setting up the detector geometry using the JSON file description.")
                raise Exception(e)

            try:
                sourceGeometry = getFieldOrNone(jsonDict, "geometry", "source")
                if sourceGeometry != None:
                    self._source.setupFromGeometryDefinition(sourceGeometry)
                else:
                    raise Exception("JSON file does not contain a valid \"source\" section in \"geometry\".")
            except Exception as e:
                log("Something went wrong when setting up the source geometry using the JSON file description.")
                raise Exception(e)

            # Calculate additional parameters:

            ## Detector Tilt Angle:
            sourceDetectorCentreConnection = copy.deepcopy(self._detector._centre)
            sourceDetectorCentreConnection.subtract(self._source._centre)

            self._detectorTotalTilt = self._detector._vectorW.angle(sourceDetectorCentreConnection)

            ## Distance between source centre and detector centre:
            dist = self._detector._centre.distance(self._source._centre)

            self._SDD = dist * math.cos(self._detectorTotalTilt)

            ## Brightest Spot in World Coordinate System:
            self._brightestSpotWorld = copy.deepcopy(self._detector._vectorW)
            self._brightestSpotWorld.scale(self._SDD)
            self._brightestSpotWorld.add(self._source._centre)

            ## Brightest Spot in Detector Coordinate System:
            self._brightestSpotDetector = copy.deepcopy(self._brightestSpotWorld)
            self._brightestSpotDetector.subtract(self._detector._centre)
            
            pxU = self._brightestSpotDetector.dot(self._detector._vectorU) / self._detector.pitchU() + self._detector.cols()/2.0
            pxV = self._brightestSpotDetector.dot(self._detector._vectorV) / self._detector.pitchV() + self._detector.rows()/2.0

            self._brightestSpotDetector = Vector(pxU, pxV, 0)
        else:
            raise Exception("JSON scenario file not available.")
        

    def info(self):
        txt  = "Detector\n"
        txt += "===========================================================\n"
        txt += "Center:          {}\n".format(self._detector._centre)
        txt += "u:               {}\n".format(self._detector._vectorU)
        txt += "v:               {}\n".format(self._detector._vectorV)
        txt += "w:               {}\n".format(self._detector._vectorW)
        txt += "Pixels:          {cols} x {rows}\n".format(cols=self._detector.cols(), rows=self._detector.rows())
        txt += "Pitch:           {pitchU} mm x {pitchV} mm\n".format(pitchU=self._detector.pitchU(), pitchV=self._detector.pitchV())
        txt += "Physical Size:   {width} mm x {height} mm\n".format(width=self._detector.physWidth(), height=self._detector.physHeight())
        txt += "Total Tilt:      {tiltRad} rad = {tiltDeg} deg\n".format(tiltRad=self._detectorTotalTilt, tiltDeg=180.0*self._detectorTotalTilt/math.pi)
        txt += "Center Distance: {} mm\n".format(self._detector._centre.distance(self._source._centre))

        # Source - Detector distance (SDD) defined by shortest distance between source and detector:
        txt += "SDD:             {} mm\n".format(self._SDD)
        txt += "Brightest Spot:\n"
        txt += "  World:         {}\n".format(self._brightestSpotWorld)
        txt += "  Pixels:        {}\n".format(self._brightestSpotDetector)

        txt += "\n"
        txt += "Source:\n"
        txt += "===========================================================\n"
        txt += "Center:          {}\n".format(self._source._centre)
        txt += "u:               {}\n".format(self._source._vectorU)
        txt += "v:               {}\n".format(self._source._vectorV)
        txt += "w:               {}\n".format(self._source._vectorW)

        return txt

    def createDetectorFlatField_rays(self):
        """ Calculate an analytical free beam intensity distribution
            picture for the given detector, to be used for an
            ideal flat field correction. """
        width      = self._detector.cols()
        height     = self._detector.rows()
        pixelSizeU = self._detector.pitchU()
        pixelSizeV = self._detector.pitchV()

        if(width == None):
            raise Exception("The detector width (in pixels) must be provided through a valid CTSimU JSON file.")
        if(height == None):
            raise Exception("The detector height (in pixels) must be provided through a valid CTSimU JSON file.")
        if(pixelSizeU == None):
            raise Exception("The pixel size (in mm) in u direction must be provided through a valid CTSimU JSON file.")
        if(pixelSizeV == None):
            raise Exception("The pixel size (in mm) in v direction must be provided through a valid CTSimU JSON file.")

        flatField = Image()
        flatField.shape(width, height, 0, flatField.getInternalDataType())

        # Positions of detector and source center:
        dx = self._detector._centre.x()
        dy = self._detector._centre.y()
        dz = self._detector._centre.z()

        sx = self._source._centre.x()
        sy = self._source._centre.y()
        sz = self._source._centre.z()

        # Vectors of the detector coordinate system:
        ux = self._detector._vectorU.x()
        uy = self._detector._vectorU.y()
        uz = self._detector._vectorU.z()

        vx = self._detector._vectorV.x()
        vy = self._detector._vectorV.y()
        vz = self._detector._vectorV.z()

        wx = self._detector._vectorW.x()
        wy = self._detector._vectorW.y()
        wz = self._detector._vectorW.z()


        # Angle 'alpha' between detector normal and connection line [detector centre -- source]:
        connectionLine = Vector(dx-sx, dy-sy, dz-sz)

        alpha = abs(self._detector._vectorW.angle(connectionLine))
        if alpha > (math.pi/2):
            alpha = math.pi - alpha

        # Distance between source centre and detector centre:
        dist = self._detector._centre.distance(self._source._centre)

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
                        pixel = self._detector.pixelVectorUpperLeft(x+(gx+1)*stepSize, y+(gy+1)*stepSize)

                        # Grid with no margin:
                        #if gridSize > 1:
                        #    stepSize = 1.0 / (gridSize-1)
                        #    pixel = self._detector.pixelVectorUpperLeft(x+gx*stepSize, y+gy*stepSize)
                        #else:
                        #    pixel = self._detector.pixelVectorCenter(x, y)

                        distToSource = self._source._centre.distance(pixel)

                        # Angle of incident rays:
                        vecSourceToPixel = Vector(pixel.x()-sx, pixel.y()-sy, pixel.z()-sz)
                        incidenceAngle = abs(self._detector._vectorW.angle(vecSourceToPixel))
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

        if len(polygon._points) >= 3:
            # Start at first point
            p1 = polygon._points[0]

            area = 0

            for i in range(1, len(polygon._points)-1):
                p2 = polygon._points[i]
                p3 = polygon._points[i+1]

                area += self.triangleAreaOnUnitSphere(p1, p2, p3)

            return area
        else:
            return 0

    """
    def createDetectorFlatField_sphere_old(self, clippingPolygon=None):
        # Positions of detector and source center:
        dx = self._detector._centre.x()
        dy = self._detector._centre.y()
        dz = self._detector._centre.z()

        sx = self._source._centre.x()
        sy = self._source._centre.y()
        sz = self._source._centre.z()

        # Angle 'alpha' between detector normal and connection line [detector centre -- source]:
        connectionLine = Vector(dx-sx, dy-sy, dz-sz)

        alpha = abs(self._detector._vectorW.angle(connectionLine))
        if alpha > (math.pi/2):
            alpha = math.pi - alpha

        # Distance between source centre and detector centre:
        dist = self._detector._centre.distance(self._source._centre)

        # Source - Detector distance (SDD) defined by shortest distance between source and detector,
        # or distance between source and spot of highest intensity on detector.
        SDD = dist * math.cos(alpha)

        # Create a new detector in a coordinate system where source is at (0, 0, 0):
        det = copy.deepcopy(self._detector)
        translationVector = Vector(-sx, -sy, -sz)
        det.translate(translationVector)
        det.computeGeometryParameters()

        # Calculate the area of the theoretical "brightest pixel" on the unit sphere:
        hpu = 0.5*det.pitchU()
        hpv = 0.5*det.pitchV()
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
                A = det.pixelVectorUpperLeft(x,   y)
                B = det.pixelVectorUpperLeft(x+1, y)
                C = det.pixelVectorUpperLeft(x+1, y+1)
                D = det.pixelVectorUpperLeft(x,   y+1)

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
        incidenceAngle = abs(self._detector._vectorW.angle(maxCenter))

        #print("Brightest Pixel: {}, {}".format(maxX, maxY))
        print("  Vector: {}, {}, {}".format(maxCenter.x(), maxCenter.y(), maxCenter.z()))
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
        D = copy.deepcopy(self._detector)
        S = copy.deepcopy(self._source)
        world = GeometryObject()  # will be initialized as world

        S.changeReferenceFrame(world, D)
        D.changeReferenceFrame(world, D)
        D.computeGeometryParameters()

        # Source - Detector distance (SDD) defined by shortest distance between source and detector,
        # or distance between source and spot of highest intensity on detector.
        SDD = abs(S._centre.z())

        # Calculate the area of the theoretical "brightest pixel" on the unit sphere:
        pu = D.pitchU()
        pv = D.pitchV()
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
        p00 = D.pixelVectorUpperLeft(0, 0)

        stepRight      = Vector(pu, 0,  0)
        stepDown       = Vector(0,  pv, 0)
        stepRightDown  = Vector(pu, pv, 0)

        # Move the clipping polygon to a coordinate system
        # where source is centered:
        for coverPolygon in coverPolygons:
            for p in range(len(coverPolygon._points)):
                coverPolygon._points[p] = coverPolygon._points[p] - S._centre

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
                pA = pA - S._centre
                pB = pB - S._centre
                pC = pC - S._centre
                pD = pD - S._centre

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

        width  = self._detector.cols()
        height = self._detector.rows()
        pitchU = self._detector.pitchU()
        pitchV = self._detector.pitchV()

        if(width == None):
            raise Exception("The detector width (in pixels) must be provided through a valid CTSimU JSON file.")
        if(height == None):
            raise Exception("The detector height (in pixels) must be provided through a valid CTSimU JSON file.")
        if(pitchU == None):
            raise Exception("The pixel size (in mm) in u direction must be provided through a valid CTSimU JSON file.")
        if(pitchV == None):
            raise Exception("The pixel size (in mm) in v direction must be provided through a valid CTSimU JSON file.")

        flatField = Image()
        flatField.shape(width, height, 0, flatField.getInternalDataType())

        # Positions of detector and source center:
        dx = self._detector._centre.x()
        dy = self._detector._centre.y()
        dz = self._detector._centre.z()

        sx = self._source._centre.x()
        sy = self._source._centre.y()
        sz = self._source._centre.z()

        # Vectors of the detector coordinate system:
        ux = self._detector._vectorU.x()
        uy = self._detector._vectorU.y()
        uz = self._detector._vectorU.z()

        vx = self._detector._vectorV.x()
        vy = self._detector._vectorV.y()
        vz = self._detector._vectorV.z()

        wx = self._detector._vectorW.x()
        wy = self._detector._vectorW.y()
        wz = self._detector._vectorW.z()


        # Angle 'alpha' between detector normal and connection line [detector centre -- source]:
        connectionLine = Vector(dx-sx, dy-sy, dz-sz)

        alpha = abs(self._detector._vectorW.angle(connectionLine))
        if alpha > (math.pi/2):
            alpha = math.pi - alpha

        # Distance between source centre and detector centre:
        dist = self._detector._centre.distance(self._source._centre)

        # Source - Detector distance (SDD) defined by shortest distance between source and detector:
        SDD = dist * math.cos(alpha)

        maxIntensity = 0
        maxX = 0
        maxY = 0

        # Create a new detector in a coordinate system where source is at (0, 0, 0):
        det = copy.deepcopy(self._detector)
        translationVector = Vector(-sx, -sy, -sz)
        det.translate(translationVector)
        det.computeGeometryParameters()

        upperLeft_u = det.pixelVectorUpperLeft(0, 0).dot(self._detector._vectorU)
        upperLeft_v = det.pixelVectorUpperLeft(0, 0).dot(self._detector._vectorV)
        upperLeft_w = det.pixelVectorUpperLeft(0, 0).dot(self._detector._vectorW)

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
        hpu = 0.5*det.pitchU()
        hpv = 0.5*det.pitchV()
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