# -*- coding: UTF-8 -*-

import numpy
import os    # File and path handling
import json
import math
import copy
import pkgutil
from datetime import datetime

from .geoprimitives import *
from .image import Image  # To create detector flat field

def basisTransformMatrix(fromCS, toCS):
    T = Matrix(3, 3)

    # Row 1:
    T.value[0][0] = toCS.vectorU.dot(fromCS.vectorU)
    T.value[0][1] = toCS.vectorU.dot(fromCS.vectorV)
    T.value[0][2] = toCS.vectorU.dot(fromCS.vectorW)

    # Row 2:
    T.value[1][0] = toCS.vectorV.dot(fromCS.vectorU)
    T.value[1][1] = toCS.vectorV.dot(fromCS.vectorV)
    T.value[1][2] = toCS.vectorV.dot(fromCS.vectorW)

    # Row 3:
    T.value[2][0] = toCS.vectorW.dot(fromCS.vectorU)
    T.value[2][1] = toCS.vectorW.dot(fromCS.vectorV)
    T.value[2][2] = toCS.vectorW.dot(fromCS.vectorW)

    return T

class GeometryObject:
    """ An object according to the CTSimU scenario specification,
        containing a centre coordinate and an orientation in 3D space. 
        So far, geometrical objects can be source, stage or detector.
        Samples would need additional attention due to
        possible stage coordinate system (instead of world)."""

    def __init__(self):
        """ Initialize all vectors as world CS: """
        self.centre  = Vector(0, 0, 0)
        self.vectorU = Vector(1, 0, 0)
        self.vectorV = Vector(0, 1, 0)
        self.vectorW = Vector(0, 0, 1)

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

            self.centre.set(cx, cy, cz)
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
        self.centre  = Vector(centerX, centerY, centerZ)
        self.vectorU = Vector(uX, uY, uZ)
        self.vectorW = Vector(wX, wY, wZ)
        
        self.vectorV = self.vectorW.cross(self.vectorU)

        self.makeUnitCS()

    def makeUnitCS(self):
        """ Make coordinate system base unit vectors. """
        self.vectorU.makeUnitVector()
        self.vectorV.makeUnitVector()
        self.vectorW.makeUnitVector()

    def translate(self, translationVector):
        """ Move object in space. """
        self.centre.add(translationVector)

    def translateX(self, dx):
        """ Move object in x direction. """
        self.centre.setx(self.centre.x() + dx)

    def translateY(self, dy):
        """ Move object in y direction. """
        self.centre.sety(self.centre.y() + dy)

    def translateZ(self, dz):
        """ Move object in z direction. """
        self.centre.setz(self.centre.z() + dz)

    def rotateAroundU(self, angleInRad):
        """ Rotate object around its u axis by given angle [rad]. """
        self.vectorV.rotate(self.vectorU, angleInRad)
        self.vectorW.rotate(self.vectorU, angleInRad)

    def rotateAroundV(self, angleInRad):
        """ Rotate object around its v axis by given angle [rad]. """
        self.vectorU.rotate(self.vectorV, angleInRad)
        self.vectorW.rotate(self.vectorV, angleInRad)

    def rotateAroundW(self, angleInRad):
        """ Rotate object around its w axis by given angle [rad]. """
        self.vectorU.rotate(self.vectorW, angleInRad)
        self.vectorV.rotate(self.vectorW, angleInRad)

    def rotate(self, axis, angleInRad):
        """ Rotate coordinate system around a given axis by angle [rad]. """
        self.vectorU.rotate(axis, angleInRad)
        self.vectorV.rotate(axis, angleInRad)
        self.vectorW.rotate(axis, angleInRad)

    def changeReferenceFrame(self, fromCS, toCS):
        # Rotate basis vectors into toCS:
        T = basisTransformMatrix(fromCS, toCS)
        self.vectorU = T * self.vectorU
        self.vectorV = T * self.vectorV
        self.vectorW = T * self.vectorW

        # Move center to toCS:
        centreDiff = fromCS.centre - toCS.centre
        newRelativeCentreInFrom = self.centre + centreDiff
        self.centre = T * newRelativeCentreInFrom


class Detector(GeometryObject):
    """ Detector class to get pixel coordinates, etc. """

    def __init__(self):
        GeometryObject.__init__(self)

        self.pixelsU     = None  # Detector pixels in u direction
        self.pixelsV     = None  # Detector pixels in v direction
        self.pixelPitchU = None
        self.pixelPitchV = None
        self.physWidth    = 0    # Physical width in units of pitch U
        self.physHeight   = 0    # Physical height in units of pitch V

        self.upperLeftX = 0
        self.upperLeftY = 0
        self.upperLeftZ = 0

    def setSize(self, nPixelsU=None, nPixelsV=None, pitchU=None, pitchV=None):
        self.pixelsU = nPixelsU
        self.pixelsV = nPixelsV
        self.pixelPitchU = pitchU
        self.pixelPitchV = pitchV

        self.computeGeometryParameters()

    def computeGeometryParameters(self):
        if not(self.pixelsU is None or self.pixelsV is None or self.pixelPitchU is None or self.pixelPitchV is None):
            # Physical width and height:
            self.physWidth  = self.pixelsU * self.pixelPitchU
            self.physHeight = self.pixelsV * self.pixelPitchV

            # Vectors of the detector coordinate system:
            ux = self.vectorU.x()
            uy = self.vectorU.y()
            uz = self.vectorU.z()
            vx = self.vectorV.x()
            vy = self.vectorV.y()
            vz = self.vectorV.z()

            # World coordinates of corner (0, 0) of detector CS:
            self.upperLeftX = self.centre.x() - 0.5*(ux*self.physWidth + vx*self.physHeight)
            self.upperLeftY = self.centre.y() - 0.5*(uy*self.physWidth + vy*self.physHeight)
            self.upperLeftZ = self.centre.z() - 0.5*(uz*self.physWidth + vz*self.physHeight)

    def cols(self):
        return self.pixelsU

    def rows(self):
        return self.pixelsV

    def physWidth(self):
        return self.physWidth

    def physHeight(self):
        return self.physHeight

    def pitchU(self):
        return self.pixelPitchU

    def pitchV(self):
        return self.pixelPitchV

    def pixelVectorUpperLeft(self, x, y):
        # x, y are coordinates in pixel coordinates system
        px = self.upperLeftX + self.vectorU.x()*x*self.pixelPitchU + self.vectorV.x()*y*self.pixelPitchV
        py = self.upperLeftY + self.vectorU.y()*x*self.pixelPitchU + self.vectorV.y()*y*self.pixelPitchV
        pz = self.upperLeftZ + self.vectorU.z()*x*self.pixelPitchU + self.vectorV.z()*y*self.pixelPitchV
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
        self.detector    = Detector()
        self.source      = GeometryObject()
        self.stage       = GeometryObject()

        self.detectorTotalTilt = None
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
            pixelPitchU = getFieldOrNone(jsonDict, "detector", "pixel_pitch", "u", "value")
            pixelPitchV = getFieldOrNone(jsonDict, "detector", "pixel_pitch", "v", "value")

            try:
                detectorGeometry = getFieldOrNone(jsonDict, "geometry", "detector")
                if detectorGeometry != None:
                    self.detector.setupFromGeometryDefinition(detectorGeometry)
                    self.detector.setSize(pixelsU, pixelsV, pixelPitchU, pixelPitchV)
                else:
                    raise Exception("JSON file does not contain a valid \"detector\" section in \"geometry\".")
            except Exception as e:
                log("Something went wrong when setting up the detector geometry using the JSON file description.")
                raise Exception(e)

            try:
                sourceGeometry = getFieldOrNone(jsonDict, "geometry", "source")
                if sourceGeometry != None:
                    self.source.setupFromGeometryDefinition(sourceGeometry)
                else:
                    raise Exception("JSON file does not contain a valid \"source\" section in \"geometry\".")
            except Exception as e:
                log("Something went wrong when setting up the source geometry using the JSON file description.")
                raise Exception(e)

            try:
                stageGeometry = getFieldOrNone(jsonDict, "geometry", "stage")
                if stageGeometry != None:
                    self.stage.setupFromGeometryDefinition(stageGeometry)
                else:
                    raise Exception("JSON file does not contain a valid \"stage\" section in \"geometry\".")
            except Exception as e:
                log("Something went wrong when setting up the stage geometry using the JSON file description.")
                raise Exception(e)

            self.update()            
        else:
            raise Exception("JSON scenario file not available.")
        
    def update(self):
        # Calculate additional parameters:

        # SOD, SDD, ODD
        world = GeometryObject()
        source_from_detector = copy.deepcopy(self.source)
        stage_from_detector  = copy.deepcopy(self.stage)

        source_from_detector.changeReferenceFrame(world, self.detector)
        stage_from_detector.changeReferenceFrame(world, self.detector)

        self.SDD = abs(source_from_detector.centre.z())
        self.ODD = abs(stage_from_detector.centre.z())
        self.SOD = self.source.centre.distance(self.stage.centre)

        ## Brightest Spot in World Coordinate System:
        self.brightestSpotWorld = copy.deepcopy(self.detector.vectorW)
        self.brightestSpotWorld.scale(self.SDD)
        self.brightestSpotWorld.add(self.source.centre)

        ## Brightest Spot in Detector Coordinate System:
        self.brightestSpotDetector = copy.deepcopy(self.brightestSpotWorld)
        self.brightestSpotDetector.subtract(self.detector.centre)
        
        pxU = self.brightestSpotDetector.dot(self.detector.vectorU) / self.detector.pitchU() + self.detector.cols()/2.0
        pxV = self.brightestSpotDetector.dot(self.detector.vectorV) / self.detector.pitchV() + self.detector.rows()/2.0

        self.brightestSpotDetector = Vector(pxU, pxV, 0)


    def info(self):
        txt  = "Detector\n"
        txt += "===========================================================\n"
        txt += "Center:          {}\n".format(self.detector.centre)
        txt += "u:               {}\n".format(self.detector.vectorU)
        txt += "v:               {}\n".format(self.detector.vectorV)
        txt += "w:               {}\n".format(self.detector.vectorW)
        txt += "Pixels:          {cols} x {rows}\n".format(cols=self.detector.cols(), rows=self.detector.rows())
        txt += "Pitch:           {pitchU} mm x {pitchV} mm\n".format(pitchU=self.detector.pitchU(), pitchV=self.detector.pitchV())
        txt += "Physical Size:   {width} mm x {height} mm\n".format(width=self.detector.physWidth(), height=self.detector.physHeight())
        txt += "Total Tilt:      {tiltRad} rad = {tiltDeg} deg\n".format(tiltRad=self.detectorTotalTilt, tiltDeg=180.0*self.detectorTotalTilt/math.pi)
        txt += "Center Distance: {} mm\n".format(self.detector.centre.distance(self.source.centre))

        # Source - Detector distance (SDD) defined by shortest distance between source and detector:
        txt += "SDD:             {} mm\n".format(self.SDD)
        txt += "Brightest Spot:\n"
        txt += "  World:         {}\n".format(self.brightestSpotWorld)
        txt += "  Pixels:        {}\n".format(self.brightestSpotDetector)

        txt += "\n"
        txt += "Source:\n"
        txt += "===========================================================\n"
        txt += "Center:          {}\n".format(self.source.centre)
        txt += "u:               {}\n".format(self.source.vectorU)
        txt += "v:               {}\n".format(self.source.vectorV)
        txt += "w:               {}\n".format(self.source.vectorW)

        return txt

    def projectionMatrix(self, mode="openCT"):
        validModes = ["openCT", "CERA"]
        if mode in validModes:
            world    = GeometryObject()
            detector = copy.deepcopy(self.detector)
            source   = copy.deepcopy(self.source)
            stage    = copy.deepcopy(self.stage)

            # Scale of the detector CS in units of the world CS: (e.g. mm->px)
            scale_u =  1.0
            scale_v =  1.0
            scale_w = -1.0

            if mode == "CERA":
                # CERA's detector CS has its origin in the lower left corner instead of the centre.
                # Let's move there.

                # Move the detector by half its width and half its height.
                halfWidth  = detector.physWidth()  / 2.0
                halfHeight = detector.physHeight() / 2.0

                detector.centre -= detector.vectorU.scaled(halfWidth)
                detector.centre += detector.vectorV.scaled(halfHeight)

                # The v axis points up instead of down, this also turns the (irrelevant) w normal axis:
                detector.vectorV.scale(-1)
                detector.vectorW.scale(-1)

                # The CERA detector also has a pixel CS instead of a mm CS:
                scale_u = 1.0 / detector.pitchU()
                scale_v = 1.0 / detector.pitchV()
                scale_w = 1.0

            # Save a source CS as seen from the detector CS. This is convenient to
            # later get the SDD, ufoc and vfoc:
            source_from_detector = copy.deepcopy(source)
            source_from_detector.changeReferenceFrame(world, detector)

            # Make the stage CS the new world CS:
            source.changeReferenceFrame(world, stage)
            detector.changeReferenceFrame(world, stage)
            stage.changeReferenceFrame(world, stage)

            # Translation vector from stage to source:
            rfoc = source.centre - stage.centre
            xfoc = rfoc.x()
            yfoc = rfoc.y()
            zfoc = rfoc.z()

            # Focus point on detector: principal, perpendicular ray.
            # In the detector coordinate system, ufoc and vfoc are the u and v coordinates
            # of the source center; SDD (perpendicular to detector plane) is source w coordinate.
            ufoc = source_from_detector.centre.x()
            vfoc = source_from_detector.centre.y()
            wfoc = source_from_detector.centre.z()
            SDD  = abs(wfoc)

            if mode == "CERA":
                # mm -> px
                ufoc = ufoc*scale_u - 0.5
                vfoc = vfoc*scale_v - 0.5

            # Mirror volume:
            M = Matrix(values=[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]])

            # Translation matrix: stage -> source:
            F = Matrix(values=[[1, 0, 0, xfoc], [0, 1, 0, yfoc], [0, 0, 1, zfoc]])

            # Rotations:
            R = basisTransformMatrix(stage, detector)

            # Projection onto detector:
            D = Matrix(values=[[-SDD*scale_u, 0, 0], [0, -SDD*scale_v, 0], [0, 0, scale_w]])

            # Shift in detector CS: (ufoc and vfoc must be in scaled units)
            V = Matrix(values=[[1, 0, ufoc], [0, 1, vfoc], [0, 0, 1]])

            # Multiply all together:
            P = V * (D * (R * (F*M)))

            # Renormalize:
            lower_right = P.get(col=3, row=2)
            if lower_right != 0:
                P.scale(1.0/lower_right)

            return P
        else:
            raise Exception("{} is not a valid mode for the projection matrix computation. Valid modes are: {}".format(mode, validModes))

    def createDetectorFlatField_rays(self):
        """ Calculate an analytical free beam intensity distribution
            picture for the given detector, to be used for an
            ideal flat field correction. """
        width      = self.detector.cols()
        height     = self.detector.rows()
        pixelSizeU = self.detector.pitchU()
        pixelSizeV = self.detector.pitchV()

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
        dx = self.detector.centre.x()
        dy = self.detector.centre.y()
        dz = self.detector.centre.z()

        sx = self.source.centre.x()
        sy = self.source.centre.y()
        sz = self.source.centre.z()

        # Vectors of the detector coordinate system:
        ux = self.detector.vectorU.x()
        uy = self.detector.vectorU.y()
        uz = self.detector.vectorU.z()

        vx = self.detector.vectorV.x()
        vy = self.detector.vectorV.y()
        vz = self.detector.vectorV.z()

        wx = self.detector.vectorW.x()
        wy = self.detector.vectorW.y()
        wz = self.detector.vectorW.z()


        # Angle 'alpha' between detector normal and connection line [detector centre -- source]:
        connectionLine = Vector(dx-sx, dy-sy, dz-sz)

        alpha = abs(self.detector.vectorW.angle(connectionLine))
        if alpha > (math.pi/2):
            alpha = math.pi - alpha

        # Distance between source centre and detector centre:
        dist = self.detector.centre.distance(self.source.centre)

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
                        pixel = self.detector.pixelVectorUpperLeft(x+(gx+1)*stepSize, y+(gy+1)*stepSize)

                        # Grid with no margin:
                        #if gridSize > 1:
                        #    stepSize = 1.0 / (gridSize-1)
                        #    pixel = self.detector.pixelVectorUpperLeft(x+gx*stepSize, y+gy*stepSize)
                        #else:
                        #    pixel = self.detector.pixelVectorCenter(x, y)

                        distToSource = self.source.centre.distance(pixel)

                        # Angle of incident rays:
                        vecSourceToPixel = Vector(pixel.x()-sx, pixel.y()-sy, pixel.z()-sz)
                        incidenceAngle = abs(self.detector.vectorW.angle(vecSourceToPixel))
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
        dx = self.detector.centre.x()
        dy = self.detector.centre.y()
        dz = self.detector.centre.z()

        sx = self.source.centre.x()
        sy = self.source.centre.y()
        sz = self.source.centre.z()

        # Angle 'alpha' between detector normal and connection line [detector centre -- source]:
        connectionLine = Vector(dx-sx, dy-sy, dz-sz)

        alpha = abs(self.detector.vectorW.angle(connectionLine))
        if alpha > (math.pi/2):
            alpha = math.pi - alpha

        # Distance between source centre and detector centre:
        dist = self.detector.centre.distance(self.source.centre)

        # Source - Detector distance (SDD) defined by shortest distance between source and detector,
        # or distance between source and spot of highest intensity on detector.
        SDD = dist * math.cos(alpha)

        # Create a new detector in a coordinate system where source is at (0, 0, 0):
        det = copy.deepcopy(self.detector)
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
        incidenceAngle = abs(self.detector.vectorW.angle(maxCenter))

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
        D = copy.deepcopy(self.detector)
        S = copy.deepcopy(self.source)
        world = GeometryObject()  # will be initialized as world

        S.changeReferenceFrame(world, D)
        D.changeReferenceFrame(world, D)
        D.computeGeometryParameters()

        # Source - Detector distance (SDD) defined by shortest distance between source and detector,
        # or distance between source and spot of highest intensity on detector.
        SDD = abs(S.centre.z())

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
            for p in range(len(coverPolygon.points)):
                coverPolygon.points[p] = coverPolygon.points[p] - S.centre

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
                pA = pA - S.centre
                pB = pB - S.centre
                pC = pC - S.centre
                pD = pD - S.centre

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
        pitchU = self.detector.pitchU()
        pitchV = self.detector.pitchV()

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
        dx = self.detector.centre.x()
        dy = self.detector.centre.y()
        dz = self.detector.centre.z()

        sx = self.source.centre.x()
        sy = self.source.centre.y()
        sz = self.source.centre.z()

        # Vectors of the detector coordinate system:
        ux = self.detector.vectorU.x()
        uy = self.detector.vectorU.y()
        uz = self.detector.vectorU.z()

        vx = self.detector.vectorV.x()
        vy = self.detector.vectorV.y()
        vz = self.detector.vectorV.z()

        wx = self.detector.vectorW.x()
        wy = self.detector.vectorW.y()
        wz = self.detector.vectorW.z()


        # Angle 'alpha' between detector normal and connection line [detector centre -- source]:
        connectionLine = Vector(dx-sx, dy-sy, dz-sz)

        alpha = abs(self.detector.vectorW.angle(connectionLine))
        if alpha > (math.pi/2):
            alpha = math.pi - alpha

        # Distance between source centre and detector centre:
        dist = self.detector.centre.distance(self.source.centre)

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

        upperLeft_u = det.pixelVectorUpperLeft(0, 0).dot(self.detector.vectorU)
        upperLeft_v = det.pixelVectorUpperLeft(0, 0).dot(self.detector.vectorV)
        upperLeft_w = det.pixelVectorUpperLeft(0, 0).dot(self.detector.vectorW)

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

    voxelSizeXY = geo.detector.pitchU() * geo.SOD / geo.SDD
    voxelSizeZ  = geo.detector.pitchV() * geo.SOD / geo.SDD

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
    psu=geo.detector.pitchU(),
    psv=geo.detector.pitchV(),
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
    detectorSizeX=geo.detector.physWidth(),
    detectorSizeY=geo.detector.physHeight(),
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