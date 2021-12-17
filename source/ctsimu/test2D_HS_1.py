# -*- coding: UTF-8 -*-
from .test import *
from .geometry import Vector, inMM
from .general import *
import json
import pkgutil

class ballCollector():
    """ We will add the circles to the ball collector. It decides which circle is which. """

    def __init__(self):
        # Prepare a list to contain the circle vectors:
        self._n_circles_expected = 10
        self._circlesUnordered = []
        self._circlesInOrder = [None] * self._n_circles_expected
        #self._distancesToCorners = [None] * self._n_circles_expected

        self._greatCircleIndex = 0
        self._maxCircleRadius = 0


    def addCircle(self, coords, R):  # tuple: coords = (cx, cy)
        if self._maxCircleRadius < R:  # Remember the biggest circle
            self._maxCircleRadius = R
            self._greatCircleIndex = len(self._circlesUnordered)

        circle2Dpos = Vector(coords[0], coords[1], 0)
        self._circlesUnordered.append(circle2Dpos)

    def sortCircles(self):
        if len(self._circlesUnordered) == self._n_circles_expected:
            # Find all connection vectors from the great hole to the other holes:
            greatCircle = self._circlesUnordered.pop(self._greatCircleIndex)
            self._circlesInOrder[0] = greatCircle

            # circlesUnordered now only contains the 9 small circles, and must remain so (not to mess up indices).

            ## Create list of connection vectors:
            connectionVectors = []
            for i in range(len(self._circlesUnordered)):
                connectionVectors.append(Vector.connection(vecFrom=greatCircle, vecTo=self._circlesUnordered[i]))

            ## Find vector pair with smallest and biggest angle:
            small = None
            hole9idx = None
            
            smallestAngle = math.pi
            smallestDistance = connectionVectors[0].length()

            for i in range(len(connectionVectors)):
                # Hole 9 is the hole with the smallest distance to the great hole:
                if smallestDistance > connectionVectors[i].length():
                    hole9idx = i
                    smallestDistance = connectionVectors[i].length()

                # Find angles between connection vectors of each circle pair:
                for j in range(i+1, len(connectionVectors)):
                    angle = connectionVectors[i].angle(connectionVectors[j])
                    if smallestAngle > abs(angle):
                        smallestAngle = abs(angle)
                        small = (i, j)
                    

            # Hole 1 is the hole from the connection vector pair with the smallest angle
            # that is also the farthest away from the great circle.
            # Hole 5 is the one that is closer to the great circle.
            connection0to1 = None

            if self._circlesUnordered[small[0]].distance(greatCircle) > self._circlesUnordered[small[1]].distance(greatCircle):
                self._circlesInOrder[1] = self._circlesUnordered[small[0]]
                self._circlesInOrder[5] = self._circlesUnordered[small[1]]
            else:
                self._circlesInOrder[5] = self._circlesUnordered[small[0]]
                self._circlesInOrder[1] = self._circlesUnordered[small[1]]

            # Hole 9 is one of the circles with the smallest distance to the big hole:
            self._circlesInOrder[9] = self._circlesUnordered[hole9idx]

            # Diagonal from hole 9 to hole 1:
            diagonal91 = Vector.connection(self._circlesInOrder[9], self._circlesInOrder[1])

            remainingCirclesAndDiagonalAngles = []
            for i in range(len(connectionVectors)):
                if not ((i in small) or (i == hole9idx)):
                    angle = diagonal91.angle(connectionVectors[i])
                    remainingCirclesAndDiagonalAngles.append((i, angle))

            # Sort list by angles:
            remainingCirclesAndDiagonalAngles = sorted(remainingCirclesAndDiagonalAngles, key=lambda x: x[1])
            i = 0
            for circle in remainingCirclesAndDiagonalAngles:
                # The sign of the cross product's z component tells us
                # which side of the diagonal the circle lies on:
                idx = circle[0]
                angle = circle[1]

                cross = diagonal91.cross(connectionVectors[idx])

                if i < 2: # First vector pair with smallest angle towards diagonal
                    if cross.z() > 0:  # To the "right"
                        self._circlesInOrder[2] = self._circlesUnordered[idx]
                    else:  # To the "left"
                        self._circlesInOrder[4] = self._circlesUnordered[idx]
                elif i < 4: # Second vector pair with medium-sized angle towards diagonal
                    if cross.z() > 0:  # To the "right"
                        self._circlesInOrder[3] = self._circlesUnordered[idx]
                    else:  # To the "left"
                        self._circlesInOrder[7] = self._circlesUnordered[idx]
                elif i < 6: # Second vector pair with medium-sized angle towards diagonal
                    if cross.z() > 0:  # To the "right"
                        self._circlesInOrder[6] = self._circlesUnordered[idx]
                    else:  # To the "left"
                        self._circlesInOrder[8] = self._circlesUnordered[idx]

                i += 1
        else:
            raise Exception("A total of {} holes is needed. Found: {} holes.".format(self._n_circles_expected, len(self._circlesUnordered)))
    

    def nCircles(self):
        return len(self._circlesInOrder)

    def nCirclesFound(self):
        n = 0
        for c in self._circlesInOrder:
            if c != None:
                n = n+1

        return n

    def getCircle(self, i):
        if i >= 0 and i < self.nCircles():
            return self._circlesInOrder[i]
        else:
            raise Exception("Circle with index {i} does not exist.".format(i=i))

    def scale(self, factor, aroundHole):
        for i in range(len(self._circlesInOrder)):
            if i != aroundHole:
                connection = Vector.connection(self._circlesInOrder[aroundHole], self._circlesInOrder[i])
                connection.scale(factor)
                connection.add(self._circlesInOrder[aroundHole])
                self._circlesInOrder[i] = connection

    def rotate(self, angle, aroundHole):
        axis = Vector(0, 0, -1)
        for i in range(len(self._circlesInOrder)):
            if i != aroundHole:
                connection = Vector.connection(self._circlesInOrder[aroundHole], self._circlesInOrder[i])
                connection.rotate(axis, angle)
                connection.add(self._circlesInOrder[aroundHole])
                self._circlesInOrder[i] = connection

    def distance(self, a, b):
        if (a != None) and (b != None):
            if (a >= 0) and (a < self.nCircles()):
                if (b >= 0) and (b < self.nCircles()):
                    if self._circlesInOrder[a] != None:
                        if self._circlesInOrder[b] != None:
                            return self._circlesInOrder[a].distance(self._circlesInOrder[b])
                        else:
                            raise Exception("Error: Circle {} not found. Cannot continue with test.".format(b))
                    else:
                        raise Exception("Error: Circle {} not found. Cannot continue with test.".format(a))

        raise Exception("Requested Distance: Index for circle a ({a}) or b ({b}) out of bounds.".format(a=a, b=b))

    def printResults(self):
        text = ""
        for i in range(self.nCircles()):
            text += "Hole {:02d}: ".format(i)
            if self._circlesInOrder[i] == None:
                text += "Not found!"
            else:
                text += "   cx = {:.3f} px,\tcy = {:.3f} px".format(self._circlesInOrder[i].x(), self._circlesInOrder[i].y())

            text += "\n"

        return text

    def printResults_tabular(self):
        text = "Hole\tcenter_x [px]\tcenter_y [px]\n"
        for i in range(self.nCircles()):
            text += "{:02d}\t".format(i)
            if self._circlesInOrder[i] == None:
                text += "Not found!"
            else:
                text += "{:.3f}\t{:.3f}".format(self._circlesInOrder[i].x(), self._circlesInOrder[i].y())

            text += "\n"

        return text


class Test2D_HS_1(generalTest):
    """ CTSimU test 2D-HS-1: object positioning. """

    def __init__(self, resultFileDirectory=".", name=None, rawOutput=False):
        
        generalTest.__init__(
            self,
            testName="2D-HS-1",
            name=name,
            nExpectedRuns=1,
            resultFileDirectory=resultFileDirectory,
            rawOutput=rawOutput)

        self._success = False

        # Future circle collector (class ballCollector):
        self._circleCollection = None

        self._nominalCircles = ballCollector()

        # Detector parameters:
        self._pixelSizeX = None
        self._pixelSizeY = None
        self._detectorWidthPx  = None
        self._detectorHeightPx = None
        self._detectorWidthMM  = None
        self._detectorHeightMM = None

        # Results of the evaluation:
        self._meanTranslationVector = None
        self._inPlaneRotation = None
        self._scaleFactor = None

    def worldToDetector(self, yWorld, zWorld):
        # mirror x and y:
        x = -yWorld
        y = -zWorld

        x += self._detectorWidthMM/2
        y += self._detectorHeightMM/2

        x /= self._pixelSizeX
        y /= self._pixelSizeY

        return (x, y)

    def prepare(self):
        """ Preparations before the test will be run with the images from the pipeline. """
        self._jsonScenarioFile = "scenarios/2D-HS-1_2021-03-24v02r00dp.json"
        self._prepared = True

        # Get Imax and the SNRmax from the JSON file:
        jsonText = pkgutil.get_data(__name__, self._jsonScenarioFile).decode()

        if jsonText != None:
            jsonDict = json.loads(jsonText)

            self._pixelSizeX = inMM(getFieldOrNone(jsonDict, "detector", "pixel_pitch", "u"))
            self._pixelSizeY = inMM(getFieldOrNone(jsonDict, "detector", "pixel_pitch", "v"))
            self._detectorWidthPx  = getFieldOrNone(jsonDict, "detector", "columns", "value")
            self._detectorHeightPx = getFieldOrNone(jsonDict, "detector", "rows", "value")

            self._detectorWidthMM  = self._detectorWidthPx  * self._pixelSizeX
            self._detectorHeightMM = self._detectorHeightPx * self._pixelSizeY

            self._nominalCircles.addCircle(self.worldToDetector(yWorld=-160, zWorld=-140), R=133)
            self._nominalCircles.addCircle(self.worldToDetector(yWorld= 205, zWorld= 205), R=66)
            self._nominalCircles.addCircle(self.worldToDetector(yWorld=   0, zWorld= 205), R=66)
            self._nominalCircles.addCircle(self.worldToDetector(yWorld=-205, zWorld= 205), R=66)
            self._nominalCircles.addCircle(self.worldToDetector(yWorld= 205, zWorld=   0), R=66)
            self._nominalCircles.addCircle(self.worldToDetector(yWorld=   0, zWorld=   0), R=66)
            self._nominalCircles.addCircle(self.worldToDetector(yWorld=-205, zWorld=   0), R=66)
            self._nominalCircles.addCircle(self.worldToDetector(yWorld= 205, zWorld=-205), R=66)
            self._nominalCircles.addCircle(self.worldToDetector(yWorld=   0, zWorld=-205), R=66)
            self._nominalCircles.addCircle(self.worldToDetector(yWorld=-205, zWorld=-205), R=66)
            
            self._nominalCircles.sortCircles()

    def run(self, image):
        self.prepare()
        self._currentRun += 1

        self._success = False

        log("Running test {}".format(self._testName))
        n_circles_expected = 10

        holeImg = copy.deepcopy(image)

        # Thresholding splits image into black and white, to remove artifacts:
        thresh = holeImg.max() - 0.1*(holeImg.max()-holeImg.min())
        log("Binarization threshold: {:.3f}".format(thresh))
        holeImg.applyThreshold(threshold=thresh, lower=0, upper=60000)

        if self._rawOutput:
            holeImg.saveRAW("{dir}/{name}_01_binarized.raw".format(dir=self._resultFileDirectory, name=self._name), dataType="uint16", addInfo=True)
        else: # TIFF
            holeImg.save("{dir}/{name}_01_binarized.tif".format(dir=self._resultFileDirectory, name=self._name), dataType="uint16")

        # Clean big patch(es) that are not circles:
        nPatches, nCleaned, nRemaining, patchGeometry = holeImg.cleanPatches(min_patch_area=(5*5), max_patch_area=None, remove_border_patches=True)

        if self._rawOutput:
            holeImg.saveRAW("{dir}/{name}_02_cleaned.raw".format(dir=self._resultFileDirectory, name=self._name), dataType="uint16", addInfo=True)
        else: # TIFF
            holeImg.save("{dir}/{name}_02_cleaned.tif".format(dir=self._resultFileDirectory, name=self._name), dataType="uint16")

        log("Looking for {} circles, found {} structures.".format(n_circles_expected, nRemaining))

        if(nRemaining != n_circles_expected):
            log("Error: Cannot continue with test.")
        else:
            log("That's good. Trying to fit circles to the following structures at center coordinates (cx, cy).")

            # Find edges:
            holeImg.filter_edges(mode="sobel")
            thresh = (holeImg.max()+holeImg.min())/2.0
            log("Binarization threshold: {:.3f}".format(thresh))
            holeImg.applyThreshold(threshold=thresh, lower=0, upper=60000)
            if self._rawOutput:
                holeImg.saveRAW("{dir}/{name}_03_edges.raw".format(dir=self._resultFileDirectory, name=self._name), dataType="uint16", addInfo=True)
            else:
                holeImg.save("{dir}/{name}_03_edges.tif".format(dir=self._resultFileDirectory, name=self._name), dataType="uint16")

            i=0
            greatRadius = 0
            greatCircle = -1
            circles = []

            log("Running the circle fits...")

            for cx, cy, width, height in patchGeometry:
                #log("  Structure {:02d}: cx={}, cy={}, width={}, height={}".format(i, cx, cy, width, height))
                circleImg = copy.deepcopy(holeImg)

                x = int(cx)
                y = int(cy)
                w = int(1.8*width)
                h = int(1.8*height)

                # Crop circle:
                left, right, top, bottom = circleImg.cropROIaroundPoint(x, y, w, h)
                if self._rawOutput:
                    circleImg.saveRAW("{dir}/{name}_04_circle_{n:02d}.raw".format(dir=self._resultFileDirectory, name=self._name, n=i), dataType="uint16", addInfo=True)
                else:
                    circleImg.save("{dir}/{name}_04_circle_{n:02d}.tif".format(dir=self._resultFileDirectory, name=self._name, n=i), dataType="uint16")

                # Fit circle:
                #log("    Fitting circle...")
                centerX, centerY, radius, meanDifference, minDifference, maxDifference = circleImg.fitCircle()
                centerX += left
                centerY += top

                log("    Found: cx={:.4f}, cy={:.4f}, R={:.4f}".format(centerX, centerY, radius))

                # Find the biggest circle. This will be the symmetry breaker.
                if radius > greatRadius:
                    greatRadius = radius
                    greatCircle = i

                circles.append((centerX, centerY, radius))

                i += 1

            if greatCircle >= 0:
                greatX, greatY, greatR = circles[greatCircle]
                self._circleCollection = ballCollector()

                for i in range(len(circles)):
                    cx, cy, R = circles[i]
                    self._circleCollection.addCircle((cx, cy), R)

                self._success = True

            else:
                log("No circles found. Cannot continue with test.")

        return image

    def followUp(self):
        if self._success:
            log("Writing evaluation results...")
            self._circleCollection.sortCircles()

            summaryText = "# Evaluation of test {}\n".format(self._testName)

            summaryText += "#\n"
            summaryText += "# Hole Coordinates:\n"

            log("Hole locations found:")
            log(self._circleCollection.printResults())

            summaryText += self._circleCollection.printResults_tabular()

            if self._circleCollection.nCirclesFound() < self._circleCollection._n_circles_expected:
                log("Only {n} of the expected {N} holes were correctly identified. Test cannot continue.".format(n=self._circleCollection.nCirclesFound(), N=self._circleCollection._n_circles_expected))
                return
            else:
                log("Okay, all circles identified.")

            # Mean Scale:
            scaleAgainstCircle = 5
            summaryText += "\n"
            summaryText += "Scale factor,\n"
            summaryText += "from hole distances: scale = 1 - dist(ideal)/dist(real)\n"
            summaryText += "===============================================================\n"
            scaleDevList = []
            for i in range(self._nominalCircles.nCircles()-1):
                for j in range(i+1, self._nominalCircles.nCircles()):
                    length_ideal = self._nominalCircles.getCircle(i).distance(self._nominalCircles.getCircle(j))
                    length_real = self._circleCollection.getCircle(i).distance(self._circleCollection.getCircle(j))
                    scale = 1.0 - length_ideal / length_real

                    summaryText += "Scale dev. ({i}, {j}): {space}{scale:.7f}\n".format(i=i, j=j, space=" "*(scale>=0), scale=scale)
                    
                    scaleDevList.append(scale)

            meanScaleDeviation, stdDevScaleDeviation = listMeanAndStdDev(scaleDevList)

            summaryText += "---------------------------------------------------------------\n"
            summaryText += "Mean scale deviation:   {space}{value:.7f}\n".format(space=" "*(meanScaleDeviation>=0), value=meanScaleDeviation)
            summaryText += "StdDev scale deviation: {space}{value:.7f}\n".format(space=" "*(stdDevScaleDeviation>=0), value=stdDevScaleDeviation)

            ## Scale the real circle collection back by that scale factor...
            self._circleCollection.scale(factor=1-meanScaleDeviation, aroundHole=scaleAgainstCircle)


            # Mean rotation angle:
            rotateAroundCircle = 5
            summaryText += "\n"
            summaryText += "In-plane rotation\n"
            summaryText += "===============================================================\n"
            meanAngle = 0
            stdDevAngle = 0
            angleList = []

            for i in range(self._nominalCircles.nCircles()-1):
                for j in range(i+1, self._nominalCircles.nCircles()):
                    connection_ideal = Vector.connection(self._nominalCircles.getCircle(i), self._nominalCircles.getCircle(j))
                    connection_real  = Vector.connection(self._circleCollection.getCircle(i), self._circleCollection.getCircle(j))
                    angle = connection_ideal.angle(connection_real)
                    
                    # angular orientation:
                    cross = connection_ideal.cross(connection_real)
                    if (cross.z() < 0) and (abs(angle) < (math.pi - 0.0000001)):
                        angle = -abs(angle)
                    else:
                        angle = abs(angle)

                    summaryText += "Tilt angle ({i}, {j}): {space}{angle:.7f} rad = {space}{angleDeg:.7f} deg\n".format(i=i, j=j, space=" "*(angle>=0), angle=angle, angleDeg=(angle*180.0/math.pi))

                    angleList.append(angle)

            meanAngle, stdDevAngle = listMeanAndStdDev(angleList)

            summaryText += "---------------------------------------------------------------\n"
            summaryText += "Mean rotation angle:   {space}{value:.7f} rad = {space}{angleDeg:.7f} deg\n".format(space=" "*(meanAngle>=0), value=meanAngle, angleDeg=(meanAngle*180.0/math.pi))
            summaryText += "StdDev rotation angle: {space}{value:.7f} rad = {space}{angleDeg:.7f} deg\n".format(space=" "*(stdDevAngle>=0), value=stdDevAngle, angleDeg=(stdDevAngle*180.0/math.pi))
            summaryText += "\n"

            ## Rotate the real circle collection back by that angle...
            self._circleCollection.rotate(meanAngle, rotateAroundCircle)


            # Mean translation vector:
            summaryText += "Translation analysis\n"
            summaryText += "===============================================================\n"
            translationsX = []
            translationsY = []

            for i in range(self._nominalCircles.nCircles()):
                translation = Vector.connection(self._nominalCircles.getCircle(i), self._circleCollection.getCircle(i))
                summaryText += "Translation hole {i:02d} [px]: {vec}\n".format(i=i, vec=translation)
                
                translationsX.append(translation.x())
                translationsY.append(translation.y())

            meanTranslationX, stdDevTranslationX = listMeanAndStdDev(translationsX)
            meanTranslationY, stdDevTranslationY = listMeanAndStdDev(translationsY)

            meanTranslation = Vector(meanTranslationX, meanTranslationY, 0)
            stdDevTranslation = Vector(stdDevTranslationX, stdDevTranslationY, 0)
            
            summaryText += "---------------------------------------------------------------\n"
            summaryText += "Mean translation vector [px]: {}\n".format(meanTranslation)
            summaryText += "StdDev translation [px]:      {}".format(stdDevTranslation)


            # Special Summary Line:
            #summaryText += "\n\n"
            #summaryText += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(meanScaleDeviation, stdDevScaleDeviation, meanAngle, stdDevAngle, meanTranslation.x(), stdDevTranslation.x(), meanTranslation.y(), stdDevTranslation.y())


            resultFileName = "{}/{}_summary.txt".format(self._resultFileDirectory, self._name)
            with open(resultFileName, 'w') as resultFile:
                resultFile.write(summaryText)
                resultFile.close()

            self.plotResults()


    def plotResults(self):
        return