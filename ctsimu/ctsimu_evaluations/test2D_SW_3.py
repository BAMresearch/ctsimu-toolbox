# -*- coding: UTF-8 -*-
from ..test import *
from ..helpers import *
import math

class Test2D_SW_3_results:
    """ Results for one sub test. """

    def __init__(self):
        # Horizontal Profile Data:
        self.profiles = []  # Horizontal grey value profile for each line


class Test2D_SW_3(generalTest):
    """ CTSimU test 2D-SW-3: Boolean Models. """

    def __init__(self, resultFileDirectory=".", name=None, rawOutput=False):
        generalTest.__init__(
            self,
            testName="2D-SW-3",
            name=name,
            resultFileDirectory=resultFileDirectory,
            rawOutput=rawOutput)

        self.results = []

        self.gv_deviation_limit = 25  # maximum expected in-material grey value deviation

        self.shrink = 0
        self.leftOffset = 450 - self.shrink
        self.nPixels = 100 + 2*self.shrink

        self.steps = [
            ImageROI(self.leftOffset, 818, self.leftOffset+self.nPixels, 910),
            ImageROI(self.leftOffset, 727, self.leftOffset+self.nPixels, 818),
            ImageROI(self.leftOffset, 636, self.leftOffset+self.nPixels, 727),
            ImageROI(self.leftOffset, 545, self.leftOffset+self.nPixels, 636),
            ImageROI(self.leftOffset, 454, self.leftOffset+self.nPixels, 545),
            ImageROI(self.leftOffset, 363, self.leftOffset+self.nPixels, 454),
            ImageROI(self.leftOffset, 272, self.leftOffset+self.nPixels, 363),
            ImageROI(self.leftOffset, 181, self.leftOffset+self.nPixels, 272),
            ImageROI(self.leftOffset,  90, self.leftOffset+self.nPixels, 181),
            ImageROI(self.leftOffset,   0, self.leftOffset+self.nPixels,  90)  # free beam, not treated.
        ]

        # Step thicknesses:
        self.thicknesses = [41.0, 35.89, 30.79, 25.68, 20.58, 15.47, 10.36, 5.26, 0.15, 0.0]

        # Those pixel lines are between the steps of the wedge:
        self.forbidden_lines = [90, 181, 272, 363, 454, 545, 636, 727, 818, 909]

        # Use Beer-Lambert to calculate expected attenuation at 80 keV:
        # Data source: NIST https://physics.nist.gov/PhysRefData/XrayMassCoef/tab3.html
        muAl = 0.0544638 # /mm
        muTi = 0.18234   # /mm
        self.attenuationsAl = []
        self.attenuationsTi = []

        for d in self.thicknesses:
            attAl = math.exp(-d*muAl)
            attTi = math.exp(-d*muTi)

            self.attenuationsAl.append(attAl)
            self.attenuationsTi.append(attTi)


    def prepare(self):
        """ Preparations before the test will be run with the images from the pipeline. """
        self.prepared = True

    def prepareRun(self, i):
        if i < len(self.subtests):
            if self.subtests[i] == "Al_Al":
                pass
            elif self.subtests[i] == "Al_Ti":
                pass
            else:
                raise Exception("{key} is not a valid subtest identifier for test scenario {test}".format(key=self.subtests[i], test=self.testName))

            self.results.append(Test2D_SW_3_results())
        else:
            if len(self.subtests) == 0:
                raise Exception("Please provide keywords that identify which metadata file belongs to which subtest. Test {testname} accepts the following keywords: 'Al_Al' and 'Al_Ti'.".format(testname=self.testName))
            else:
                raise Exception("Number of provided image metadata files exceeds number of test runs ({expected}).".format(expected=len(self.subtests)))

    def is_near(self, expected, other):
        if expected != 0:  # relative deviation no more than 5%
            return ((abs(other-expected)/abs(expected)) < 0.05)   # 5% maximum deviation from expectation value
        else:
            return (abs(other) < self.gv_deviation_limit)   # If the expected difference is 0, the maximum abs. GV deviation should be within the expected limit.

    def run(self, image):
        self.prepare()
        self.prepareRun(self.currentRun)
        i = self.currentRun
        subtestName = self.subtests[i]

        detectorSize = 1000

        air = self.steps[-1]  # The last ROI is the air ROI.
        airStats = image.stats(air)
        gvAir = airStats["mean"]

        print("Air grey value: {}".format(gvAir))

        # Calculate a derivative image:
        derivative = copy.deepcopy(image)
        derivative.px[:,:-1] = derivative.px[:,1:] - derivative.px[:,:-1]
        derivative.crop(0, 0, 999, 1000)

        # Clean using the 0.04% deviation threshold.
        self.gv_deviation_limit = gvAir * 4e-4  # (0.04% of 60000 = 24)
        print("GV deviation limit: {}".format(self.gv_deviation_limit))

        derivative.px[numpy.where(numpy.absolute(derivative.px) < self.gv_deviation_limit)] = 0

        derivative.saveTIFF(filename="{dir}/{name}_{subname}_derivative_cleaned.tif".format(dir=self.resultFileDirectory, name=self.name, subname=subtestName), dataType="float32")


        # Calculate expected grey values:
        gvAl = [gvAir*att for att in self.attenuationsAl]
        gvTi = [gvAir*att for att in self.attenuationsTi]
        gvRight = gvAl  # grey values for the wedge material to the right (either Al or Ti)

        if subtestName == "Al_Ti":
            gvRight = gvTi

        for step in range(0, 9):  # only treat the first nine (true) steps, ignore free beam
            roi = self.steps[step]
            y0 = roi.y0
            y1 = roi.y1

            for y in range(y1-1, y0-1, -1):
                # -- LEFT BOUNDARY
                """
                # Calculate GV drop for the left boundary at x=89,90:
                sum_across_left_boundary = numpy.sum(derivative.px[y,89:91])

                # The sum over these two columns should represent a drop from GV(air) to GV(Al)
                if self.is_near(-diff_air_al, sum_across_left_boundary):
                    # Set both of them to zero so that they are not counted as an anomaly:
                    derivative.px[y,89:91] = 0
                elif y in self.forbidden_lines:
                    # If this is a forbidden line, only check if all slopes point downward.
                    if not any(diff > 0 for diff in derivative.px[y,89:91]):
                        derivative.px[y,89:91] = 0
                    else:
                        # Set one of the pixels to the expectation value, the others to zero. -> introduces anomaly
                        derivative.px[y,89] = 0
                        derivative.px[y,90] = -diff_air_al
                else:
                    # Set one of the pixels to the expectation value, the others to zero. -> introduces anomaly
                    derivative.px[y,89] = 0
                    derivative.px[y,90] = -diff_air_al
                

                # -- RIGHT BOUNDARY
                x_bevel = 883 - y
                x_right = 908
                padding = 2
                if x_bevel < 0:
                    x_right += x_bevel
                    padding = 5

                x0 = x_right - padding
                x1 = x_right + padding
                sum_across_right_boundary = numpy.sum(derivative.px[y,x0:x1+1])
                if self.is_near(diff_air_right, sum_across_right_boundary):
                    derivative.px[y,x0:x1+1] = 0
                elif y in self.forbidden_lines:
                    # If this is a forbidden line, only check if all slopes point upward.
                    if not any(diff < 0 for diff in derivative.px[y,x0:x1+1]):
                        derivative.px[y,x0:x1+1] = 0
                    else:
                        # Set one of the pixels to the expectation value, the others to zero. -> introduces anomaly
                        derivative.px[y,x0:x1+1] = 0
                        derivative.px[y,x0] = diff_air_right
                else:
                    # Set one of the pixels to the expectation value, the others to zero. -> introduces anomaly
                    derivative.px[y,x0:x1+1] = 0
                    derivative.px[y,x0] = diff_air_right
                """

                # We do not treat the outer boundaries anymore.
                # They are not subject to Boolean models or contiguous surfaces.
                # Instead, set the derivatives in these regions to zero.
                derivative.px[y,80:100] = 0
                derivative.px[y,870:920] = 0



                # Width of air gap in the upper third of the wedge, for upper and lower y-line boundary:
                w_gap_upper = (-11*5/(3*1000))*(y - (4*1000)/11)
                w_gap_lower = (-11*5/(3*1000))*((y+1) - (4*1000)/11)

                # Position of material interface:
                x_mat_upper = (-11*10/(3*1000))*(y - (7*1000)/11)
                x_mat_lower = (-11*10/(3*1000))*((y+1) - (7*1000)/11)

                diff_air_al = gvAir - gvAl[step]
                diff_air_right = gvAir - gvRight[step]
                diff_right_al  = gvRight[step] - gvAl[step]
                

                # -- CENTRAL BOUNDARIES
                if w_gap_lower > 0:  # air gap, steps 7, 8, 9
                    # Right x locations of the gap boundary (left is always at 500).
                    # In the variable names, left and right refer to the leftmost
                    # and rightmost pixel to consider for the sloped edge.
                    x_gap_right0 = 500 + int(math.trunc(w_gap_lower)) - 1
                    x_gap_right1 = 500 + int(math.trunc(w_gap_upper))

                    # If the gap is far enough to separate Al, air and Al/Ti:
                    if w_gap_lower > 1:
                        # Check if there is a GV rise to air at x=498,499
                        sum_across_boundary_to_air = numpy.sum(derivative.px[y,498:500])
                        if self.is_near(diff_air_al, sum_across_boundary_to_air):
                            # Set both of them to zero so that they are not counted as an anomaly:
                            derivative.px[y,498:500] = 0
                        elif (y in self.forbidden_lines) and (not any(diff < 0 for diff in derivative.px[y,498:500])):
                            derivative.px[y,498:500] = 0
                        else:
                            # Set one of the pixels to the expectation value, the others to zero. -> introduces anomaly
                            derivative.px[y,498:500] = 0
                            derivative.px[y,499] = diff_air_al

                        # Check if there is a GV drop to the material at the right of the gap (Al/Ti):
                        sum_across_boundary_to_right_material = numpy.sum(derivative.px[y,x_gap_right0:x_gap_right1+1])
                        if self.is_near(-diff_air_right, sum_across_boundary_to_right_material):
                            # Set right gap boundary pixels to zero so that they are not counted as an anomaly:
                            derivative.px[y,x_gap_right0:x_gap_right1+1] = 0
                        elif (y in self.forbidden_lines) and (not any(diff > 0 for diff in derivative.px[y,x_gap_right0:x_gap_right1+1])):
                            derivative.px[y,x_gap_right0:x_gap_right1+1] = 0
                        else:
                            # Set one of the pixels to the expectation value, the others to zero. -> introduces anomaly
                            derivative.px[y,x_gap_right0:x_gap_right1+1] = 0
                            derivative.px[y,x_gap_right0] = -diff_air_right
                    else:  # if the gap is too narrow to separate left and right materials from the air gap
                        # Only the total sum over the gap has to match the material GV difference:
                        sum_across_complete_gap = numpy.sum(derivative.px[y,498:x_gap_right1+1])
                        if self.is_near(diff_right_al, sum_across_complete_gap) or (y in self.forbidden_lines):
                            derivative.px[y,498:x_gap_right1+1] = 0
                        else:
                            # Set one of the pixels to the expectation value, the others to zero. -> introduces anomaly
                            derivative.px[y,498:x_gap_right1+1] = 0
                            derivative.px[y,499] = diff_right_al
                elif subtestName == "Al_Ti":
                    # We can skip Al_Al because the expected boundary difference is zero.

                    x_boundary_left  = 499 - 1
                    x_boundary_right = 499 + 1

                    if x_mat_lower < 0:  # lower wedge in steps 1, 2, 3
                        x_boundary_left  += int(math.trunc(x_mat_lower))
                        x_boundary_right += int(math.trunc(x_mat_upper))

                    # Check if sum over central boundary matches expected grey value drop from left to right material:
                    sum_across_boundary = numpy.sum(derivative.px[y,x_boundary_left:x_boundary_right+1])
                    if self.is_near(diff_right_al, sum_across_boundary):
                        derivative.px[y,x_boundary_left:x_boundary_right+1] = 0
                    elif (y in self.forbidden_lines) and (not any(diff > 0 for diff in derivative.px[y,x_boundary_left:x_boundary_right+1])):
                        derivative.px[y,x_boundary_left:x_boundary_right+1] = 0
                    else:
                        # Set one of the pixels to the expectation value, the others to zero. -> introduces anomaly
                        derivative.px[y,x_boundary_left:x_boundary_right+1] = 0
                        derivative.px[y,x_boundary_right] = diff_right_al



        derivative.saveTIFF(filename="{dir}/{name}_{subname}_anomalies.tif".format(dir=self.resultFileDirectory, name=self.name, subname=subtestName), dataType="float32")

        absDevSum = numpy.sum(numpy.absolute(derivative.px))
        anomalies = numpy.where(numpy.absolute(derivative.px) > 0)
        anomalyCount = len(anomalies[0])
        anomalyMean = "-"
        if anomalyCount > 0:
            anomalyMean = absDevSum / anomalyCount

        print("Number of detected grey value anomalies: {}".format(anomalyCount))
        print("Mean grey value difference of all anomalies: {}".format(anomalyMean))

        summaryText  = "# Evaluation of Test {name}, {subname}:\n".format(name=self.name, subname=subtestName)
        summaryText += "# \n"
        summaryText += "# Number of detected grey value anomalies: {} \n".format(anomalyCount)
        summaryText += "# Mean anomaly grey value: {} \n".format(anomalyMean)

        if anomalyCount > 0:
            summaryText += "# \n"
            if anomalyCount <= 60000:
                summaryText += "# Anomaly positions:\n"
                summaryText += "# x\ty\n"
                for i in range(len(anomalies[0])):
                    summaryText += "{x}\t{y}\n".format(x=anomalies[1][i], y=anomalies[0][i])
            else:
                summaryText += "# Too many anomalies for a text file. Please use the anomaly image file to inspect anomalies.\n"

        summaryFileName = "{dir}/{name}_{subname}_summary.txt".format(dir=self.resultFileDirectory, name=self.name, subname=subtestName)
        with open(summaryFileName, 'w') as summaryFile:
            summaryFile.write(summaryText)
            summaryFile.close()

        #self.plotResults()

        self.currentRun += 1
        return image

    def followUp(self):
        pass

    def plotResults(self):
        pass