# -*- coding: UTF-8 -*-
from .test import *
from .general import *
from .geoprimitives import *
from .isrb import Interpolation, OrderError, ResultError
import json
import pkgutil

class Test2D_DW_1_results:
    """ Results for one sub test. """

    def __init__(self):
        self._lineProfileGV  = None   # Profile positions
        self._lineProfilePos = None   # Profile values

        # For the modulation depth vs. wire spacing diagram:
        self._wireSpacings     = None
        self._modulationDepths = None

        self._interpolation_wireSpacings = None
        self._interpolation_modulationDepths = None

        self._interpolation_a = 0
        self._interpolation_b = 0
        self._interpolation_c = 0

        self._pixelSize   = None
        self._SDD         = None
        self._SOD         = None

        self._criticalIndex = None   # Last wire pair with dip > 20%.
        self._SRb_nominal  = None
        self._SRb_interpolated = None


class Test2D_DW_1(generalTest):
    """ CTSimU test 2D-DW-1: Detector Unsharpness. """

    def __init__(self, resultFileDirectory=".", name=None, rawOutput=False):

        generalTest.__init__(
            self,
            testName="2D-DW-1",
            name=name,
            nExpectedRuns=2,
            resultFileDirectory=resultFileDirectory,
            rawOutput=rawOutput)

        # Results for each sub test:
        self._results = []

        # ROI of imaged duplex wire: (622, 312) -- (1630, 678)

        duplexAngle = 3 * (math.pi/180.0)     # 3 deg duplex wire rotation
        duplexCenter = Vector2D(1126, 496)    # Center of duplex wire ROI

        self._profileLength = 1050
        self._profileRes    = 0.1
        self._profileWidth  = 20

        duplexDirection = Vector2D(self._profileLength/2, 0)   # points "right"
        duplexDirection.rotate(duplexAngle)   # Rotate by angle of duplex wire arrangement

        self._p0 = duplexCenter - duplexDirection
        self._p1 = duplexCenter + duplexDirection

        self._currentNominalSRb = 0
        self._currentPixelSize  = 0
        self._currentSDD        = 0
        self._currentSOD        = 0


    def prepare(self):
        """ Preparations before the test will be run with the images from the pipeline. """       
        self._prepared = True

    def prepareRun(self, i):
        if i < len(self._subtests):
            self._jsonScenarioFile = "scenarios/2D-DW-1_Detektor1_2021-05-26v01r01dp.json"
            if self._subtests[i] == "SR75":
                self._jsonScenarioFile = "scenarios/2D-DW-1_Detektor1_2021-05-26v01r01dp.json"
            elif self._subtests[i] == "SR150":
                self._jsonScenarioFile = "scenarios/2D-DW-1_Detektor2_2021-05-26v01r01dp.json"
            #else:
            #    raise Exception("{key} is not a valid subtest identifier for test scenario {test}.".format(key=self._subtests[i], test=self._testName))

            jsonText = pkgutil.get_data(__name__, self._jsonScenarioFile).decode()

            if jsonText != None:
                jsonDict = json.loads(jsonText)

                self._currentNominalSRb = inMM(getFieldOrNone(jsonDict, "detector", "sharpness", "basic_spatial_resolution"))

                if self._currentNominalSRb == None:
                    raise Exception("Test {name}: Cannot find 'basic_spatial_resolution' value in JSON scenario description: {json}".format(name=self._name, json=self._jsonScenarioFile))

                self._currentPixelSize = inMM(getFieldOrNone(jsonDict, "detector", "pixel_pitch", "u"))

                if self._currentPixelSize == None:
                    raise Exception("Test {name}: Cannot find 'pixel_pitch/u' value in JSON scenario description: {json}".format(name=self._name, json=self._jsonScenarioFile))

                self._currentSDD = inMM(getFieldOrNone(jsonDict, "geometry", "detector", "centre", "x"))

                if self._currentSDD == None:
                    raise Exception("Test {name}: Cannot find 'geometry/detector/centre/x' value in JSON scenario description: {json}".format(name=self._name, json=self._jsonScenarioFile))

                self._currentSOD = inMM(getFieldOrNone(jsonDict, "geometry", "stage", "centre", "x"))

                if self._currentSOD == None:
                    raise Exception("Test {name}: Cannot find 'geometry/stage/centre/x' value in JSON scenario description: {json}".format(name=self._name, json=self._jsonScenarioFile))

            else:
                raise Exception("Test {name}: Cannot open JSON scenario description: {json}".format(name=self._name, json=self._jsonScenarioFile))

            # Prepare this test run:
            self.prepare()
            self._results.append(Test2D_DW_1_results())
            self._results[i]._SRb_nominal = self._currentNominalSRb
            self._results[i]._pixelSize = self._currentPixelSize
            self._results[i]._SDD = self._currentSDD
            self._results[i]._SOD = self._currentSOD
           
        else:
            if len(self._subtests) == 0:
                raise Exception("Please provide keywords that identify which metadata file belongs to which subtest. Test {testname} accepts two keywords: 'SR50' for 50um basic spatial resolution and 'SR100' for 100um basic spatial resolution.".format(testname=self._testName))
            else:
                raise Exception("Number of provided image metadata files exceeds number of test runs ({expected}).".format(expected=len(self._subtests)))

    def writeResultFile(self, subname, results):
        pos  = results._lineProfilePos
        prof = results._lineProfileGV

        profileText  = "# Profile data\n"
        profileText += "# s [px]\tGV"
        profileText += "\n"

        for j in range(len(pos)):
            profileText += "{pos:.2f}\t{GV}".format(pos=pos[j], GV=prof[j])
            profileText += "\n"
     
        profileFileName = "{dir}/{name}_{subname}_profile.txt".format(dir=self._resultFileDirectory, name=self._name, subname=subname)
        with open(profileFileName, 'w') as profileFile:
            profileFile.write(profileText)
            profileFile.close()



        resultText  = "# Results\n"
        resultText += "# \n"
        resultText += "# iSRb: {:.5f} mm\n".format(results._SRb_interpolated)
        resultText += "# \n"
        resultText += "# Interpolation function: f(x) = ax²+bx+c\n"
        resultText += "# a: {:.5f}\n".format(results._interpolation_a)
        resultText += "# b: {:.5f}\n".format(results._interpolation_b)
        resultText += "# c: {:.5f}\n".format(results._interpolation_c)
        resultText += "# \n"
        resultText += "# Wire spacing [mm]\tModulation depth [%]\n"

        for j in range(len(results._wireSpacings)):
            resultText += "{ws:.3f}\t{md:.3f}".format(ws=results._wireSpacings[j], md=results._modulationDepths[j])
            resultText += "\n"
     
        resultTextFileName = "{dir}/{name}_{subname}_results.txt".format(dir=self._resultFileDirectory, name=self._name, subname=subname)
        with open(resultTextFileName, 'w') as resultFile:
            resultFile.write(resultText)
            resultFile.close()


    def run(self, image):
        self.prepareRun(self._currentRun)
        i = self._currentRun
        subtestName = self._subtests[i]

        self._results[i]._lineProfileGV, self._results[i]._lineProfilePos, stepsize = image.lineProfile(x0=self._p0._x, y0=self._p0._y, x1=self._p1._x, y1=self._p1._y, width=self._profileWidth, resolution=self._profileRes)

        print("Pixel Size: {}, SOD: {}, SDD: {}".format(self._currentPixelSize, self._currentSOD, self._currentSDD))

        # Use Bendix' iSRb tool to evaluate the basic spatial resolution
        tool = Interpolation(im=image._px, pixelsize=self._currentPixelSize, SOD=self._currentSOD, SDD=self._currentSDD, wire_length=15, wire_spacing=[0.8, 0.63, 0.5, 0.4, 0.32, 0.25, 0.2, 0.16, 0.13, 0.1, 0.08, 0.063, 0.05])

        # The line profile has already been measured. Provide it to the tool:
        tool.measure = self._results[i]._lineProfileGV
        tool.ind     = self._results[i]._lineProfilePos
        tool.calc_dips(prominence=100)
        tool.interpolate()

        if tool.dip20 is not None:
            self._results[i]._SRb_interpolated = tool.dip20
        else:
            self._results[i]._SRb_interpolated = 0
        print("iSRb: {:.5f} mm".format(self._results[i]._SRb_interpolated))

        self._results[i]._criticalIndex = tool.criticalIndex

        self._results[i]._wireSpacings = tool.wire_spacing
        self._results[i]._modulationDepths = tool.dips[0:len(tool.wire_spacing)]

        self._results[i]._interpolation_wireSpacings = tool.inter[0]
        self._results[i]._interpolation_modulationDepths = tool.inter[1]
        self._results[i]._interpolation_a = tool.a
        self._results[i]._interpolation_b = tool.b
        self._results[i]._interpolation_c = tool.c

        self.writeResultFile(subtestName, self._results[i])
        self.plotResults()
        self._currentRun += 1


        return image

    def followUp(self):
        pass
        
    def plotResults(self):
        i = self._currentRun
        subtestName = self._subtests[i]
        try:
            import matplotlib
            import matplotlib.pyplot
            from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

            matplotlib.use("agg")

            fig, (ax1, ax2) = matplotlib.pyplot.subplots(nrows=2, ncols=1, figsize=(9, 7.5))
            
            # Grey Value Profile:
            ax1.plot(self._results[i]._lineProfilePos, self._results[i]._lineProfileGV, linewidth=1.5, label="Measured", color='#1f77b4')

            ax1.set_xlabel("Horizontal distance in px")
            ax1.set_ylabel("Grey Value")
            #ax1.set_xlim([-3*self._results[i]._nominalGaussianSigmaPX, 3*self._results[i]._nominalGaussianSigmaPX])
            ax1.set_title("Grey Value Profile")
            #ax1.xaxis.set_ticklabels([])
            ax1.grid(b=True, which='major', axis='both', color='#d9d9d9', linestyle='dashed')
            ax1.grid(b=True, which='minor', axis='both', color='#e7e7e7', linestyle='dotted')
            ax1.legend(loc='best')


            # Contrast vs. Wire Spacing
            ax2.plot(self._results[i]._interpolation_wireSpacings, self._results[i]._interpolation_modulationDepths, linewidth=1.5, label="Interpolation", color='#ffaa00')

            ax2.plot(self._results[i]._wireSpacings, self._results[i]._modulationDepths, 'o', markersize=4.0, label="Measured modulation depth", color='#1f77b4')
            if self._results[i]._criticalIndex is not None:
                ci = self._results[i]._criticalIndex
                ax2.plot(self._results[i]._wireSpacings[ci], self._results[i]._modulationDepths[ci], 'o', markersize=9.0, label="Critical Pair", markerfacecolor='none', markeredgecolor='black')
            ax2.plot(self._results[i]._SRb_interpolated, 20, 'x', markersize=7.0, label="iSRb", color='#ff0000')

            ax2.set_xlabel("Wire spacing in mm")
            ax2.set_ylabel("Modulation depth in %")
            ax2.set_xlim([0.9, 0])
            ax2.set_title("SRb Interpolation")
            #ax2.xaxis.set_ticklabels([])
            ax2.grid(b=True, which='major', axis='both', color='#d9d9d9', linestyle='dashed')
            ax2.grid(b=True, which='minor', axis='both', color='#e7e7e7', linestyle='dotted')
            ax2.legend(loc='best')

            ax22 = ax2.twiny()
            ax22.set_xlim(ax2.get_xlim())
            ax22.set_xticks(numpy.array([0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0]))
            ax22.set_xticklabels(["0.56", "0.625", "0.714", "0.83", "1", "1.25", "1.67", "2.5", "5", "10", "∞"])
            ax22.set_xlabel("Spatial resolution in lp/mm")


            fig.tight_layout(pad=2.5)

            plotFilename = "{dir}/{name}_{subname}_results.png".format(dir=self._resultFileDirectory, name=self._name, subname=subtestName)
            matplotlib.pyplot.savefig(plotFilename)
            fig.clf()
            matplotlib.pyplot.close('all')
        except:
            log("Warning: Error plotting results for test {name}, {subname} using matplotlib.".format(name=self._name, subname=subtestName))
