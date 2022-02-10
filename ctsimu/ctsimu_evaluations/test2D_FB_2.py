# -*- coding: UTF-8 -*-
from ..test import *

class Test2D_FB_2(generalTest):
    """ CTSimU test 2D-FB-2: 1/r^2 law, cos(alpha) intensity on pixel. """

    def __init__(self, resultFileDirectory=".", name=None, rawOutput=False):

        generalTest.__init__(
            self,
            testName="2D-FB-2",
            name=name,
            nExpectedRuns=1,
            resultFileDirectory=resultFileDirectory,
            rawOutput=rawOutput)
        
        self.geometry = None
        self.analyticalIntensityProfileImage = None  # stores the image after preparation

        # Horizontal profile data:
        self.lineNr = 0
        self.profile_analytical = None
        self.profile_measured   = None

        # Maximum absolute grey value difference between analytical image and mesured image: 
        self.maxDiff = 0

        # RMS of GV differences:
        self.rmsGVDiff = 0


    def prepare(self):
        """ Preparations before the test will be run with the images from the pipeline. """
        if not isinstance(self.pipe, Pipeline):
            self.prepared = False
            raise Exception("Step must be part of a processing pipeline before it can prepare. Current pipeline: {}".format(self.pipe))

        if not self.prepared:
            self.jsonScenarioFile = "ctsimu_evaluations/scenarios/2D-FB-2_2021-03-24v06r00dp.json"

            if(self.jsonScenarioFile != None):
                self.geometry = Geometry(jsonFileFromPkg=self.jsonScenarioFile)
                self.analyticalIntensityProfileImage = self.geometry.createDetectorFlatField_analytical()

                # Raise normalized image to maximum grey value of 60000,
                # as demanded by test 2D-FB-2:
                self.analyticalIntensityProfileImage.renormalize(newMin=0, newMax=60000.0, currentMin=0)

                self.prepared = True
            else:
                raise Exception("Test 2D-FB-2: Please provide a JSON scenario description.")


    def run(self, image):
        self.prepare()
        self.currentRun += 1

        self.lineNr = int(round(self.geometry.brightestSpotDetector.y))
        log("Getting horizontal profile for detector row {l}.".format(l=self.lineNr))

        # Horizontal profile of the analytical image along the central line:
        self.profile_analytical = self.analyticalIntensityProfileImage.horizontalProfile(self.lineNr)

        # Horizontal profile of the image to be evaluated:
        self.profile_measured = image.horizontalProfile(self.lineNr)

        return image


    def followUp(self):
        log("Writing evaluation results...")
        nPixels = len(self.profile_analytical)
        if len(self.profile_analytical) == len(self.profile_measured):
            csvText = "# x [px]\tAnalytical GV\tMeasured GV\tDifference\tRel. Deviation [%]\n"
            squaresSum = 0
            self.maxDiff = 0
            for i in range(nPixels):
                gv0 = self.profile_analytical[i]
                gv1 = self.profile_measured[i]
                delta = gv0 - gv1
                relativeDeviation = 0.0
                if gv0 != 0:
                    relativeDeviation = (100.0*delta)/gv0

                squaresSum += delta*delta

                if abs(delta) > self.maxDiff:
                    self.maxDiff = abs(delta)

                csvText += "{i}\t{gvAnalytical:.3f}\t{gvMeasured:.3f}\t{delta:.3f}\t{relDev:.3f}\n".format(i=i, gvAnalytical=gv0, gvMeasured=gv1, delta=delta, relDev=relativeDeviation)

            self.rmsGVDiff = math.sqrt(squaresSum / nPixels)

            # Write evaluation text file:
            header = "# Evaluation of Test {name}:\n# Max Absolute GV Difference: {maxAbsDiff:.5f}\n# RMS GV Difference:  {rms:.5f}\n# \n# Horizontal profile along y={line}:\n".format(name=self.name, maxAbsDiff=self.maxDiff, rms=self.rmsGVDiff, line=self.lineNr)
            csvText = header + csvText

            csvFileName = "{dir}/{name}_summary.txt".format(dir=self.resultFileDirectory, name=self.name)
            with open(csvFileName, 'w') as csvFile:
                csvFile.write(csvText)
                csvFile.close()


            # Write analytical image:
            if self.rawOutput:
                self.analyticalIntensityProfileImage.saveRAW("{dir}/{name}_analytical.raw".format(dir=self.resultFileDirectory, name=self.name), dataType="float32", addInfo=True)
            else: # TIFF
                self.analyticalIntensityProfileImage.save("{dir}/{name}_analytical.tif".format(dir=self.resultFileDirectory, name=self.name), dataType="float32")

            self.plotResults()

            log("Evaluation data for test {name} written to {dir}.".format(name=self.name, dir=self.resultFileDirectory))

        else:
            raise Exception("Profile width mismatch between analytical result ({}) and measured result ({}).".format(len(profile_analytical), len(profile_measured)))

    def plotResults(self):
        try:
            log("Plotting evaluation results...")
            import matplotlib
            import matplotlib.pyplot
            from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

            xValues = numpy.linspace(0, len(self.profile_analytical), len(self.profile_analytical), endpoint=False)
            absDev = self.profile_analytical - self.profile_measured
            relDev = 100.0*absDev / self.profile_analytical

            matplotlib.use("agg")

            fig, (axUpper, axLower) = matplotlib.pyplot.subplots(nrows=2, ncols=1, figsize=(9, 9))
            
            # Grey Value Profile:
            axUpper.plot(xValues, self.profile_analytical, linewidth=8.0, label="Analytical", color='#ffaa00')
            axUpper.plot(xValues, self.profile_measured,   linewidth=2.0, label="Measured", color='#1f77b4')
            axUpper.set_xlabel("Horizontal distance in px")
            axUpper.set_ylabel("Grey value")
            axUpper.set_title("Grey value line profile along y = {line} px".format(line=self.lineNr))
            axUpper.xaxis.set_major_locator(MultipleLocator(100))
            axUpper.xaxis.set_major_formatter(FormatStrFormatter('%d'))
            axUpper.xaxis.set_minor_locator(MultipleLocator(50))
            axUpper.grid(b=True, which='major', axis='both', color='#d9d9d9', linestyle='dashed')
            axUpper.grid(b=True, which='minor', axis='both', color='#e7e7e7', linestyle='dotted')
            axUpper.legend()

            # Absolute Deviation:
            line_absolute = axLower.plot(xValues, absDev, linewidth=1.25, label="Absolute", color='#ffaa00')
            axLower.set_xlabel("Horizontal distance in px")
            axLower.set_ylabel("Absolute deviation in grey values")
            axLower.set_title("Deviation = Analytical - Measured")
            axLower.xaxis.set_major_locator(MultipleLocator(100))
            axLower.xaxis.set_major_formatter(FormatStrFormatter('%d'))
            axLower.xaxis.set_minor_locator(MultipleLocator(50))
            axLower.grid(b=True, which='major', axis='both', color='#d9d9d9', linestyle='dashed')
            axLower.grid(b=True, which='minor', axis='both', color='#e7e7e7', linestyle='dotted')
            
            # Relative Deviation:
            axLower2 = axLower.twinx()
            axLower2.axhline(0, linestyle='-', linewidth=1.25, color='#000000')
            line_relative = axLower2.plot(xValues, relDev, linewidth=1.25, label="Relative", color='#1f77b4')
            axLower2.set_ylabel("Relative deviation in %")
            maxDev = max(abs(relDev.min()), abs(relDev.max()))
            axLower2.set_ylim([-1.5*maxDev, 1.5*maxDev])

            # Add all relevant lines to legend:
            lines_all = line_absolute + line_relative
            labels = [l.get_label() for l in lines_all]
            axLower.legend(lines_all, labels, loc=0)      

            fig.tight_layout(pad=5.0)

            plotFilename = "{dir}/{name}_results.png".format(dir=self.resultFileDirectory, name=self.name)
            matplotlib.pyplot.savefig(plotFilename)
            fig.clf()
            matplotlib.pyplot.close('all')
        except:
            log("Warning: Error plotting results for test {name} using matplotlib.".format(name=self.name, subname=subtestName))