# -*- coding: UTF-8 -*-
from ..test import *
from ..helpers import *
import json
import pkgutil

class Test2D_FB_1_results:
    """ Results for one sub test of the noise scenario. """

    def __init__(self):
        self.mean = 0
        self.stdDev = 0
        self.snr = 0
        self.fwhm = 0

        # Nominal values from JSON file:
        self.SNRmax = None
        self.nominalMean = None
        self.nominalFWHM = None
        self.nominalSigma = None

        # Grey Value Histogram:
        self.gvHistogram = None  # stores the GV histogram
        self.gvHistNormalized = None  # area-normalized GV histogram


class Test2D_FB_1(generalTest):
    """ CTSimU test 2D-FB-1: Noise. """

    def __init__(self, resultFileDirectory=".", name=None, rawOutput=False):

        generalTest.__init__(
            self,
            testName="2D-FB-1",
            name=name,
            nExpectedRuns=2,
            resultFileDirectory=resultFileDirectory,
            rawOutput=rawOutput)

        self.imax = None    # maximum (noise-free) grey value in free beam, as required by the scenario
        self.gvMax = 70000  # maximum grey value for the histograms
        self.xValues = None # x grey values for the histograms

        # Results for each sub test:
        self.results = []

    def prepare(self):
        """ Preparations before the test will be run with the images from the pipeline. """       
        if not isinstance(self.pipe, Pipeline):
            self.prepared = False
            raise Exception("Step must be part of a processing pipeline before it can prepare. Current pipeline: {}".format(self.pipe))

        if not self.prepared:
            self.results = []
            self.xValues = numpy.linspace(0, self.gvMax, self.gvMax, endpoint=False)

            if self.jsonScenarioFile != None:
                self.prepared = True
            else:
                raise Exception("Test {name}: Please provide a JSON scenario description.".format(name=self.name))

            self.prepared = True

    def prepareRun(self, i):
        if i < len(self.subtests):
            if self.subtests[i] == "SNR100":
                self.jsonScenarioFile = "scenarios/2D-FB-1_Detektor1_SNR100_2021-05-25v06r00dp.json"
            elif self.subtests[i] == "SNR250":
                self.jsonScenarioFile = "scenarios/2D-FB-1_Detektor2_SNR250_2021-05-25v06r00dp.json"
            else:
                raise Exception("{key} is not a valid subtest identifier for test scenario {test}".format(key=self.subtests[i], test=self.testName))

            # Get Imax and the SNRmax from the JSON file:
            jsonText = pkgutil.get_data(__name__, self.jsonScenarioFile).decode()

            if jsonText != None:
                jsonDict = json.loads(jsonText)

                self.imax = getFieldOrNone(jsonDict, "detector", "grey_value", "imax", "value")
                snrMax = getFieldOrNone(jsonDict, "detector", "noise", "snr_at_imax", "value")

                if self.imax is None:
                    raise Exception("Test {name}: Cannot find 'imax' value in JSON scenario description: {json}".format(name=self.name, json=self.jsonScenarioFile))

                if snrMax is None:
                    raise Exception("Test {name}: Cannot find 'snr_at_imax' value in JSON scenario description: {json}".format(name=self.name, json=self.jsonScenarioFile))

                nominalSigma = self.imax / snrMax
                nominalFWHM  = 2.0*math.sqrt(2.0*math.log(2.0))*self.imax / snrMax

            else:
                raise Exception("Test {name}: Cannot open JSON scenario description: {json}".format(name=self.name, json=self.jsonScenarioFile))

            # Prepare this test run:
            self.prepare()
            self.results.append(Test2D_FB_1_results())
            self.results[i].gvHistogram = numpy.zeros(self.gvMax, dtype=numpy.uint32)
            self.results[i].SNRmax = snrMax
            self.results[i].nominalMean = self.imax # after ff-correction
            self.results[i].nominalFWHM = nominalFWHM
            self.results[i].nominalSigma = nominalSigma
           
        else:
            if len(self.subtests) == 0:
                raise Exception("Please provide keywords that identify which metadata file belongs to which subtest. Test {testname} accepts two keywords: 'SNR100' and 'SNR250'.".format(testname=self.testName))
            else:
                raise Exception("Number of provided image metadata files exceeds number of test runs ({expected}).".format(expected=len(self.subtests)))


    def run(self, image):
        if self.currentRun < self.nExpectedRuns:
            self.prepareRun(self.currentRun)
            i = self.currentRun
            subtestName = self.subtests[i]  

            # Create Grey Value Histogram:
            for x in range(image.getWidth()):
                print("\rCalculating grey value histogram... {:.1f}%".format(100*x/image.getWidth()), end="")
                for y in range(image.getHeight()):
                    gv = int(image.getPixel(x, y))
                    if gv >= 0 and gv < self.gvMax:
                        self.results[i].gvHistogram[gv] += 1

            print("\rCalculating grey value histogram... 100%  ")

            self.results[i].mean   = image.mean()
            self.results[i].stdDev = image.stdDev()
            self.results[i].snr    = self.results[i].mean / self.results[i].stdDev
            self.results[i].fwhm   = 2.0*math.sqrt(2.0*math.log(2.0)) * self.results[i].stdDev

            self.results[i].gvHistNormalized = self.results[i].gvHistogram / numpy.sum(self.results[i].gvHistogram)
            
 

            log("Writing evaluation results...")

            # Write evaluation text files:
            result = self.results[i]

            summary  = "# Evaluation of test {name}, {subname}:\n".format(name=self.name, subname=subtestName)
            summary += "#\n"
            summary += "# Signal-to-noise ratio (SNR)\n"
            summary += "# -----------------------------------------\n"
            summary += "# Overall projection SNR:   {snr:.3f}\n".format(snr=result.snr)
            summary += "# Nominal SNR:              {nominalSNR:.3f}\n".format(nominalSNR=result.SNRmax)
            summary += "# Relative SNR deviation:   {reldev:.5f}\n".format(reldev=(result.snr-result.SNRmax)/result.SNRmax)
            summary += "#\n"

            summary += "# Properties of grey value distribution\n"
            summary += "# -----------------------------------------\n"
            summary += "# Mean grey value:          {meanGV:.3f}\n".format(meanGV=result.mean)
            summary += "# Nominal mean grey value:  {nominalMeanGV:.3f}\n".format(nominalMeanGV=result.nominalMean)
            summary += "# Relative mean deviation:  {reldev:.5f}\n".format(reldev=(result.mean-result.nominalMean)/self.imax)
            summary += "#\n"
            summary += "# Measured RMSD:            {stdDev:.3f}\n".format(stdDev=result.stdDev)
            summary += "# Nominal RMSD:             {nominalStdDev:.3f}\n".format(nominalStdDev=result.nominalSigma)
            summary += "# Relative RMSD deviation:  {reldev:.5f}\n".format(reldev=(result.stdDev-result.nominalSigma)/result.nominalSigma)
            summary += "#\n"
            summary += "# Measured noise FWHM:      {fwhm:.3f}\n".format(fwhm=result.fwhm)
            summary += "# Nominal FWHM:             {nominalFWHM:.3f}\n".format(nominalFWHM=result.nominalFWHM)
            summary += "# Relative FWHM deviation:  {reldev:.5f}\n".format(reldev=(result.fwhm-result.nominalFWHM)/result.nominalFWHM)

            resultFileName = "{dir}/{name}_{subname}_summary.txt".format(dir=self.resultFileDirectory, name=self.name, subname=subtestName)
            with open(resultFileName, 'w') as resultFile:
                resultFile.write(summary)
                resultFile.close()


            histogram  = "# Grey value histogram\n"
            histogram += "# -----------------------------------------\n"
            histogram += "# GV\tcounts\tprobability density\n"

            for i in range(len(result.gvHistogram)):
                histogram += "{GV}\t{counts}\t{probDensity}\n".format(GV=i, counts=result.gvHistogram[i], probDensity=result.gvHistNormalized[i])

            histogramFileName = "{dir}/{name}_{subname}_histogram.txt".format(dir=self.resultFileDirectory, name=self.name, subname=subtestName)
            with open(histogramFileName, 'w') as histogramFile:
                histogramFile.write(histogram)
                histogramFile.close()

            # Make graphs:
            try:
                log("Plotting evaluation results...")
                import matplotlib
                import matplotlib.pyplot
                #from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

                # Display only an interval of +-4sigma
                maxNominalSigma = result.nominalSigma

                xStart = self.imax - 4*maxNominalSigma
                xStop  = self.imax + 4*maxNominalSigma

                matplotlib.use("agg")

                fig, ax = matplotlib.pyplot.subplots(nrows=1, ncols=1, figsize=(9, 7))
                
                ax.plot(self.xValues, result.gvHistNormalized, 'o', markersize=1.0, label="Measured (Subtest {subtest})".format(subtest=subtestName), color='#1f77b4')

                ax.set_xlabel("Grey Value")
                ax.set_ylabel("Probability Density")
                ax.set_title("Test 2D-FB-1, {subtest}: Grey Value Distribution".format(subtest=subtestName))
                #ax.xaxis.set_major_locator(MultipleLocator(100))
                #ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
                #ax.xaxis.set_minor_locator(MultipleLocator(50))
                ax.set_xlim(xStart, xStop)
                ax.grid(b=True, which='major', axis='both', color='#d9d9d9', linestyle='dashed')
                ax.grid(b=True, which='minor', axis='both', color='#e7e7e7', linestyle='dotted')
                ax.legend()

                fig.tight_layout(pad=5.0)

                plotFilename = "{dir}/{name}_{subname}_histogram.png".format(dir=self.resultFileDirectory, name=self.name, subname=subtestName)
                matplotlib.pyplot.savefig(plotFilename)
                fig.clf()
                matplotlib.pyplot.close('all')
            except:
                log("Warning: Error plotting results for test {name}, {subname} using matplotlib.".format(name=self.name, subname=subtestName))

            log("Evaluation data for test {name}, {subtest} written to {dir}.".format(name=self.name, subtest=subtestName, dir=self.resultFileDirectory))


            self.currentRun += 1

            # Return image to the pipeline
            return image
        else:
            raise Exception("Number of provided image metadata files exceeds expected number of test runs ({expected}).".format(expected=self.nExpectedRuns))

    def followUp(self):
        pass
        