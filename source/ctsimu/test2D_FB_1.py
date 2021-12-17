# -*- coding: UTF-8 -*-
from .test import *
from .general import *
import json
import pkgutil

class Test2D_FB_1_results:
    """ Results for one sub test of the noise scenario. """

    def __init__(self):
        self._mean = 0
        self._stdDev = 0
        self._snr = 0
        self._fwhm = 0

        # Nominal values from JSON file:
        self._SNRmax = None
        self._nominalMean = None
        self._nominalFWHM = None
        self._nominalSigma = None

        # Grey Value Histogram:
        self._gvHistogram = None  # stores the GV histogram
        self._gvHistNormalized = None  # area-normalized GV histogram


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

        self._imax = None    # maximum (noise-free) grey value in free beam, as required by the scenario
        self._gvMax = 70000  # maximum grey value for the histograms
        self._xValues = None # x grey values for the histograms

        # Results for each sub test:
        self._results = []

    def prepare(self):
        """ Preparations before the test will be run with the images from the pipeline. """       
        if not isinstance(self._pipe, ProcessingPipeline):
            self._prepared = False
            raise Exception("Step must be part of a processing pipeline before it can prepare. Current pipeline: {}".format(self._pipe))

        if not self._prepared:
            self._results = []
            self._xValues = numpy.linspace(0, self._gvMax, self._gvMax, endpoint=False)

            if self._jsonScenarioFile != None:
                self._prepared = True
            else:
                raise Exception("Test {name}: Please provide a JSON scenario description.".format(name=self._name))

            self._prepared = True

    def prepareRun(self, i):
        if i < len(self._subtests):
            if self._subtests[i] == "SNR100":
                self._jsonScenarioFile = "scenarios/2D-FB-1_Detektor1_SNR100_2021-05-25v06r00dp.json"
            elif self._subtests[i] == "SNR250":
                self._jsonScenarioFile = "scenarios/2D-FB-1_Detektor2_SNR250_2021-05-25v06r00dp.json"
            else:
                raise Exception("{key} is not a valid subtest identifier for test scenario {test}".format(key=self._subtests[i], test=self._testName))

            # Get Imax and the SNRmax from the JSON file:
            jsonText = pkgutil.get_data(__name__, self._jsonScenarioFile).decode()

            if jsonText != None:
                jsonDict = json.loads(jsonText)

                self._imax = getFieldOrNone(jsonDict, "detector", "grey_value", "imax", "value")
                snrMax = getFieldOrNone(jsonDict, "detector", "noise", "snr_at_imax", "value")

                if self._imax == None:
                    raise Exception("Test {name}: Cannot find 'imax' value in JSON scenario description: {json}".format(name=self._name, json=self._jsonScenarioFile))

                if snrMax == None:
                    raise Exception("Test {name}: Cannot find 'snr_at_imax' value in JSON scenario description: {json}".format(name=self._name, json=self._jsonScenarioFile))

                nominalSigma = self._imax / snrMax
                nominalFWHM  = 2.0*math.sqrt(2.0*math.log(2.0))*self._imax / snrMax

            else:
                raise Exception("Test {name}: Cannot open JSON scenario description: {json}".format(name=self._name, json=self._jsonScenarioFile))

            # Prepare this test run:
            self.prepare()
            self._results.append(Test2D_FB_1_results())
            self._results[i]._gvHistogram = numpy.zeros(self._gvMax, dtype=numpy.uint32)
            self._results[i]._SNRmax = snrMax
            self._results[i]._nominalMean = self._imax # after ff-correction
            self._results[i]._nominalFWHM = nominalFWHM
            self._results[i]._nominalSigma = nominalSigma
           
        else:
            if len(self._subtests) == 0:
                raise Exception("Please provide keywords that identify which metadata file belongs to which subtest. Test {testname} accepts two keywords: 'SNR100' and 'SNR250'.".format(testname=self._testName))
            else:
                raise Exception("Number of provided image metadata files exceeds number of test runs ({expected}).".format(expected=len(self._subtests)))


    def run(self, image):
        if self._currentRun < self._nExpectedRuns:
            self.prepareRun(self._currentRun)
            i = self._currentRun
            subtestName = self._subtests[i]  

            # Create Grey Value Histogram:
            for x in range(image.getWidth()):
                print("\rCalculating grey value histogram... {:.1f}%".format(100*x/image.getWidth()), end="")
                for y in range(image.getHeight()):
                    gv = int(image.getPixel(x, y))
                    if gv >= 0 and gv < self._gvMax:
                        self._results[i]._gvHistogram[gv] += 1

            print("\rCalculating grey value histogram... 100%  ")

            self._results[i]._mean   = image.mean()
            self._results[i]._stdDev = image.stdDev()
            self._results[i]._snr    = self._results[i]._mean / self._results[i]._stdDev
            self._results[i]._fwhm   = 2.0*math.sqrt(2.0*math.log(2.0)) * self._results[i]._stdDev

            self._results[i]._gvHistNormalized = self._results[i]._gvHistogram / numpy.sum(self._results[i]._gvHistogram)
            
 

            log("Writing evaluation results...")

            # Write evaluation text files:
            result = self._results[i]

            summary  = "# Evaluation of test {name}, {subname}:\n".format(name=self._name, subname=subtestName)
            summary += "#\n"
            summary += "# Signal-to-noise ratio (SNR)\n"
            summary += "# -----------------------------------------\n"
            summary += "# Overall projection SNR:   {snr:.3f}\n".format(snr=result._snr)
            summary += "# Nominal SNR:              {nominalSNR:.3f}\n".format(nominalSNR=result._SNRmax)
            summary += "# Relative SNR deviation:   {reldev:.5f}\n".format(reldev=(result._snr-result._SNRmax)/result._SNRmax)
            summary += "#\n"

            summary += "# Properties of grey value distribution\n"
            summary += "# -----------------------------------------\n"
            summary += "# Mean grey value:          {meanGV:.3f}\n".format(meanGV=result._mean)
            summary += "# Nominal mean grey value:  {nominalMeanGV:.3f}\n".format(nominalMeanGV=result._nominalMean)
            summary += "# Relative mean deviation:  {reldev:.5f}\n".format(reldev=(result._mean-result._nominalMean)/self._imax)
            summary += "#\n"
            summary += "# Measured RMSD:            {stdDev:.3f}\n".format(stdDev=result._stdDev)
            summary += "# Nominal RMSD:             {nominalStdDev:.3f}\n".format(nominalStdDev=result._nominalSigma)
            summary += "# Relative RMSD deviation:  {reldev:.5f}\n".format(reldev=(result._stdDev-result._nominalSigma)/result._nominalSigma)
            summary += "#\n"
            summary += "# Measured noise FWHM:      {fwhm:.3f}\n".format(fwhm=result._fwhm)
            summary += "# Nominal FWHM:             {nominalFWHM:.3f}\n".format(nominalFWHM=result._nominalFWHM)
            summary += "# Relative FWHM deviation:  {reldev:.5f}\n".format(reldev=(result._fwhm-result._nominalFWHM)/result._nominalFWHM)

            resultFileName = "{dir}/{name}_{subname}_summary.txt".format(dir=self._resultFileDirectory, name=self._name, subname=subtestName)
            with open(resultFileName, 'w') as resultFile:
                resultFile.write(summary)
                resultFile.close()


            histogram  = "# Grey value histogram\n"
            histogram += "# -----------------------------------------\n"
            histogram += "# GV\tcounts\tprobability density\n"

            for i in range(len(result._gvHistogram)):
                histogram += "{GV}\t{counts}\t{probDensity}\n".format(GV=i, counts=result._gvHistogram[i], probDensity=result._gvHistNormalized[i])

            histogramFileName = "{dir}/{name}_{subname}_histogram.txt".format(dir=self._resultFileDirectory, name=self._name, subname=subtestName)
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
                maxNominalSigma = result._nominalSigma

                xStart = self._imax - 4*maxNominalSigma
                xStop  = self._imax + 4*maxNominalSigma

                matplotlib.use("agg")

                fig, ax = matplotlib.pyplot.subplots(nrows=1, ncols=1, figsize=(9, 7))
                
                ax.plot(self._xValues, result._gvHistNormalized, 'o', markersize=1.0, label="Measured (Subtest {subtest})".format(subtest=subtestName), color='#1f77b4')

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

                plotFilename = "{dir}/{name}_{subname}_histogram.png".format(dir=self._resultFileDirectory, name=self._name, subname=subtestName)
                matplotlib.pyplot.savefig(plotFilename)
                fig.clf()
                matplotlib.pyplot.close('all')
            except:
                log("Warning: Error plotting results for test {name}, {subname} using matplotlib.".format(name=self._name, subname=subtestName))

            log("Evaluation data for test {name}, {subtest} written to {dir}.".format(name=self._name, subtest=subtestName, dir=self._resultFileDirectory))


            self._currentRun += 1

            # Return image to the pipeline
            return image
        else:
            raise Exception("Number of provided image metadata files exceeds expected number of test runs ({expected}).".format(expected=self._nExpectedRuns))

    def followUp(self):
        pass
        