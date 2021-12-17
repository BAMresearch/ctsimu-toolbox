# -*- coding: UTF-8 -*-
from .test import *
from .general import *

import pkgutil
import io

class Test2D_SW_2_results:
    """ Results for one sub test of the filtering scenario. """

    def __init__(self):
        self._longName = ""

        # Grey value means and ratios (per wedge step):
        self._means             = []     # Mean value for each step
        self._ratios            = None   # Grey value ratios between the steps

        # reference values from Monte-Carlo simulation
        self._means_mc_total  = None   # accounting for scatter radiation
        self._means_mc_primary  = None   # accounting only for primary radiation
        self._stddev_mc_total = None   # standard deviations from Monte-Carlo simulations
        self._stddev_mc_primary = None   # standard deviations from Monte-Carlo simulations

        # reference grey value ratios
        self._ratios_mc_total = None
        self._ratios_mc_primary = None
        self._ratios_stddev_mc_total = []
        self._ratios_stddev_mc_primary = []

    def loadReference(self, name):
        dataText = pkgutil.get_data(__name__, "data/2D-SW-2_scenario{name}.txt".format(name=name)).decode()
        dataIO = io.StringIO(dataText)
        allData = numpy.loadtxt(dataIO, delimiter='\t')  # ignore free beam
        
        means_rosi_total  = allData[:,1]
        means_mcray_total = allData[:,4]
        means_total = (means_rosi_total + means_mcray_total) * 0.5
        delta_total = numpy.absolute(means_rosi_total - means_mcray_total)

        stddev_rosi_total  = allData[:,2]
        stddev_mcray_total = allData[:,5]

        err_total = 0.5*(delta_total + numpy.sqrt(numpy.square(stddev_rosi_total) + numpy.square(stddev_mcray_total)))


        means_rosi_primary  = allData[:,9]
        means_mcray_primary = allData[:,12]
        means_primary = (means_rosi_primary + means_mcray_primary) * 0.5
        delta_primary = numpy.absolute(means_rosi_primary - means_mcray_primary)

        stddev_rosi_primary  = allData[:,10]
        stddev_mcray_primary = allData[:,13]

        err_primary = 0.5*(delta_primary + numpy.sqrt(numpy.square(stddev_rosi_primary) + numpy.square(stddev_mcray_primary)))
        
        self._means_mc_total = means_total
        self._means_mc_primary = means_primary
        self._stddev_mc_total = err_total
        self._stddev_mc_primary = err_primary

        dataIO.close()

        self._ratios_mc_total = ratios(self._means_mc_total)
        self._ratios_mc_primary = ratios(self._means_mc_primary)

        # calculate Gaussian uncertainties of MC ratios:
        for v in range(1, len(self._ratios_mc_primary)+1):
            c, uncertainty_total = divideAndError(
                muA = self._means_mc_total[v-1],
                muB = self._means_mc_total[v],
                sigmaA = self._stddev_mc_total[v],
                sigmaB = self._stddev_mc_total[v-1]
                )
            self._ratios_stddev_mc_total.append(uncertainty_total)

            c, uncertainty_primary = divideAndError(
                muA = self._means_mc_primary[v-1],
                muB = self._means_mc_primary[v],
                sigmaA = self._stddev_mc_primary[v],
                sigmaB = self._stddev_mc_primary[v-1]
                )
            self._ratios_stddev_mc_primary.append(uncertainty_primary)

class Test2D_SW_2(generalTest):
    """ CTSimU test 2D-SW-1: detector models / scintillators.
        CTSimU test 2D-SW-2: spectral filtering. """

    def __init__(self, resultFileDirectory=".", name=None, rawOutput=False):
        generalTest.__init__(
            self,
            testName="2D-SW-2",
            name=name,
            resultFileDirectory=resultFileDirectory,
            rawOutput=rawOutput)

        self._results = []

        self._shrink = 35
        self._leftOffset = 510 - self._shrink
        self._nPixels = 20 + 2*self._shrink

        # Absolute step definitions. Will be shrunk to accept tolerance border,
        # but these absolute definitions are needed for grey value rescaling and clipping.
        self._steps = [
            ImageROI(self._leftOffset, 819, self._leftOffset+self._nPixels, 910),
            ImageROI(self._leftOffset, 728, self._leftOffset+self._nPixels, 819),
            ImageROI(self._leftOffset, 637, self._leftOffset+self._nPixels, 728),
            ImageROI(self._leftOffset, 546, self._leftOffset+self._nPixels, 637),
            ImageROI(self._leftOffset, 455, self._leftOffset+self._nPixels, 546),
            ImageROI(self._leftOffset, 364, self._leftOffset+self._nPixels, 455),
            ImageROI(self._leftOffset, 273, self._leftOffset+self._nPixels, 364),
            ImageROI(self._leftOffset, 182, self._leftOffset+self._nPixels, 273),
            ImageROI(self._leftOffset,  91, self._leftOffset+self._nPixels, 182),
            ImageROI(self._leftOffset,   0, self._leftOffset+self._nPixels,  91)  # free beam
        ]

    def prepare(self):
        """ Preparations before the test will be run with the images from the pipeline. """
        self._prepared = True

    def prepareRun(self, i):
        if i < len(self._subtests):
            results = Test2D_SW_2_results()

            if self._subtests[i] == "Al":
                results.loadReference("01_Al_noFilter_200kV-poly_Ideal")
                results._longName = "Al wedge, no filter"
            elif self._subtests[i] == "AlCu":
                results.loadReference("02_Al_CuFilter_200kV-poly_Ideal")
                results._longName = "Al wedge, 2 mm Cu filter"
            elif self._subtests[i] == "Fe":
                results.loadReference("03_Fe_noFilter_200kV-poly_Ideal")
                results._longName = "Fe wedge, no filter"
            elif self._subtests[i] == "FeCu":
                results.loadReference("04_Fe_CuFilter_200kV-poly_Ideal")
                results._longName = "Fe wedge, 2 mm Cu filter"
            else:
                raise Exception("{key} is not a valid subtest identifier for test scenario {test}".format(key=self._subtests[i], test=self._testName))

            self._name = self._testName
            self._results.append(results)
        else:
            if len(self._subtests) == 0:
                raise Exception("Please provide keywords that identify which metadata file belongs to which subtest. Test {testname} accepts the following keywords: 'Al', 'AlCu', 'Fe' and 'FeCu'.".format(testname=self._testName))
            else:
                raise Exception("Number of provided image metadata files exceeds number of test runs ({expected}).".format(expected=len(self._subtests)))

    def run(self, image):
        self.prepare()
        self.prepareRun(self._currentRun)
        i = self._currentRun
        subtestName = self._subtests[i]  

        # Grey value summary        
        statsText = "# Evaluation of Test {name}, {subname}:\n".format(name=self._name, subname=subtestName)
        statsText += "# {longDesc}\n".format(longDesc=self._results[i]._longName)
        statsText += "# \n"        
        statsText += "# ROI mean grey value per step\n"
        statsText += "# step\tx0\ty0\tx1\ty1\twidth [px]\theight [px]\tarea [px]\tmean [GV]\n"

        step = 0
        for roi in self._steps:
            step += 1
            smallerROI = copy.deepcopy(roi)
            smallerROI.grow(-self._shrink)
            stats = image.stats(smallerROI)

            statsText += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\n".format(step, smallerROI.x0(), smallerROI.y0(), smallerROI.x1(), smallerROI.y1(), stats["width"], stats["height"], stats["area"], stats["mean"])

            self._results[i]._means.append(stats["mean"])

        self._results[i]._ratios = ratios(self._results[i]._means)      

        statsFileName = "{dir}/{name}_{subname}_grey_values.txt".format(dir=self._resultFileDirectory, name=self._name, subname=subtestName)
        with open(statsFileName, 'w') as statsFile:
            statsFile.write(statsText)
            statsFile.close()

        # Ratio summary        
        ratioText = "# Evaluation of Test {name}, {subname}:\n".format(name=self._name, subname=subtestName)
        ratioText += "# {longDesc}\n".format(longDesc=self._results[i]._longName)
        ratioText += "# \n"    

        ratioText += "# Grey value ratios\n"
        ratioText += "# step A\tstep B\tratio (A/B)\treference ratio (primary)\tuncertainty (ref. primary)\trel. deviation (primary)\treference ratio (total)\tuncertainty (ref. total)\trel. deviation (total)\n"

        for r in range(len(self._results[i]._ratios)):
            ratioText += "{A}\t{B}\t{ratio:.5f}\t{refPrimary:.5f}\t{stddevPrimary:.5f}\t{devFacPrimary:.5f}\t{reftotal:.5f}\t{stddevtotal:.5f}\t{devFacSctr:.5f}\n".format(
                    A = (r+2),
                    B = (r+1),
                    ratio = self._results[i]._ratios[r],
                    refPrimary = self._results[i]._ratios_mc_primary[r],
                    stddevPrimary = self._results[i]._ratios_stddev_mc_primary[r],
                    devFacPrimary = (self._results[i]._ratios[r] / self._results[i]._ratios_mc_primary[r] - 1),
                    reftotal = self._results[i]._ratios_mc_total[r],
                    stddevtotal = self._results[i]._ratios_stddev_mc_total[r],
                    devFacSctr = (self._results[i]._ratios[r] / self._results[i]._ratios_mc_total[r] - 1)
                )


        ratioFileName = "{dir}/{name}_{subname}_ratios.txt".format(dir=self._resultFileDirectory, name=self._name, subname=subtestName)
        with open(ratioFileName, 'w') as ratiosFile:
            ratiosFile.write(ratioText)
            ratiosFile.close()

        self.plotResults()

        self._currentRun += 1
        return image

    def followUp(self):
        pass

    def plotResults(self):
        i = self._currentRun
        subtestName = self._subtests[i]
        xValues = numpy.linspace(0, len(self._results[i]._ratios)-1, len(self._results[i]._ratios), endpoint=True)
        xLabels = ("2/1", "3/2", "4/3", "5/4", "6/5", "7/6", "8/7", "9/8", "10/9")

        try:
            import matplotlib
            import matplotlib.pyplot
            from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

            matplotlib.use("agg")

            for mode in ("primary", "total"):             
                fig, ax = matplotlib.pyplot.subplots(nrows=1, ncols=1, figsize=(9, 7))
                
                # Grey Value Profile:
                if mode == "primary":
                    modeDescription = "primary (without scatter radiation)"
                    ax.errorbar(xValues, self._results[i]._ratios_mc_primary, xerr=None, yerr=self._results[i]._ratios_stddev_mc_primary, linewidth=0, elinewidth=2.0, ecolor='#ffaa00')
                    ax.plot(xValues, self._results[i]._ratios_mc_primary, 'x', markersize=7.0, label="Monte-Carlo reference ± 1u", color='#ffaa00')
                else:
                    modeDescription = "total (accounting for scatter radiation)"
                    ax.errorbar(xValues, self._results[i]._ratios_mc_total, xerr=None, yerr=self._results[i]._ratios_stddev_mc_total, linewidth=0, elinewidth=2.0, ecolor='#ffaa00')
                    ax.plot(xValues, self._results[i]._ratios_mc_total, 'x', markersize=7.0, label="Monte-Carlo reference ± 1u", color='#ffaa00')

                ax.plot(xValues, self._results[i]._ratios, 'o', markersize=4.0, label="measured", color='#1f77b4')

                ax.set_xlabel("step pair")
                ax.set_ylabel("grey value ratio")
                ax.set_title("2D-SW-2, {sub} ({details})".format(sub=self._results[i]._longName, details=modeDescription))
                ax.set_xticks(xValues)
                ax.xaxis.set_ticklabels(xLabels)
                ax.grid(b=True, which='major', axis='both', color='#d9d9d9', linestyle='dashed')
                ax.grid(b=True, which='minor', axis='both', color='#e7e7e7', linestyle='dotted')
                ax.legend(loc='lower left')

                fig.tight_layout(pad=2.5)

                plotFilename = "{dir}/{name}_{subname}_ratios_{mode}.png".format(dir=self._resultFileDirectory, name=self._name, subname=subtestName, mode=mode)
                matplotlib.pyplot.savefig(plotFilename)
                fig.clf()
                matplotlib.pyplot.close('all')

        except:
            log("Warning: Error plotting results for test {name}, {subname} using matplotlib.".format(name=self._name, subname=subtestName))
