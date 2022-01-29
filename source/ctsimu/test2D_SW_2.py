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

        # Grey values from Monte-Carlo simuation:
        self._means_rosi_total    = None   # accounting for scatter radiation
        self._means_rosi_primary  = None   # accounting only for primary radiation
        self._means_mcray_total    = None   # accounting for scatter radiation
        self._means_mcray_primary  = None   # accounting only for primary radiation

        # Reference values for Monte-Carlo simulation (calculated in loadReference())
        self._means_mc_total    = None   # accounting for scatter radiation
        self._means_mc_primary  = None   # accounting only for primary radiation
        self._error_mc_total_upper   = None
        self._error_mc_total_lower   = None
        self._error_mc_primary_upper = None
        self._error_mc_primary_lower = None

        # reference grey value ratios
        self._ratios_mc_total = None
        self._ratios_mc_primary = None
        self._error_ratios_mc_total_upper = []
        self._error_ratios_mc_total_lower = []
        self._error_ratios_mc_primary_upper = []
        self._error_ratios_mc_primary_lower = []

    def loadReference(self, name):
        dataText = pkgutil.get_data(__name__, "data/2D-SW-2_scenario{name}.txt".format(name=name)).decode()
        dataIO = io.StringIO(dataText)
        allData = numpy.loadtxt(dataIO, delimiter='\t')  # ignore free beam
        
        # --- total radiation
        means_rosi_total  = allData[:,1]
        means_mcray_total = allData[:,4]
        means_total = (means_rosi_total + means_mcray_total) * 0.5
        delta_total = numpy.absolute(means_rosi_total - means_mcray_total)

        stddev_rosi_total  = allData[:,2]
        stddev_mcray_total = allData[:,5]

        #err_total = 0.5*(delta_total + numpy.sqrt(numpy.square(stddev_rosi_total) + numpy.square(stddev_mcray_total)))
        err_total_rosi  = numpy.fmax(stddev_rosi_total,  stddev_mcray_total-delta_total) + 0.5*delta_total
        err_total_mcray = numpy.fmax(stddev_mcray_total, stddev_rosi_total-delta_total)  + 0.5*delta_total
        err_total_upper = numpy.zeros_like(err_total_rosi)
        err_total_lower = numpy.zeros_like(err_total_rosi)
        for i in range(len(err_total_upper)):
            if means_rosi_total[i] > means_total[i]:  # ROSI is upper bound
                err_total_upper[i] = err_total_rosi[i]
                err_total_lower[i] = err_total_mcray[i]
            else:   # McRay is upper bound
                err_total_upper[i] = err_total_mcray[i]
                err_total_lower[i] = err_total_rosi[i]


        # --- primary radiation
        means_rosi_primary  = allData[:,9]
        means_mcray_primary = allData[:,12]
        means_primary = (means_rosi_primary + means_mcray_primary) * 0.5
        delta_primary = numpy.absolute(means_rosi_primary - means_mcray_primary)

        stddev_rosi_primary  = allData[:,10]
        stddev_mcray_primary = allData[:,13]

        #err_primary = 0.5*(delta_primary + numpy.sqrt(numpy.square(stddev_rosi_primary) + numpy.square(stddev_mcray_primary)))
        err_primary_rosi  = numpy.fmax(stddev_rosi_primary,  stddev_mcray_primary-delta_primary) + 0.5*delta_primary
        err_primary_mcray = numpy.fmax(stddev_mcray_primary, stddev_rosi_primary-delta_primary)  + 0.5*delta_primary
        err_primary_upper = numpy.zeros_like(err_primary_rosi)
        err_primary_lower = numpy.zeros_like(err_primary_rosi)
        for i in range(len(err_primary_upper)):
            if means_rosi_primary[i] > means_primary[i]:  # ROSI is upper bound
                err_primary_upper[i] = err_primary_rosi[i]
                err_primary_lower[i] = err_primary_mcray[i]
            else:   # McRay is upper bound
                err_primary_upper[i] = err_primary_mcray[i]
                err_primary_lower[i] = err_primary_rosi[i]


        self._means_rosi_total    = means_rosi_total
        self._means_rosi_primary  = means_rosi_primary
        self._means_mcray_total    = means_mcray_total
        self._means_mcray_primary  = means_mcray_primary 
        
        self._means_mc_total         = means_total
        self._means_mc_primary       = means_primary
        self._error_mc_total_upper   = err_total_upper
        self._error_mc_total_lower   = err_total_lower
        self._error_mc_primary_upper = err_primary_upper
        self._error_mc_primary_lower = err_primary_lower

        dataIO.close()

        self._ratios_mc_total   = ratios(self._means_mc_total)
        self._ratios_mc_primary = ratios(self._means_mc_primary)

        # calculate maximum uncertainties of MC ratios:
        for v in range(1, len(self._ratios_mc_primary)+1):
            c, uncertainty_total_upper = divideAndError(
                muA = self._means_mc_total[v-1],
                muB = self._means_mc_total[v],
                errA = self._error_mc_total_upper[v],
                errB = self._error_mc_total_upper[v-1]
                )
            c, uncertainty_primary_upper = divideAndError(
                muA = self._means_mc_primary[v-1],
                muB = self._means_mc_primary[v],
                errA = self._error_mc_primary_upper[v],
                errB = self._error_mc_primary_upper[v-1]
                )

            c, uncertainty_total_lower = divideAndError(
                muA = self._means_mc_total[v-1],
                muB = self._means_mc_total[v],
                errA = self._error_mc_total_lower[v],
                errB = self._error_mc_total_lower[v-1]
                )
            c, uncertainty_primary_lower = divideAndError(
                muA = self._means_mc_primary[v-1],
                muB = self._means_mc_primary[v],
                errA = self._error_mc_primary_lower[v],
                errB = self._error_mc_primary_lower[v-1]
                )

            self._error_ratios_mc_total_upper.append(uncertainty_total_upper)
            self._error_ratios_mc_primary_upper.append(uncertainty_primary_upper)
            self._error_ratios_mc_total_lower.append(uncertainty_total_lower)
            self._error_ratios_mc_primary_lower.append(uncertainty_primary_lower)

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
        ratioText += "# step A\tstep B\tratio (A/B)\treference ratio (primary)\terror_lower (ref. primary)\terror_upper (ref. primary)\trel. deviation (primary)\treference ratio (total)\terror_lower (ref. total)\terror_upper (ref. total)\trel. deviation (total)\n"

        for r in range(len(self._results[i]._ratios)):
            ratioText += "{A}\t{B}\t{ratio:.5f}\t{refPrimary:.5f}\t{errorPrimaryLower:.5f}\t{errorPrimaryUpper:.5f}\t{devFacPrimary:.5f}\t{reftotal:.5f}\t{errorTotalLower:.5f}\t{errorTotalUpper:.5f}\t{devFacTotal:.5f}\n".format(
                    A = (r+2),
                    B = (r+1),
                    ratio = self._results[i]._ratios[r],
                    refPrimary = self._results[i]._ratios_mc_primary[r],
                    errorPrimaryLower = self._results[i]._error_ratios_mc_primary_lower[r],
                    errorPrimaryUpper = self._results[i]._error_ratios_mc_primary_upper[r],
                    devFacPrimary = (self._results[i]._ratios[r] / self._results[i]._ratios_mc_primary[r] - 1),
                    reftotal = self._results[i]._ratios_mc_total[r],
                    errorTotalLower = self._results[i]._error_ratios_mc_total_lower[r],
                    errorTotalUpper = self._results[i]._error_ratios_mc_total_upper[r],
                    devFacTotal = (self._results[i]._ratios[r] / self._results[i]._ratios_mc_total[r] - 1)
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
        xLabels = ("2 by 1", "3 by 2", "4 by 3", "5 by 4", "6 by 5", "7 by 6", "8 by 7", "9 by 8", "10 by 9")

        try:
            import matplotlib
            import matplotlib.pyplot
            from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

            matplotlib.use("agg")

            for mode in ("primary", "total"):             
                fig, ax = matplotlib.pyplot.subplots(nrows=1, ncols=1, figsize=(8, 5))
                
                # Grey Value Profile:
                if mode == "primary":
                    modeDescription = "primary radiation"
                    ax.errorbar(xValues, self._results[i]._ratios_mc_primary, xerr=None, yerr=[self._results[i]._error_ratios_mc_primary_lower, self._results[i]._error_ratios_mc_primary_upper], linewidth=0, elinewidth=2.0, ecolor='#fe6100')
                    ax.plot(xValues, self._results[i]._ratios_mc_primary, '_', markersize=11.0, label="Monte-Carlo reference", color='#fe6100')
                else:
                    modeDescription = "total radiation (scatter+primary)"
                    ax.errorbar(xValues, self._results[i]._ratios_mc_total, xerr=None, yerr=[self._results[i]._error_ratios_mc_total_lower, self._results[i]._error_ratios_mc_total_upper], linewidth=0, elinewidth=2.0, ecolor='#fe6100')
                    ax.plot(xValues, self._results[i]._ratios_mc_total, '_', markersize=11.0, label="Monte-Carlo reference", color='#fe6100')

                ax.plot(xValues, self._results[i]._ratios, 'o', markersize=5.0, label="measured", color='#1f77b4')

                ax.set_xlabel("step pair division")
                ax.set_ylabel("grey value ratio")
                ax.set_title("2D-SW-2, {sub}, {details}".format(sub=self._results[i]._longName, details=modeDescription), loc="left", fontsize=10)
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
