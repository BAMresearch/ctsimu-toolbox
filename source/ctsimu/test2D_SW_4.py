# -*- coding: UTF-8 -*-
from .test import *
from .general import *

import pkgutil
import io

class Test2D_SW_4_results:
    """ Results for one sub test of the scenario. """

    def __init__(self):
        self._profiles            = []
        self._profiles_nominal    = []
        self._profiles_fitpars    = []  # fit parameters
        self._profiles_nominalfit = []
        self._profiles_pos        = None
        

class Test2D_SW_4(generalTest):
    """ CTSimU test 2D-SW-4: scattering. """

    def __init__(self, resultFileDirectory=".", name=None, rawOutput=False):
        generalTest.__init__(
            self,
            testName="2D-SW-4",
            name=name,
            resultFileDirectory=resultFileDirectory,
            rawOutput=rawOutput)

        # Profile Data:
        self._results = Test2D_SW_4_results()

        # Definitions of the central and right ROI for the vertical grey value profiles:
        self._steps = [
            ImageROI(101, 0, 111, 200),
            ImageROI(190, 0, 200, 200)
        ]

    def prepare(self):
        """ Preparations before the test will be run with the images from the pipeline. """
        if not self._prepared:
            # Load the nominal profile data:
            profilesText = pkgutil.get_data(__name__, "data/2D-SW-4_profiles_monte-carlo_reference.txt").decode()

            profilesIO = io.StringIO(profilesText)


            profilesData = numpy.loadtxt(profilesIO, delimiter='\t')
            for i in range(len(self._steps)):
                self._results._profiles_nominal.append(profilesData[:,(i+1)])

            profilesIO.close()


            # Vertical profile positions:
            nP = len(self._results._profiles_nominal[0])
            self._results._profiles_pos = numpy.linspace(0, nP, nP, endpoint=False)


            # Fit a fourth order polynomial to the MC data to smooth it:
            for i in range(len(self._steps)):
                popt, pcov = optimize.curve_fit(f=poly4, xdata=self._results._profiles_pos, ydata=self._results._profiles_nominal[i], p0=None)
                a = popt[0]
                b = popt[1]
                c = popt[2]
                d = popt[3]
                e = popt[4]
                self._results._profiles_nominalfit.append(poly4(x=self._results._profiles_pos, a=a, b=b, c=c, d=d, e=e))

        self._prepared = True

    def writeProfiles(self, title, filename, pos, profiles):
        profilesText  = "# {}\n".format(title)
        profilesText += "# Mean vertical grey value profiles\n"
        profilesText += "# pos [px]\tcentral [GV]\tright [GV]\n"

        for i in range(len(pos)):
            profilesText += "{}".format(int(pos[i]))

            for p in profiles:
                profilesText += "\t{}".format(p[i])

            profilesText += "\n"

        with open(filename, 'w') as profilesFile:
            profilesFile.write(profilesText)
            profilesFile.close()

    def run(self, image):
        self.prepare()

        # Vertical profiles:
        for roi in self._steps:
            profile = image.verticalROIProfile(roi)
            self._results._profiles.append(profile)   

        # Write profiles
        profileFileName = "{dir}/{name}_profiles_measured.txt".format(dir=self._resultFileDirectory, name=self._name)
        self.writeProfiles(title="Measured profiles", filename=profileFileName, pos=self._results._profiles_pos, profiles=self._results._profiles)

        profileFileName = "{dir}/{name}_profiles_monte-carlo.txt".format(dir=self._resultFileDirectory, name=self._name)
        self.writeProfiles(title="Monte Carlo results as basis for reference values", filename=profileFileName, pos=self._results._profiles_pos, profiles=self._results._profiles_nominal)
        
        profileFileName = "{dir}/{name}_profiles_reference.txt".format(dir=self._resultFileDirectory, name=self._name)
        self.writeProfiles(title="Reference values (fourth order polynomial fit to Monte-Carlo measurements)", filename=profileFileName, pos=self._results._profiles_pos, profiles=self._results._profiles_nominalfit)

        log("Evaluation data for test {name} written to {dir}.".format(name=self._name, dir=self._resultFileDirectory))

        self.plotResults()

        self._currentRun += 1
        return image

    def followUp(self):
        pass

    def plotResults(self):
        try:
            import matplotlib
            import matplotlib.pyplot
            from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

            xValues = self._results._profiles_pos
            tolerance = 50

            matplotlib.use("agg")

            fig, ((ax1, ax2)) = matplotlib.pyplot.subplots(nrows=2, ncols=1, figsize=(8, 8))
          
            # Grey Value Profiles:
            #ax1.fill_between(xValues, self._results._profiles_nominalfit[0]+tolerance, self._results._profiles_nominalfit[0]-tolerance, facecolor='#ffcc00')
            ax1.plot(xValues, self._results._profiles_nominalfit[0], linewidth=5.0, linestyle="dotted", label="", color='#c18100')
            ax1.plot(xValues, self._results._profiles[0], linewidth=2.0, label="", color='#1f77b4')
            ax1.set_xlabel("Vertical distance in px")
            ax1.set_ylabel("Grey value")
            ax1.set_title("Central profile")
            ax1.xaxis.set_major_locator(MultipleLocator(20))
            ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
            ax1.xaxis.set_minor_locator(MultipleLocator(10))
            ax1.grid(b=True, which='major', axis='both', color='#d9d9d9', linestyle='dashed')
            ax1.grid(b=True, which='minor', axis='both', color='#e7e7e7', linestyle='dotted')

            #ax2.fill_between(xValues, self._results._profiles_nominalfit[1]+tolerance, self._results._profiles_nominalfit[1]-tolerance, facecolor='#ffcc00')
            ax2.plot(xValues, self._results._profiles_nominalfit[1], linewidth=5.0, linestyle="dotted", label="", color='#c18100')
            ax2.plot(xValues, self._results._profiles[1], linewidth=2.0, label="", color='#1f77b4')
            ax2.set_xlabel("Vertical distance in px")
            ax2.set_ylabel("Grey value")
            ax2.set_title("Right border profile")
            ax2.xaxis.set_major_locator(MultipleLocator(20))
            ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))
            ax2.xaxis.set_minor_locator(MultipleLocator(10))
            ax2.grid(b=True, which='major', axis='both', color='#d9d9d9', linestyle='dashed')
            ax2.grid(b=True, which='minor', axis='both', color='#e7e7e7', linestyle='dotted')

            fig.tight_layout()

            plotFilename = "{dir}/{name}_profiles.png".format(dir=self._resultFileDirectory, name=self._name)
            matplotlib.pyplot.savefig(plotFilename)
            fig.clf()
            matplotlib.pyplot.close('all')
        except:
            log("Warning: Error plotting results for test {name} using matplotlib.".format(name=self._name))
