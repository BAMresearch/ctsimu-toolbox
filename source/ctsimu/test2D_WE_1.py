# -*- coding: UTF-8 -*-
from .test import *
from .general import *
from .geoprimitives import *
import json
import pkgutil
from scipy import optimize

class Test2D_WE_1_results:
    """ Results for one sub test. """

    def __init__(self):
        # Profile Data:
        self._lineProfilePos      = []
        self._lineProfileGV       = []  # this will be the edge spread function (ESF)
        self._lineProfileSmoothed = []
        self._LSF                 = []  # line spread function; derivative of ESF.
        self._LSF_windowed        = []  # von-Hann window applied before FFT
        self._LSF_gaussian_fit    = []
        self._MTF                 = []
        self._MTFpos              = []
        self._fnyq                = 0   # Nyquist frequency

        self._ideal_ESF_Convolution = None      # Convolution of Gaussian with SampledESF

        self._ideal_LSF_Gaussian    = None
        self._ideal_LSF_Convolution = None

        self._ideal_MTF_Gaussian       = None   # MTF of Gaussian LSF
        self._ideal_MTF_ConvolutionFFT = None   # Ideal MTF by FFT of Convolution (Gaussian and SampledLSF)
        self._ideal_MTF_Multiplication = None   # Ideal MTF from SampledMTF * GaussianMTF; same result.

        # Fit results [px]
        self._nominalGaussianSigmaPX       = 0

        self._idealGaussianSigmaPX         = 0
        self._idealGaussianSigmaPXerror    = 0
        self._idealGaussianMuPX            = 0
        self._idealGaussianMuPXerror       = 0
        self._idealGaussianAmpPX           = 0
        self._idealGaussianAmpPXerror      = 0

        self._measuredGaussianSigmaPX      = 0
        self._measuredGaussianSigmaPXerror = 0
        self._measuredGaussianMuPX         = 0
        self._measuredGaussianMuPXerror    = 0
        self._measuredGaussianAmpPX        = 0
        self._measuredGaussianAmpPXerror   = 0

        # Fit results [mm]
        self._nominalGaussianSigmaMM       = 0

        self._idealGaussianSigmaMM         = 0
        self._idealGaussianSigmaMMerror    = 0
        self._idealGaussianMuMM            = 0
        self._idealGaussianMuMMerror       = 0
        self._idealGaussianAmpMM           = 0
        self._idealGaussianAmpMMerror      = 0

        self._measuredGaussianSigmaMM      = 0
        self._measuredGaussianSigmaMMerror = 0
        self._measuredGaussianMuMM         = 0
        self._measuredGaussianMuMMerror    = 0
        self._measuredGaussianAmpMM        = 0
        self._measuredGaussianAmpMMerror   = 0

        # MTF evaluation results
        self._MTF20freq_ideal    = 0
        self._MTF20freq_measured = 0

        #self._MTF10freq_ideal    = 0
        #self._MTF10freq_measured = 0

        self._pixelSize = 0

    def MTF20freq_absDev(self):
        return (self._MTF20freq_measured - self._MTF20freq_ideal)

    def MTF20freq_relDev(self):
        if self._MTF20freq_ideal != 0:
            return (self._MTF20freq_measured - self._MTF20freq_ideal) / self._MTF20freq_ideal

        return 0

    #def MTF10freq_absDev(self):
    #    return (self._MTF10freq_measured - self._MTF10freq_ideal)

    #def MTF10freq_relDev(self):
    #    if self._MTF10freq_ideal != 0:
    #        return (self._MTF10freq_measured - self._MTF10freq_ideal) / self._MTF10freq_ideal

        return 0

    def sigma_absDev(self):
        return (self._measuredGaussianSigmaPX - self._idealGaussianSigmaPX)

    def sigma_relDev(self):
        if self._idealGaussianSigmaPX != 0:
            return (self._measuredGaussianSigmaPX - self._idealGaussianSigmaPX) / self._idealGaussianSigmaPX

        return 0


class Test2D_WE_1(generalTest):
    """ CTSimU test 2D-WE-1: focal spot size. """

    def __init__(self, resultFileDirectory=".", name=None, rawOutput=False):
        generalTest.__init__(
            self,
            testName="2D-WE-1",
            name=name,
            resultFileDirectory=resultFileDirectory,
            rawOutput=rawOutput)

        self._geometry = None
        self._analyticalIntensityProfileImage = None  # analytical flat field
        self._analyticalEdgeImage             = None  # stores the analytically computed edge image (sharp image)
        self._analyticalEdgeImageFF           = None  # stores the analytically computed edge image, flat-field corrected (sharp image)

        self._results = []
        self._results_idealPointSource = Test2D_WE_1_results()

        # Start point (x0,y0) and end point (x1,y1) for line profile:
        edgeAngle = 3 * (math.pi/180.0)    # 3 deg edge rotation

        # point of origin for edge profile (i.e. profile center):
        edgeOrigin = Vector2D(0, 1)                          # start with unit vector pointing "down"
        edgeOrigin.rotate(edgeAngle)                         # rotate by edge angle
        edgeOrigin.scale(500.5 / math.cos(edgeAngle) / 2.0)  # scale to half of visible edge length
        edgeOrigin = edgeOrigin + Vector2D(500.5, 500.5)     # move to detector centre

        #print("Edge origin:")
        #print(edgeOrigin)

        self._profileLength = 100.1   # pixels

        edgeDirection = Vector2D(self._profileLength/2, 0)
        edgeDirection.rotate(edgeAngle)

        self._p0 = edgeOrigin - edgeDirection
        self._p1 = edgeOrigin + edgeDirection

        #print("Line From p0:")
        #print(self._p0)
        #print("Line To p1:")
        #print(self._p1)

        self._profileWidth = 200  # pixels
        self._profileRes   = 0.1  # pixels



        # Also, prepare the clipping rectangle for the analytical
        # calculation of the ideal edge image. In pixel coordinates.
        A = Vector(   0,    0, 0)
        B = Vector(   0,  300, 0)
        C = Vector(-300,  300, 0)
        D = Vector(-300,    0, 0)

        rotAxis = Vector(0, 0, 1)
        #A.rotate(rotAxis, edgeAngle)
        B.rotate(rotAxis, edgeAngle)
        C.rotate(rotAxis, edgeAngle)
        D.rotate(rotAxis, edgeAngle)

        self._clippingRectangle = Polygon(A, B, C, D)


    def prepare(self):
        """ Preparations before the test will be run with the images from the pipeline. """
        if not isinstance(self._pipe, ProcessingPipeline):
            self._prepared = False
            raise Exception("Step must be part of a processing pipeline before it can prepare. Current pipeline: {}".format(self._pipe))

        if len(self._subtests) > 0:
            if not self._prepared:
                # It doesn't matter which of the sub-scenarios we take here.
                # They all have the same geometry.
                self._jsonScenarioFile = "scenarios/2D-WE-1_Spot1_2021-10-07v02r02rtr.json"

                if(self._jsonScenarioFile != None):
                    self._geometry = Geometry(jsonFileFromPkg=self._jsonScenarioFile)

                    print("Computing an analytical image for an ideal point source...")
                    self._analyticalIntensityProfileImage, self._analyticalEdgeImage = self._geometry.createDetectorFlatField_sphere(self._clippingRectangle)
                    self._analyticalEdgeImageFF = copy.deepcopy(self._analyticalEdgeImage)
                    self._analyticalEdgeImageFF.applyFlatfield(ref=self._analyticalIntensityProfileImage, rescaleFactor=60000.0)

                    # Raise analytical images to maximum grey value of 60000 before saving them.
                    # This rescaling does not affect the previous FF correction.
                    self._analyticalIntensityProfileImage.renormalize(newMin=0.0, newMax=60000.0, currentMin=0.0, currentMax=1.0)
                    self._analyticalEdgeImage.renormalize(newMin=0.0, newMax=60000.0, currentMin=0.0, currentMax=1.0)


                    # Write analytical images:
                    if self._rawOutput:
                        self._analyticalIntensityProfileImage.saveRAW("{dir}/{name}_ideal_flat.raw".format(dir=self._resultFileDirectory, name=self._name), dataType="float32", addInfo=True)
                        self._analyticalEdgeImage.saveRAW("{dir}/{name}_ideal_edge.raw".format(dir=self._resultFileDirectory, name=self._name), dataType="float32", addInfo=True)
                        self._analyticalEdgeImageFF.saveRAW("{dir}/{name}_ideal_edge_corrected.raw".format(dir=self._resultFileDirectory, name=self._name), dataType="float32", addInfo=True)
                    else: # TIFF
                        self._analyticalIntensityProfileImage.save("{dir}/{name}_ideal_flat.tif".format(dir=self._resultFileDirectory, name=self._name), dataType="float32")
                        self._analyticalEdgeImage.save("{dir}/{name}_ideal_edge.tif".format(dir=self._resultFileDirectory, name=self._name), dataType="float32")
                        self._analyticalEdgeImageFF.save("{dir}/{name}_ideal_edge_corrected.tif".format(dir=self._resultFileDirectory, name=self._name), dataType="float32")


                    print("Calculating the fundamental LSF and MTF from the ideal edge image...")
                    # Calculate the edge spread function (using a line profile across the edge):
                    self._results_idealPointSource._lineProfileGV, self._results_idealPointSource._lineProfilePos, stepsize = self._analyticalEdgeImageFF.lineProfile(x0=self._p0._x, y0=self._p0._y, x1=self._p1._x, y1=self._p1._y, width=self._profileWidth, resolution=self._profileRes)

                    # Nyquist frequency determines center of MTF frequency range:
                    self._results_idealPointSource._fnyq = 1.0 / (2.0*stepsize)
                    nSamples = len(self._results_idealPointSource._lineProfileGV)
                    #self._results_idealPointSource._MTFpos = numpy.linspace(start=0, stop=2*self._results_idealPointSource._fnyq, num=nSamples, endpoint=False, dtype=numpy.float64)

                    # Make the profile positions symmetric:
                    self._results_idealPointSource._lineProfilePos -= self._profileLength / 2

                    # Calculate the line spread function and MTF.
                    self._results_idealPointSource._lineProfileSmoothed, self._results_idealPointSource._LSF, self._results_idealPointSource._LSF_windowed, self._results_idealPointSource._MTF, self._results_idealPointSource._MTFpos = MTF(positions=self._results_idealPointSource._lineProfilePos, ESF=self._results_idealPointSource._lineProfileGV)

                    self._results_idealPointSource._lineProfileDelta = self._results_idealPointSource._lineProfileSmoothed - self._results_idealPointSource._lineProfileGV

                    self.writeResultFile(subname="fundamental_spotPoint", results=self._results_idealPointSource)

                    self._prepared = True
                else:
                    raise Exception("Test 2D-WE-1: Please provide a JSON scenario description.")

    def prepareRun(self, i):
        if i < len(self._subtests):
            self._jsonScenarioFile = "scenarios/2D-WE-1_Spot1_2021-10-07v02r02rtr.json"
            if self._subtests[i] == "spot1":
                self._jsonScenarioFile = "scenarios/2D-WE-1_Spot1_2021-10-07v02r02rtr.json"
            elif self._subtests[i] == "spot2":
                self._jsonScenarioFile = "scenarios/2D-WE-1_Spot2_2021-10-07v02r02rtr.json"
            else:
                raise Exception("{key} is not a valid subtest identifier for test scenario {test}".format(key=self._subtests[i], test=self._testName))

            results = Test2D_WE_1_results()
            
            # Get the Gaussian spot size and the pixel size from the JSON file:
            jsonText = pkgutil.get_data(__name__, self._jsonScenarioFile).decode()

            if jsonText != None:
                jsonDict = json.loads(jsonText)

                results._nominalGaussianSigmaMM = inMM(getFieldOrNone(jsonDict, "source", "spot", "sigma", "u"))
                results._pixelSize = inMM(getFieldOrNone(jsonDict, "detector", "pixel_pitch", "u"))

                results._nominalGaussianSigmaPX = results._nominalGaussianSigmaMM / results._pixelSize
            
            if self._subtests[i] == "spotPoint":
                results._nominalGaussianSigmaMM = 0
                results._nominalGaussianSigmaPX = 0
            else:
                # Create an ideal Gaussian LSF:
                nBins = len(self._results_idealPointSource._lineProfilePos)
                sigma = results._nominalGaussianSigmaPX
                results._ideal_LSF_Gaussian = numpy.zeros_like(a=self._results_idealPointSource._lineProfilePos, dtype=numpy.dtype('float64'))
                for j in range(nBins):
                    s = self._results_idealPointSource._lineProfilePos[j]
                    sLeft  = s - self._profileRes
                    sRight = s + self._profileRes
                    results._ideal_LSF_Gaussian[j] = math.erf(sRight/(math.sqrt(2.0)*sigma)) - math.erf(sLeft/(math.sqrt(2.0)*sigma)) #gaussian((s+sNext)/2.0, 0, sigma, 1)

                # Normalize ideal LSF maximum to 1:
                idealLSFmax = numpy.amax(results._ideal_LSF_Gaussian)
                if idealLSFmax != 0:
                    results._ideal_LSF_Gaussian /= idealLSFmax

                # Calculate convolution of ideal LSF of sharp image and Gaussian LSF:
                results._ideal_LSF_Convolution = numpy.convolve(a=self._results_idealPointSource._LSF, v=results._ideal_LSF_Gaussian, mode='same')
                convMax = numpy.amax(results._ideal_LSF_Convolution)
                if convMax != 0:
                    results._ideal_LSF_Convolution /= convMax

                # Calculate ideal ESF as convolution of Gaussian with the ideal sampledESF:
                results._ideal_ESF_Convolution = numpy.convolve(a=self._results_idealPointSource._lineProfileGV, v=results._ideal_LSF_Gaussian, mode='same')
                convMax = numpy.amax(results._ideal_ESF_Convolution)
                if convMax != 0:
                    results._ideal_ESF_Convolution /= convMax

                # Fit a Gaussian to central half of the data:
                fitParametersInitial = (0, sigma, 1)
                fitPositions = self._results_idealPointSource._lineProfilePos[int(nBins//4):int(3*nBins//4)]
                #fitData = results._ideal_LSF_Gaussian #[int(nBins//4):int(3*nBins//4)]
                #popt, pcov = optimize.curve_fit(f=gaussian, xdata=fitPositions, ydata=fitData, p0=fitParametersInitial)

                #print("- Ideal Gaussian Fit Mu:    {}".format(popt[0]))
                #print("- Ideal Gaussian Fit Sigma: {}".format(popt[1]))
                #print("- Ideal Gaussian Fit A:     {}".format(popt[2]))

                fitData = results._ideal_LSF_Convolution[int(nBins//4):int(3*nBins//4)]
                popt, pcov = optimize.curve_fit(f=gaussian, xdata=fitPositions, ydata=fitData, p0=fitParametersInitial)

                perr = numpy.sqrt(numpy.diag(pcov))

                results._idealGaussianSigmaPX      = popt[1]
                results._idealGaussianSigmaPXerror = perr[1]
                results._idealGaussianMuPX         = popt[0]
                results._idealGaussianMuPXerror    = perr[0]
                results._idealGaussianAmpPX        = popt[2]
                results._idealGaussianAmpPXerror   = perr[2]

                results._idealGaussianSigmaMM      = results._pixelSize * results._idealGaussianSigmaPX
                results._idealGaussianSigmaMMerror = results._pixelSize * results._idealGaussianSigmaPXerror
                results._idealGaussianMuMM         = results._pixelSize * results._idealGaussianMuPX
                results._idealGaussianMuMMerror    = results._pixelSize * results._idealGaussianMuPXerror
                results._idealGaussianAmpMM        = results._pixelSize * results._idealGaussianAmpPX
                results._idealGaussianAmpMMerror   = results._pixelSize * results._idealGaussianAmpPXerror

                print("- Analytical LSF Fit Mu:    {} +- {}".format(popt[0], perr[0]))
                print("- Analytical LSF Fit Sigma: {} +- {}".format(popt[1], perr[1]))
                print("- Analytical LSF Fit A:     {} +- {}".format(popt[2], perr[2]))


                # Calculate MTF for Gaussian LSF and Convolution LSF:
                smoothedESF, retLSF, retLSFsmoothed, results._ideal_MTF_Gaussian, pos = MTF(positions=self._results_idealPointSource._lineProfilePos, LSF=results._ideal_LSF_Gaussian)

                smoothedESF, retLSF, retLSFsmoothed, results._ideal_MTF_ConvolutionFFT, pos = MTF(positions=self._results_idealPointSource._lineProfilePos, LSF=results._ideal_LSF_Convolution)

                results._ideal_MTF_Multiplication = self._results_idealPointSource._MTF * results._ideal_MTF_Gaussian
                if results._ideal_MTF_Multiplication[0] != 0:
                    results._ideal_MTF_Multiplication /= results._ideal_MTF_Multiplication[0]


            print("Gaussian Spot Size: {} mm = {} px".format(results._nominalGaussianSigmaMM, results._nominalGaussianSigmaPX))

            self._results.append(results)
        else:
            if len(self._subtests) == 0:
                raise Exception("Please provide keywords that identify which metadata file belongs to which subtest. Test {testname} accepts the following keywords: 'spot1' and 'spot2'.".format(testname=self._testName))
            else:
                raise Exception("Number of provided image metadata files exceeds number of test runs ({expected}).".format(expected=len(self._subtests)))

    def writeSummaryFile(self, subname, results):
        summaryText  = "Evaluation results for {testname} {subtest}\n".format(testname=self._testName, subtest=subname)
        summaryText += "=================================================\n\n"

        summaryText += "Fit results for gauss(x) = A * exp(-(x-mu)²/(2*sigma²))\n"
        summaryText += "Errors: 1 standard deviation\n"
        summaryText += "-------------------------------------------------------\n"
        summaryText += "measured:   sigma [px] = {} +- {}\n".format(results._measuredGaussianSigmaPX, results._measuredGaussianSigmaPXerror)
        summaryText += "                  [mm] = {} +- {}\n".format(results._measuredGaussianSigmaMM, results._measuredGaussianSigmaMMerror)
        summaryText += "measured:   mu    [px] = {} +- {}\n".format(results._measuredGaussianMuPX,    results._measuredGaussianMuPXerror)
        summaryText += "                  [mm] = {} +- {}\n".format(results._measuredGaussianMuMM,    results._measuredGaussianMuMMerror)
        summaryText += "measured:   A     [px] = {} +- {}\n".format(results._measuredGaussianAmpPX,   results._measuredGaussianAmpPXerror)
        summaryText += "                  [mm] = {} +- {}\n\n".format(results._measuredGaussianAmpMM, results._measuredGaussianAmpMMerror)

        summaryText += "analytical: sigma [px] = {} +- {}\n".format(results._idealGaussianSigmaPX, results._idealGaussianSigmaPXerror)
        summaryText += "                  [mm] = {} +- {}\n".format(results._idealGaussianSigmaMM, results._idealGaussianSigmaMMerror)
        summaryText += "analytical: mu    [px] = {} +- {}\n".format(results._idealGaussianMuPX,    results._idealGaussianMuPXerror)
        summaryText += "                  [mm] = {} +- {}\n".format(results._idealGaussianMuMM,    results._idealGaussianMuMMerror)
        summaryText += "analytical: A     [px] = {} +- {}\n".format(results._idealGaussianAmpPX,   results._idealGaussianAmpPXerror)
        summaryText += "                  [mm] = {} +- {}\n\n".format(results._idealGaussianAmpMM, results._idealGaussianAmpMMerror)

        summaryText += "sigma: abs. deviation (measured-analytical) [px]       = {:.5f}\n".format(results.sigma_absDev())
        summaryText += "sigma: abs. deviation (measured-analytical) [mm]       = {:.5f}\n".format(results.sigma_absDev() * results._pixelSize)
        summaryText += "sigma: rel. deviation (measured-analytical)/analytical = {:.5f}\n\n".format(results.sigma_relDev())


        #summaryText += "MTF 10% frequency\n"
        #summaryText += "-------------------------------------------------------\n"
        #summaryText += "MTF10 measured: [cycles/px] = {:.3f}\n".format(results._MTF10freq_measured)
        #summaryText += "                [cycles/mm] = {:.3f}\n\n".format(results._MTF10freq_measured / results._pixelSize)

        #summaryText += "MTF10 ideal:    [cycles/px] = {:.3f}\n".format(results._MTF10freq_ideal)
        #summaryText += "                [cycles/mm] = {:.3f}\n\n".format(results._MTF10freq_ideal / results._pixelSize)

        #summaryText += "MTF10: abs. deviation (measured-ideal) [cycles/px] = {:.5f}\n".format(results.MTF10freq_absDev())
        #summaryText += "MTF10: abs. deviation (measured-ideal) [cycles/mm] = {:.5f}\n".format(results.MTF10freq_absDev() / results._pixelSize)
        #summaryText += "MTF10: rel. deviation (measured-ideal)/ideal       = {:.5f}\n".format(results.MTF10freq_relDev())


        summaryText += "\nMTF 20% frequency\n"
        summaryText += "-------------------------------------------------------\n"
        summaryText += "MTF20 measured:   [cycles/px] = {:.3f}\n".format(results._MTF20freq_measured)
        summaryText += "                  [cycles/mm] = {:.3f}\n\n".format(results._MTF20freq_measured / results._pixelSize)

        summaryText += "MTF20 analytical: [cycles/px] = {:.3f}\n".format(results._MTF20freq_ideal)
        summaryText += "                  [cycles/mm] = {:.3f}\n\n".format(results._MTF20freq_ideal / results._pixelSize)

        summaryText += "MTF20: abs. deviation (measured-analytical) [cycles/px] = {:.5f}\n".format(results.MTF20freq_absDev())
        summaryText += "MTF20: abs. deviation (measured-analytical) [cycles/mm] = {:.5f}\n".format(results.MTF20freq_absDev() / results._pixelSize)
        summaryText += "MTF20: rel. deviation (measured-analytical)/analytical  = {:.5f}\n".format(results.MTF20freq_relDev())


        summaryFileName = "{dir}/{name}_{subname}_summary.txt".format(dir=self._resultFileDirectory, name=self._name, subname=subname)
        with open(summaryFileName, 'w') as summaryFile:
            summaryFile.write(summaryText)
            summaryFile.close()


    def writeResultFile(self, subname, results):
        pos = results._lineProfilePos

        # ESF
        ESF         = results._lineProfileGV
        ESFsmoothed = results._lineProfileSmoothed
        ESFideal    = results._ideal_ESF_Convolution

        # LSF
        LSF         = results._LSF
        LSFwindowed = results._LSF_windowed
        LSFideal    = results._ideal_LSF_Convolution

        # MTF
        MTFpos      = results._MTFpos
        MTF         = results._MTF
        MTFideal    = results._ideal_MTF_ConvolutionFFT

        profileText  = "# Profile data: edge spread function (ESF), smoothed ESF, line spread function (LSF) and windowed LSF (von-Hann window)\n"
        profileText += "# s [px]\tESF\tESF_smoothed"
        if not(ESFideal is None):
            profileText += "\tESF_ideal"

        profileText += "\tLSF\tLSF_windowed"

        if not(LSFideal is None):
            profileText += "\tLSF_ideal"

        profileText += "\n"

        mtfText  = "# Modulation Transfer Function (MTF)\n"
        mtfText += "# f [cycles/px]\tMTF_measured"
        if not(MTFideal is None):
            mtfText += "\tMTF_ideal"

        mtfText += "\n"

        for j in range(len(pos)):
            profileText += "{pos:.2f}\t{ESF}\t{ESFsmoothed}".format(pos=pos[j], ESF=ESF[j], ESFsmoothed=ESFsmoothed[j])

            if not(ESFideal is None):
                profileText += "\t{ESFideal}".format(ESFideal=ESFideal[j])

            profileText += "\t{LSF}\t{LSFwindowed}".format(LSF=LSF[j], LSFwindowed=LSFwindowed[j])

            if not(LSFideal is None):
                profileText += "\t{LSFideal}".format(LSFideal=LSFideal[j])

            profileText += "\n"


            mtfText += "{f:.3f}\t{mtf}".format(f=MTFpos[j], mtf=MTF[j])

            if not(MTFideal is None):
                mtfText += "\t{MTFideal}".format(MTFideal=MTFideal[j])


            mtfText += "\n"
     
        profileFileName = "{dir}/{name}_{subname}_ESF_LSF.txt".format(dir=self._resultFileDirectory, name=self._name, subname=subname)
        with open(profileFileName, 'w') as profileFile:
            profileFile.write(profileText)
            profileFile.close()

        mtfFileName = "{dir}/{name}_{subname}_MTF.txt".format(dir=self._resultFileDirectory, name=self._name, subname=subname)
        with open(mtfFileName, 'w') as mtfFile:
            mtfFile.write(mtfText)
            mtfFile.close()

    def run(self, image):
        self.prepare()
        self.prepareRun(self._currentRun)
        i = self._currentRun
        subtestName = self._subtests[i]  

        print("Calculating the LSF and MTF of the projection image...")

        # Calculate the edge spread function (using a line profile across the edge):
        self._results[i]._lineProfileGV, self._results[i]._lineProfilePos, stepsize = image.lineProfile(x0=self._p0._x, y0=self._p0._y, x1=self._p1._x, y1=self._p1._y, width=self._profileWidth, resolution=self._profileRes)

        ESFmax = numpy.amax(self._results[i]._lineProfileGV)
        if ESFmax != 0:
            self._results[i]._lineProfileGV /= ESFmax

        # Nyquist frequency determines center of MTF frequency range:
        self._results[i]._fnyq = 1.0 / (2.0*stepsize)
        nSamples = len(self._results[i]._lineProfileGV)
        #self._results[i]._MTFpos = numpy.linspace(start=0, stop=2*self._results[i]._fnyq, num=nSamples, endpoint=False, dtype=numpy.float64)

        # Make the profile positions symmetric:
        self._results[i]._lineProfilePos -= self._profileLength / 2

        # Calculate the line spread function and MTF.
        self._results[i]._lineProfileSmoothed, self._results[i]._LSF, self._results[i]._LSF_windowed, self._results[i]._MTF, self._results[i]._MTFpos = MTF(positions=self._results[i]._lineProfilePos, ESF=self._results[i]._lineProfileGV)

        # Calculate the ideal and the measured MTF10 and MTF20 frequency:
        #self._results[i]._MTF10freq_ideal    = MTFfreq(MTFpos=self._results[i]._MTFpos, MTF=self._results[i]._ideal_MTF_ConvolutionFFT, modulation=0.1)
        #self._results[i]._MTF10freq_measured = MTFfreq(MTFpos=self._results[i]._MTFpos, MTF=self._results[i]._MTF, modulation=0.1)
        self._results[i]._MTF20freq_ideal    = MTFfreq(MTFpos=self._results[i]._MTFpos, MTF=self._results[i]._ideal_MTF_ConvolutionFFT, modulation=0.2)
        self._results[i]._MTF20freq_measured = MTFfreq(MTFpos=self._results[i]._MTFpos, MTF=self._results[i]._MTF, modulation=0.2)


        # Fit a Gaussian to central half of the data:
        sigma = self._results[i]._nominalGaussianSigmaPX
        nBins = len(self._results[i]._lineProfilePos)
        fitParametersInitial = (0, sigma, 1)
        fitPositions = self._results[i]._lineProfilePos[int(nBins//4):int(3*nBins//4)]
        fitData = self._results[i]._LSF[int(nBins//4):int(3*nBins//4)]
        popt, pcov = optimize.curve_fit(f=gaussian, xdata=fitPositions, ydata=fitData, p0=fitParametersInitial)

        perr = numpy.sqrt(numpy.diag(pcov))

        self._results[i]._measuredGaussianSigmaPX      = popt[1]
        self._results[i]._measuredGaussianSigmaPXerror = perr[1]
        self._results[i]._measuredGaussianMuPX         = popt[0]
        self._results[i]._measuredGaussianMuPXerror    = perr[0]
        self._results[i]._measuredGaussianAmpPX        = popt[2]
        self._results[i]._measuredGaussianAmpPXerror   = perr[2]

        self._results[i]._measuredGaussianSigmaMM      = self._results[i]._pixelSize * self._results[i]._measuredGaussianSigmaPX
        self._results[i]._measuredGaussianSigmaMMerror = self._results[i]._pixelSize * self._results[i]._measuredGaussianSigmaPXerror
        self._results[i]._measuredGaussianMuMM         = self._results[i]._pixelSize * self._results[i]._measuredGaussianMuPX
        self._results[i]._measuredGaussianMuMMerror    = self._results[i]._pixelSize * self._results[i]._measuredGaussianMuPXerror
        self._results[i]._measuredGaussianAmpMM        = self._results[i]._pixelSize * self._results[i]._measuredGaussianAmpPX
        self._results[i]._measuredGaussianAmpMMerror   = self._results[i]._pixelSize * self._results[i]._measuredGaussianAmpPXerror

        print("- Measured LSF Fit Mu:    {} +- {}".format(popt[0], perr[0]))
        print("- Measured LSF Fit Sigma: {} +- {}".format(popt[1], perr[1]))
        print("- Measured LSF Fit A:     {} +- {}".format(popt[2], perr[2]))

        # Rescale LSF to maximum of fit function:
        if self._results[i]._measuredGaussianAmpPX != 0.0:
            self._results[i]._LSF /= self._results[i]._measuredGaussianAmpPX
            self._results[i]._LSF_windowed /= self._results[i]._measuredGaussianAmpPX
            self._results[i]._measuredGaussianAmpPX = 1.0

        self._results[i]._LSF_gaussian_fit = numpy.zeros_like(a=self._results[i]._lineProfilePos, dtype=numpy.dtype('float64'))
        for j in range(len(self._results[i]._LSF_gaussian_fit)):
            self._results[i]._LSF_gaussian_fit[j] = gaussian(x=self._results[i]._lineProfilePos[j], mu=self._results[i]._measuredGaussianMuPX, sigma=self._results[i]._measuredGaussianSigmaPX, A=self._results[i]._measuredGaussianAmpPX)
        
        self.writeResultFile(subname=subtestName, results=self._results[i])
        self.writeSummaryFile(subname=subtestName, results=self._results[i])

        log("Evaluation data for test {name}, {subname} written to {dir}.".format(name=self._name, subname=subtestName, dir=self._resultFileDirectory))

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

            fig, (ax1, ax2, ax3) = matplotlib.pyplot.subplots(nrows=3, ncols=1, figsize=(6, 8))
            
            # ESF:
            ax1.plot(self._results[i]._lineProfilePos, self._results[i]._lineProfileGV, linewidth=1.5, label="Measured", color='#ffaa00')
            ax1.plot(self._results[i]._lineProfilePos, self._results[i]._ideal_ESF_Convolution, linewidth=1.0, label="Analytical", color='#000000', linestyle='dotted')

            ax1.set_xlabel("Horizontal distance in px")
            ax1.set_ylabel("ESF")
            ax1.set_xlim([-3*self._results[i]._nominalGaussianSigmaPX, 3*self._results[i]._nominalGaussianSigmaPX])
            ax1.set_title("Edge Spread Function (ESF)")
            #ax1.xaxis.set_ticklabels([])
            ax1.grid(b=True, which='major', axis='both', color='#d9d9d9', linestyle='dashed')
            ax1.grid(b=True, which='minor', axis='both', color='#e7e7e7', linestyle='dotted')
            ax1.legend(loc='best')


            # LSF
            ax2.plot(self._results[i]._lineProfilePos, self._results[i]._LSF, linewidth=1.5, label="Measured", color='#ffaa00')
            ax2.plot(self._results[i]._lineProfilePos, self._results[i]._LSF_gaussian_fit, linewidth=1.0, label="Fit", color='#1f77b4', linestyle='dashed')
            ax2.plot(self._results[i]._lineProfilePos, self._results[i]._ideal_LSF_Convolution, linewidth=1.0, label="Analytical", color='#000000', linestyle='dotted')

            ax2.set_xlabel("Horizontal distance in px")
            ax2.set_ylabel("LSF")
            ax2.set_xlim([-3*self._results[i]._nominalGaussianSigmaPX, 3*self._results[i]._nominalGaussianSigmaPX])
            ax2.set_title("Line Spread Function (LSF)")
            #ax2.xaxis.set_ticklabels([])
            ax2.grid(b=True, which='major', axis='both', color='#d9d9d9', linestyle='dashed')
            ax2.grid(b=True, which='minor', axis='both', color='#e7e7e7', linestyle='dotted')
            ax2.legend(loc='best')


            # Lines for MTF 10:
            #mtf10ideal_vert_x    = numpy.array([self._results[i]._MTF10freq_ideal, self._results[i]._MTF10freq_ideal])
            #mtf10measured_vert_x = numpy.array([self._results[i]._MTF10freq_measured, self._results[i]._MTF10freq_measured])
            #mtf10_vert_y = numpy.array([-0.1, 0.1])

            #mtf10_horz_x = numpy.array([0, 3.0*max(self._results[i]._MTF10freq_measured, self._results[i]._MTF10freq_ideal)])
            #mtf10_horz_y = numpy.array([0.1, 0.1])


            # Lines for MTF 20:
            mtf20ideal_vert_x    = numpy.array([self._results[i]._MTF20freq_ideal, self._results[i]._MTF20freq_ideal])
            mtf20measured_vert_x = numpy.array([self._results[i]._MTF20freq_measured, self._results[i]._MTF20freq_measured])
            mtf20_vert_y = numpy.array([-0.1, 0.2])

            mtf20_horz_x = numpy.array([0, 3.0*max(self._results[i]._MTF20freq_measured, self._results[i]._MTF20freq_ideal)])
            mtf20_horz_y = numpy.array([0.2, 0.2])

            # MTF
            ax3.plot(self._results[i]._MTFpos, self._results[i]._MTF, linewidth=1.5, label="Measured", color='#ffaa00')
            ax3.plot(self._results[i]._MTFpos, self._results[i]._ideal_MTF_ConvolutionFFT, linewidth=1.0, label="Analytical", color='#000000', linestyle='dotted')

            #ax3.plot(mtf10_horz_x, mtf10_horz_y,         linewidth=0.5, color='#000000', label="MTF10")
            #ax3.plot(mtf10ideal_vert_x, mtf10_vert_y,    linewidth=0.5, color='#000000', label="")
            #ax3.plot(mtf10measured_vert_x, mtf10_vert_y, linewidth=0.5, color='#000000', label="")

            #ax3.plot(mtf20_horz_x, mtf20_horz_y,         linewidth=0.5, color='#808080', label="MTF20")
            ax3.plot(mtf20ideal_vert_x, mtf20_vert_y,    linewidth=0.5, color='#808080', label="", linestyle='dotted')
            ax3.plot(mtf20measured_vert_x, mtf20_vert_y, linewidth=0.5, color='#808080', label="")

            ax3.set_xlabel("Modulation frequency in cycles/px")
            ax3.set_ylabel("Modulation contrast")
            ax3.set_xlim([0, 3.0*max(self._results[i]._MTF20freq_measured, self._results[i]._MTF20freq_ideal)])
            ax3.set_ylim([-0.1, 1.1])
            ax3.set_yticks(numpy.array([0, 0.2, 0.4, 0.6, 0.8, 1]))
            ax3.set_title("Modulation Transfer Function (MTF)")
            #ax3.xaxis.set_ticklabels([])
            ax3.grid(b=True, which='major', axis='both', color='#d9d9d9', linestyle='dashed')
            ax3.grid(b=True, which='minor', axis='both', color='#e7e7e7', linestyle='dotted')
            ax3.legend(loc='best')

            fig.tight_layout(pad=2.5)

            plotFilename = "{dir}/{name}_{subname}_ESF_LSF_MTF.png".format(dir=self._resultFileDirectory, name=self._name, subname=subtestName)
            matplotlib.pyplot.savefig(plotFilename)
            fig.clf()
            matplotlib.pyplot.close('all')
        except:
            log("Warning: Error plotting results for test {name}, {subname} using matplotlib.".format(name=self._name, subname=subtestName))
