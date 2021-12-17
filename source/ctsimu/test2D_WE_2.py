# -*- coding: UTF-8 -*-
from .test import *
from .general import *
from .geoprimitives import Vector, Polygon

class Test2D_WE_2(generalTest):
    """ CTSimU test 2D-WE-2: Effect of partial pixel coverage. """

    def __init__(self, resultFileDirectory=".", name=None, rawOutput=False):

        generalTest.__init__(
            self,
            testName="2D-WE-2",
            name=name,
            nExpectedRuns=1,
            resultFileDirectory=resultFileDirectory,
            rawOutput=rawOutput)
        
        self._geometry = None
        self._analyticalIntensityProfileImage = None  # analytical flat field
        self._analyticalEdgeImage             = None  # stores the analytically computed edge image
        self._analyticalEdgeImageFF           = None  # stores the analytically computed edge image, flat-field corrected

        self._differenceImage                 = None  # difference between analytical edge image and provided image

        self._numberOfPartlyCoveredPixels_analytical = 0
        self._numberOfPartlyCoveredPixels_measured   = 0
        self._matchingPixels                         = 0  # number of partly covered pixels with matching positions
        self._partlyCovered_rmsd                     = 0  # RMSD of pixels in difference image that are partly covered in analytical image

        # Prepare the clipping rectangle for the analytical
        # calculation of the ideal edge image. In pixel coordinates.
        A = Vector(   0,    0, 0)
        B = Vector(   0,  300, 0)
        C = Vector(-300,  300, 0)
        D = Vector(-300,    0, 0)

        edgeAngle = 3 * (math.pi/180.0)  # 3 deg edge rotation
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

        if not self._prepared:
            self._jsonScenarioFile = "scenarios/2D-WE-2_2021-07-07v02r00fbo.json"

            # Calculate the analytical image of an edge covering the image:
            if(self._jsonScenarioFile != None):
                self._geometry = Geometry(jsonFileFromPkg=self._jsonScenarioFile)
                self._analyticalIntensityProfileImage, self._analyticalEdgeImage = self._geometry.createDetectorFlatField_sphere(self._clippingRectangle)
                self._analyticalEdgeImageFF = copy.deepcopy(self._analyticalEdgeImage)
                self._analyticalEdgeImageFF.applyFlatfield(ref=self._analyticalIntensityProfileImage, rescaleFactor=60000.0)

                # Raise analytical images to maximum grey value of 60000 before saving them.
                # This rescaling does not affect the previous FF correction.
                self._analyticalIntensityProfileImage.renormalize(newMin=0.0, newMax=60000.0, currentMin=0.0, currentMax=1.0)
                self._analyticalEdgeImage.renormalize(newMin=0.0, newMax=60000.0, currentMin=0.0, currentMax=1.0)


                # Write analytical images:
                if self._rawOutput:
                    self._analyticalIntensityProfileImage.saveRAW("{dir}/{name}_analytical_flat.raw".format(dir=self._resultFileDirectory, name=self._name), dataType="float32", addInfo=True)
                    self._analyticalEdgeImage.saveRAW("{dir}/{name}_analytical_edge.raw".format(dir=self._resultFileDirectory, name=self._name), dataType="float32", addInfo=True)
                    self._analyticalEdgeImageFF.saveRAW("{dir}/{name}_analytical_edge_corrected.raw".format(dir=self._resultFileDirectory, name=self._name), dataType="float32", addInfo=True)
                else: # TIFF
                    self._analyticalIntensityProfileImage.save("{dir}/{name}_analytical_flat.tif".format(dir=self._resultFileDirectory, name=self._name), dataType="float32")
                    self._analyticalEdgeImage.save("{dir}/{name}_analytical_edge.tif".format(dir=self._resultFileDirectory, name=self._name), dataType="float32")
                    self._analyticalEdgeImageFF.save("{dir}/{name}_analytical_edge_corrected.tif".format(dir=self._resultFileDirectory, name=self._name), dataType="float32")

                self._numberOfPartlyCoveredPixels_analytical = numpy.count_nonzero((self._analyticalEdgeImageFF._px > 0) & (self._analyticalEdgeImageFF._px < 60000))

                print("Number of partly covered pixels in analytical image: {}".format(self._numberOfPartlyCoveredPixels_analytical))

                self._prepared = True
            else:
                raise Exception("Test 2D-WE-2: Please provide a JSON scenario description.")


    def run(self, image):
        self.prepare()
        self._currentRun += 1

        self._numberOfPartlyCoveredPixels_measured = numpy.count_nonzero((image._px > 0) & (image._px < 60000))

        print("Number of partly covered pixels in measured image:   {}".format(self._numberOfPartlyCoveredPixels_measured))

        # Match in covered pixels:
        self._matchingPixels = numpy.count_nonzero((self._analyticalEdgeImageFF._px > 0) & (self._analyticalEdgeImageFF._px < 60000) & (image._px > 0) & (image._px < 60000))

        print("Number of pixels partly covered in both images:      {}".format(self._matchingPixels))

        self._matchRatio = self._matchingPixels / self._numberOfPartlyCoveredPixels_analytical

        print("Ratio (covered in both)/(covered in analytical):     {:.2f}%".format(100.0*self._matchRatio))


        self._differenceImage = copy.deepcopy(image)
        self._differenceImage.subtractImage(self._analyticalEdgeImageFF)    

        if self._rawOutput:
            self._differenceImage.saveRAW("{dir}/{name}_difference.raw".format(dir=self._resultFileDirectory, name=self._name), dataType="float32", addInfo=True)
        else: # TIFF
            self._differenceImage.save("{dir}/{name}_difference.tif".format(dir=self._resultFileDirectory, name=self._name), dataType="float32")


        self._differenceImage.square()
        self._partlyCovered_rmsd = math.sqrt(numpy.mean(self._differenceImage._px[numpy.nonzero((self._analyticalEdgeImageFF._px > 0) & (self._analyticalEdgeImageFF._px < 60000))]))

        print("RMSD of partly covered (analytical) pixels:          {:.2f} GV".format(self._partlyCovered_rmsd))



        # Write a CSV file with all pixel values and their differences
        csvText = "# x [px]\ty [px]\tAnalytical GV\tMeasured GV\tDifference\n"

        partlyCoveredCoordinates_analytical = numpy.nonzero((self._analyticalEdgeImageFF._px > 0) & (self._analyticalEdgeImageFF._px < 60000))

        for i in range(len(partlyCoveredCoordinates_analytical[0])):
            x = partlyCoveredCoordinates_analytical[1][i]
            y = partlyCoveredCoordinates_analytical[0][i]
            analytical = self._analyticalEdgeImageFF._px[y][x]
            measured   = image._px[y][x]
            difference = measured - analytical

            csvText += "{x}\t{y}\t{analytical:.3f}\t{measured:.3f}\t{delta:.3f}\n".format(x=x, y=y, analytical=analytical, measured=measured, delta=difference)

        csvFileName = "{dir}/{name}_pixel_list.txt".format(dir=self._resultFileDirectory, name=self._name)
        with open(csvFileName, 'w') as csvFile:
            csvFile.write(csvText)
            csvFile.close()

        self.plotResults()

        return image

    def plotResults(self):
        pass

    def followUp(self):
        log("Writing evaluation results...")

        summaryText  = "# |A| = Number of partly covered pixels in analytical image: {}\n".format(self._numberOfPartlyCoveredPixels_analytical)
        summaryText += "# |M| = Number of partly covered pixels in measured image:   {}\n".format(self._numberOfPartlyCoveredPixels_measured)
        summaryText += "# Number of pixels partly covered in both images:            {}\n".format(self._matchingPixels)
        summaryText += "# r = Ratio (covered in both)/(covered in analytical):       {:.2f}%\n".format(100.0*self._matchRatio)
        summaryText += "# RMSD [partly covered analytical vs. measured]:             {:.2f} GV".format(self._partlyCovered_rmsd)

        summaryFileName = "{dir}/{name}_summary.txt".format(dir=self._resultFileDirectory, name=self._name)
        with open(summaryFileName, 'w') as summaryFile:
            summaryFile.write(summaryText)
            summaryFile.close()
