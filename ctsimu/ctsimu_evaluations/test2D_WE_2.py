# -*- coding: UTF-8 -*-
from ..test import *
from ..helpers import *
from ..primitives import Vector, Polygon

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
        
        self.geometry = None
        self.analyticalIntensityProfileImage = None  # analytical flat field
        self.analyticalEdgeImage             = None  # stores the analytically computed edge image
        self.analyticalEdgeImageFF           = None  # stores the analytically computed edge image, flat-field corrected

        self.differenceImage                 = None  # difference between analytical edge image and provided image

        self.numberOfPartlyCoveredPixels_analytical = 0
        self.numberOfPartlyCoveredPixels_measured   = 0
        self.matchingPixels                         = 0  # number of partly covered pixels with matching positions
        self.partlyCovered_rmsd                     = 0  # RMSD of pixels in difference image that are partly covered in analytical image

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

        self.clippingRectangle = Polygon(A, B, C, D)


    def prepare(self):
        """ Preparations before the test will be run with the images from the pipeline. """
        if not isinstance(self.pipe, Pipeline):
            self.prepared = False
            raise Exception("Step must be part of a processing pipeline before it can prepare. Current pipeline: {}".format(self.pipe))

        if not self.prepared:
            self.jsonScenarioFile = "ctsimu_evaluations/scenarios/2D-WE-2_2021-07-07v02r00fbo.json"

            # Calculate the analytical image of an edge covering the image:
            if(self.jsonScenarioFile != None):
                self.geometry = Geometry(jsonFileFromPkg=self.jsonScenarioFile)
                self.analyticalIntensityProfileImage, self.analyticalEdgeImage = self.geometry.createDetectorFlatField_sphere(self.clippingRectangle)
                self.analyticalEdgeImageFF = copy.deepcopy(self.analyticalEdgeImage)
                self.analyticalEdgeImageFF.applyFlatfield(ref=self.analyticalIntensityProfileImage, rescaleFactor=60000.0)

                # Raise analytical images to maximum grey value of 60000 before saving them.
                # This rescaling does not affect the previous FF correction.
                self.analyticalIntensityProfileImage.renormalize(newMin=0.0, newMax=60000.0, currentMin=0.0, currentMax=1.0)
                self.analyticalEdgeImage.renormalize(newMin=0.0, newMax=60000.0, currentMin=0.0, currentMax=1.0)


                # Write analytical images:
                if self.rawOutput:
                    self.analyticalIntensityProfileImage.saveRAW("{dir}/{name}_analytical_flat.raw".format(dir=self.resultFileDirectory, name=self.name), dataType="float32", addInfo=True)
                    self.analyticalEdgeImage.saveRAW("{dir}/{name}_analytical_edge.raw".format(dir=self.resultFileDirectory, name=self.name), dataType="float32", addInfo=True)
                    self.analyticalEdgeImageFF.saveRAW("{dir}/{name}_analytical_edge_corrected.raw".format(dir=self.resultFileDirectory, name=self.name), dataType="float32", addInfo=True)
                else: # TIFF
                    self.analyticalIntensityProfileImage.save("{dir}/{name}_analytical_flat.tif".format(dir=self.resultFileDirectory, name=self.name), dataType="float32")
                    self.analyticalEdgeImage.save("{dir}/{name}_analytical_edge.tif".format(dir=self.resultFileDirectory, name=self.name), dataType="float32")
                    self.analyticalEdgeImageFF.save("{dir}/{name}_analytical_edge_corrected.tif".format(dir=self.resultFileDirectory, name=self.name), dataType="float32")

                self.numberOfPartlyCoveredPixels_analytical = numpy.count_nonzero((self.analyticalEdgeImageFF.px > 0) & (self.analyticalEdgeImageFF.px < 60000))

                print("Number of partly covered pixels in analytical image: {}".format(self.numberOfPartlyCoveredPixels_analytical))

                self.prepared = True
            else:
                raise Exception("Test 2D-WE-2: Please provide a JSON scenario description.")


    def run(self, image):
        self.prepare()
        self.currentRun += 1

        self.numberOfPartlyCoveredPixels_measured = numpy.count_nonzero((image.px > 0) & (image.px < 60000))

        print("Number of partly covered pixels in measured image:   {}".format(self.numberOfPartlyCoveredPixels_measured))

        # Match in covered pixels:
        self.matchingPixels = numpy.count_nonzero((self.analyticalEdgeImageFF.px > 0) & (self.analyticalEdgeImageFF.px < 60000) & (image.px > 0) & (image.px < 60000))

        print("Number of pixels partly covered in both images:      {}".format(self.matchingPixels))

        self.matchRatio = self.matchingPixels / self.numberOfPartlyCoveredPixels_analytical

        print("Ratio (covered in both)/(covered in analytical):     {:.2f}%".format(100.0*self.matchRatio))


        self.differenceImage = copy.deepcopy(image)
        self.differenceImage.subtractImage(self.analyticalEdgeImageFF)    

        if self.rawOutput:
            self.differenceImage.saveRAW("{dir}/{name}_difference.raw".format(dir=self.resultFileDirectory, name=self.name), dataType="float32", addInfo=True)
        else: # TIFF
            self.differenceImage.save("{dir}/{name}_difference.tif".format(dir=self.resultFileDirectory, name=self.name), dataType="float32")


        self.differenceImage.square()
        self.partlyCovered_rmsd = math.sqrt(numpy.mean(self.differenceImage.px[numpy.nonzero((self.analyticalEdgeImageFF.px > 0) & (self.analyticalEdgeImageFF.px < 60000))]))

        print("RMSD of partly covered (analytical) pixels:          {:.2f} GV".format(self.partlyCovered_rmsd))



        # Write a CSV file with all pixel values and their differences
        csvText = "# x [px]\ty [px]\tAnalytical GV\tMeasured GV\tDifference\n"

        partlyCoveredCoordinates_analytical = numpy.nonzero((self.analyticalEdgeImageFF.px > 0) & (self.analyticalEdgeImageFF.px < 60000))

        for i in range(len(partlyCoveredCoordinates_analytical[0])):
            x = partlyCoveredCoordinates_analytical[1][i]
            y = partlyCoveredCoordinates_analytical[0][i]
            analytical = self.analyticalEdgeImageFF.px[y][x]
            measured   = image.px[y][x]
            difference = measured - analytical

            csvText += "{x}\t{y}\t{analytical:.3f}\t{measured:.3f}\t{delta:.3f}\n".format(x=x, y=y, analytical=analytical, measured=measured, delta=difference)

        csvFileName = "{dir}/{name}_pixel_list.txt".format(dir=self.resultFileDirectory, name=self.name)
        with open(csvFileName, 'w') as csvFile:
            csvFile.write(csvText)
            csvFile.close()

        self.plotResults()

        return image

    def plotResults(self):
        pass

    def followUp(self):
        log("Writing evaluation results...")

        summaryText  = "# |A| = Number of partly covered pixels in analytical image: {}\n".format(self.numberOfPartlyCoveredPixels_analytical)
        summaryText += "# |M| = Number of partly covered pixels in measured image:   {}\n".format(self.numberOfPartlyCoveredPixels_measured)
        summaryText += "# Number of pixels partly covered in both images:            {}\n".format(self.matchingPixels)
        summaryText += "# r = Ratio (covered in both)/(covered in analytical):       {:.2f}%\n".format(100.0*self.matchRatio)
        summaryText += "# RMSD [partly covered analytical vs. measured]:             {:.2f} GV".format(self.partlyCovered_rmsd)

        summaryFileName = "{dir}/{name}_summary.txt".format(dir=self.resultFileDirectory, name=self.name)
        with open(summaryFileName, 'w') as summaryFile:
            summaryFile.write(summaryText)
            summaryFile.close()
