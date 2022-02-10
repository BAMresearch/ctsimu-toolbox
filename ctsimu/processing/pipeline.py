import os    # File and path handling
import numpy
from ..image import *
from ..helpers import *

from .step import Step

class Pipeline:
    """ Perform given processing steps on input images, save as output. """

    def __init__(self, inputFileStack=None, outputFileStack=None):
        self.inputFileStack   = None
        self.outputFileStack  = None

        self.setInputFileStack(inputFileStack)
        self.setOutputFileStack(outputFileStack)

        self.processingSteps = []

    def setInputFileStack(self, inputFileStack):
        self.inputFileStack = createImageStack(inputFileStack)

    def setOutputFileStack(self, outputFileStack):
        self.outputFileStack = createImageStack(outputFileStack)

    def addStep(self, step):
        """ Add a processing step to the pipeline. """
        if isinstance(step, Step):
            step.setPipeline(self)
            self.processingSteps.append(step)
        else:
            raise Exception("Failed to add processing step. Not a valid type of processing step: {}".format(step))

    def step(self, i):
        if i < len(self.processingSteps):
            return self.processingSteps[i]
        else:
            raise Exception("Processing Pipeline: Requested processing step id (#{i}) exceeds number of stored processing steps ({n}).".format(i=i, n=len(self.processingSteps)))

    def run(self):
        """ Run the processing pipeline with its current configuration. """
        if self.inputFileStack is not None:
            nImagesProcessed = 0

            # Build stack if it hasn't been built yet:
            if not self.inputFileStack.built:
                self.inputFileStack.buildStack()

            saveOutputFiles = False
            if isinstance(self.outputFileStack, ImageStack):
                saveOutputFiles = True

                # Set output parameters to input parameters if no others given:
                if self.outputFileStack.getFileDataType() is None:
                    self.outputFileStack.setFileDataType(self.inputFileStack.getFileDataType())

                if self.outputFileStack.getFileByteOrder() is None:
                    self.outputFileStack.setFileByteOrder(self.inputFileStack.getFileByteOrder())

            # Treat projection files
            if self.inputFileStack.nSlices > 0:
                outputFiles = None               

                if saveOutputFiles:
                    outputFiles = self.outputFileStack.files
                    generalOutputName = self.outputFileStack.getFilename()
                    outputFolder   = os.path.dirname(generalOutputName)
                    outputBasename = os.path.basename(generalOutputName)

                    if outputFolder == "" or outputFolder is None:
                        outputFolder = "."

                    if '%' in outputBasename:
                        self.outputFileStack.setVolumeChunk(False)
                        leadOutputName, trailOutputName, nDigitsExpected = self.outputFileStack.fileStackInfo(outputBasename)
                    else:
                        self.outputFileStack.setVolumeChunk(True)

                print("Running image processing pipeline.")
                for i in range(self.inputFileStack.nSlices):
                    progress = 100*(float(i+1)/float(self.inputFileStack.nSlices))

                    print("\rImage {}/{} ({:0.1f}%)  \r".format((i+1), self.inputFileStack.nSlices, progress), end='')
                    image = self.inputFileStack.getImage(index=i, outputFile=outputFiles)

                    # Run through processing steps:
                    for step in self.processingSteps:
                        image = step.run(image)
                        if image is None:
                            raise Exception("Step {i} did not return a valid image from its run() method.".format(i=step.getIdentifier()))

                    if saveOutputFiles:
                        # Append if output target is a volume chunk:
                        appendMode = False
                        if self.outputFileStack.isVolumeChunk():
                            outputPath = generalOutputName
                            if (i>0):
                                appendMode = True
                        else:
                            fileNameDigit = i
                            if len(self.inputFileStack.fileNumbers) > 0:
                                if len(self.inputFileStack.fileNumbers) > i:
                                    fileNameDigit = self.inputFileStack.fileNumbers[i]

                            outputPath = "{folder}/{lead}{digits:{fill}{nDigits}}{trail}".format(folder=outputFolder, lead=leadOutputName, digits=fileNameDigit, fill='0', nDigits=nDigitsExpected, trail=trailOutputName)
               
                        image.save(filename=outputPath, appendChunk=appendMode)

                    nImagesProcessed += 1

                print("Image {}/{} (100%)  ".format((i+1), self.inputFileStack.nSlices))

            else:
                log("No projection file(s) found that match the given name pattern.")

            if(nImagesProcessed == 0):
                log("No image files processed.      ")
            elif(nImagesProcessed == 1):
                log("1 image file processed.      ")
            else:
                log("{} image files processed.      ".format(nImagesProcessed))

            #if not saveOutputFiles:
            #    log("No output images were saved.")
        else:
            log("Nothing to do. The pipeline has no input files.")