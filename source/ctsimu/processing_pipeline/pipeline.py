import os    # File and path handling
import numpy
from ..image import *
from ..general import *

from .processing_step import ProcessingStep

class ProcessingPipeline:
    """ Perform given processing steps on input images, save as output. """

    def __init__(self, inputFileStack=None, outputFileStack=None):
        self._inputFileStack   = None
        self._outputFileStack  = None

        self.setInputFileStack(inputFileStack)
        self.setOutputFileStack(outputFileStack)

        self._processingSteps = []

    def setInputFileStack(self, inputFileStack):
        self._inputFileStack = createImageStack(inputFileStack)

    def setOutputFileStack(self, outputFileStack):
        self._outputFileStack = createImageStack(outputFileStack)

    def addStep(self, step):
        """ Add a processing step to the pipeline. """
        if isinstance(step, ProcessingStep):
            step.setPipeline(self)
            self._processingSteps.append(step)
        else:
            raise Exception("Failed to add processing step. Not a valid type of processing step: {}".format(step))

    def step(self, i):
        if i < len(self._processingSteps):
            return self._processingSteps[i]
        else:
            raise Exception("Processing Pipeline: Requested processing step id (#{i}) exceeds number of stored processing steps ({n}).".format(i=i, n=len(self._processingSteps)))

    def run(self):
        """ Run the processing pipeline with its current configuration. """
        if self._inputFileStack is not None:
            nImagesProcessed = 0

            # Build stack if it hasn't been built yet:
            if not self._inputFileStack._built:
                self._inputFileStack.buildStack()

            saveOutputFiles = False
            if isinstance(self._outputFileStack, ImageStack):
                saveOutputFiles = True

                # Set output parameters to input parameters if no others given:
                if self._outputFileStack.getFileDataType() == None:
                    self._outputFileStack.setFileDataType(self._inputFileStack.getFileDataType())

                if self._outputFileStack.getFileByteOrder() == None:
                    self._outputFileStack.setFileByteOrder(self._inputFileStack.getFileByteOrder())

            # Treat projection files
            if self._inputFileStack.nSlices() > 0:
                outputFiles = None               

                if saveOutputFiles:
                    outputFiles = self._outputFileStack._files
                    generalOutputName = self._outputFileStack.getFilename()
                    outputFolder   = os.path.dirname(generalOutputName)
                    outputBasename = os.path.basename(generalOutputName)

                    if outputFolder == "" or outputFolder == None:
                        outputFolder = "."

                    if '%' in outputBasename:
                        self._outputFileStack.setVolumeChunk(False)
                        leadOutputName, trailOutputName, nDigitsExpected = self._outputFileStack.fileStackInfo(outputBasename)
                    else:
                        self._outputFileStack.setVolumeChunk(True)

                print("Running image processing pipeline.")
                for i in range(self._inputFileStack.nSlices()):
                    progress = 100*(float(i+1)/float(self._inputFileStack.nSlices()))

                    print("\rImage {}/{} ({:0.1f}%)  \r".format((i+1), self._inputFileStack.nSlices(), progress), end='')
                    image = self._inputFileStack.getImage(index=i, outputFile=outputFiles)

                    # Run through processing steps:
                    for step in self._processingSteps:
                        image = step.run(image)
                        if image == None:
                            raise Exception("Step {i} did not return a valid image from its run() method.".format(i=step.getIdentifier()))

                    if saveOutputFiles:
                        # Append if output target is a volume chunk:
                        appendMode = False
                        if self._outputFileStack.isVolumeChunk():
                            outputPath = generalOutputName
                            if (i>0):
                                appendMode = True
                        else:
                            fileNameDigit = i
                            if len(self._inputFileStack._fileNumbers) > 0:
                                if len(self._inputFileStack._fileNumbers) > i:
                                    fileNameDigit = self._inputFileStack._fileNumbers[i]

                            outputPath = "{folder}/{lead}{digits:{fill}{nDigits}}{trail}".format(folder=outputFolder, lead=leadOutputName, digits=fileNameDigit, fill='0', nDigits=nDigitsExpected, trail=trailOutputName)
               
                        image.save(filename=outputPath, appendChunk=appendMode)

                    nImagesProcessed += 1

                print("Image {}/{} (100%)  ".format((i+1), self._inputFileStack.nSlices()))

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