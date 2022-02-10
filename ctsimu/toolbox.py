# -*- coding: UTF-8 -*-
import os
import json

from .helpers import *
from .geometry import Geometry
from .image import ImageStack
from .processing.pipeline import Pipeline
from .processing.flat_field import Step_FlatFieldCorrection
from .ctsimu_evaluations.test2D_FB_2 import Test2D_FB_2
from .ctsimu_evaluations.test2D_FB_1 import Test2D_FB_1
from .ctsimu_evaluations.test2D_SW_1 import Test2D_SW_1
from .ctsimu_evaluations.test2D_SW_2 import Test2D_SW_2
from .ctsimu_evaluations.test2D_SW_3 import Test2D_SW_3
from .ctsimu_evaluations.test2D_SW_4 import Test2D_SW_4
from .ctsimu_evaluations.test2D_DW_1 import Test2D_DW_1
from .ctsimu_evaluations.test2D_WE_1 import Test2D_WE_1
from .ctsimu_evaluations.test2D_WE_2 import Test2D_WE_2
from .ctsimu_evaluations.test2D_HS_1 import Test2D_HS_1

class Toolbox:
    """ Manages a test run, including preliminary flat field corrections, based on metadata JSON input files. """

    def __init__(self, op, *metadata, **subtests):
        if op == "info":
            for metafile in metadata:
                geo = Geometry(jsonFile=metafile)
                print(geo.info())
                
            return

        ffRescaleFactor = float(60000.0)
        ffRescaleOffset = 0
        if op == "correction":
            for key, metafile in subtests.items():
                if key == "rescaleFactor":
                    ffRescaleFactor = metafile
                elif key == "offsetAfterRescale":
                    ffRescaleOffset = metafile

        # Prepare a pipeline
        pipeline = Pipeline()
        ffCorrection = Step_FlatFieldCorrection()
        ffCorrection.setFlatFieldRescaleFactor(ffRescaleFactor, ffRescaleOffset)
        pipeline.addStep(ffCorrection)

        evaluationStep = None
        if op=="2D-FB-2":
            evaluationStep = Test2D_FB_2()
        elif op=="2D-FB-1":
            evaluationStep = Test2D_FB_1()
        elif op=="2D-SW-1":
            evaluationStep = Test2D_SW_1()
        elif op=="2D-SW-2":
            evaluationStep = Test2D_SW_2()
        elif op=="2D-SW-3":
            evaluationStep = Test2D_SW_3()
        elif op=="2D-SW-4":
            evaluationStep = Test2D_SW_4()
        elif op=="2D-DW-1":
            evaluationStep = Test2D_DW_1()
        elif op=="2D-WE-1":
            evaluationStep = Test2D_WE_1()
        elif op=="2D-WE-2":
            evaluationStep = Test2D_WE_2()
        elif op=="2D-HS-1":
            evaluationStep = Test2D_HS_1()
        elif op=="correction":
            pass
        else:
            log("Not a valid test scenario name or operation: {op}".format(op=op))

        # Evaluate keyword list and see if metadata files are provided
        # together with keywords that might identify subtest scenarios.
        metadataList = []
        if op != "correction":
            for key, metafile in subtests.items():
                if evaluationStep is not None:
                    evaluationStep.addSubtest(key)
                    metadataList.append(metafile)
                else:
                    raise Exception("'{op}' is not a valid command for the toolbox.".format(op=op))

        if len(metadataList) == 0:  # Fallback to non-keyword list. Works only for tests that require no subtest identification.
            metadataList = list(metadata)

        if evaluationStep != None:
            pipeline.addStep(evaluationStep)

        for metafile in metadataList:
            if os.path.isfile(metafile):
                metafileAbsolute = ""
                metafileAbsDir   = ""

                # Try to find the working directory with the metadata files:
                if metafile != None:
                    metafileAbsolute = os.path.abspath(metafile)
                    metafileAbsDir   = os.path.dirname(metafileAbsolute)

                resultFileDir = "{testname}-results".format(testname=op)
                if metafileAbsDir != None and len(metafileAbsDir) > 0:
                    resultFileDir = metafileAbsDir + "/" + resultFileDir

                if evaluationStep != None:
                    pipeline.step(1).setResultFileDirectory(resultFileDir)

                log("---")
                log("Metadata File: {mf}".format(mf=metafile))

                jsonFile = open(metafile, "r")
                jsonText = jsonFile.read()
                jsonFile.close()

                jsonDict = json.loads(jsonText)

                metafileAbsolute = os.path.abspath(metafile)
                metafileAbsDir   = os.path.dirname(metafileAbsolute)

                projFilename = getFieldOrNone(jsonDict, "output", "projections", "filename")
                if projFilename != None:
                    if not os.path.isabs(projFilename): # Check if an absolute path is provided
                        # If a relative path is provided, this path is relative to the location of the metadata file:
                        projFilename = "{metafileDir}/{imgfile}".format(metafileDir=metafileAbsDir, imgfile=projFilename)

                log("  Projection File(s):    {proj}".format(proj=projFilename))

                projNumber = getFieldOrNone(jsonDict, "output", "projections", "number")
                log("  Number of Projections: {nproj}".format(nproj=projNumber))

                projDataType = getFieldOrNone(jsonDict, "output", "projections", "datatype")
                log("  Data Type:             {dataType}".format(dataType=projDataType))

                projByteOrder = getFieldOrNone(jsonDict, "output", "projections", "byteorder")
                log("  Byte Order:            {byteOrder}".format(byteOrder=projByteOrder))

                projHeaderSizeFile = getFieldOrNone(jsonDict, "output", "projections", "headersize", "file")
                log("  File Header Size:      {fileHeaderSize}".format(fileHeaderSize=projHeaderSizeFile))

                projHeaderSizeImage = getFieldOrNone(jsonDict, "output", "projections", "headersize", "image")
                log("  Image Header Size:     {imageHeaderSize}".format(imageHeaderSize=projHeaderSizeImage))

                width  = getFieldOrNone(jsonDict, "output", "projections", "dimensions", "x", "value")
                log("  Projection Width:      {w}".format(w=width))

                height = getFieldOrNone(jsonDict, "output", "projections", "dimensions", "y", "value")
                log("  Projection Height:     {h}".format(h=height))

                darkFilename  = getFieldOrNone(jsonDict, "output", "projections", "dark_field", "filename")
                if darkFilename != None:
                    if not os.path.isabs(darkFilename): # Check if an absolute path is provided
                        # If a relative path is provided, this path is relative to the location of the metadata file:
                        darkFilename = "{metafileDir}/{imgfile}".format(metafileDir=metafileAbsDir, imgfile=darkFilename)

                log("  Dark Field(s):         {df}".format(df=darkFilename))

                darkNumber    = getFieldOrNone(jsonDict, "output", "projections", "dark_field", "number")
                log("  Number of Dark Fields: {n}".format(n=darkNumber))

                darkCorrected = getFieldOrNone(jsonDict, "output", "projections", "dark_field", "projections_corrected")
                log("  DF correction applied: {dfapplied}".format(dfapplied=darkCorrected))

                flatFilename  = getFieldOrNone(jsonDict, "output", "projections", "flat_field", "filename")
                if flatFilename != None:
                    if not os.path.isabs(flatFilename): # Check if an absolute path is provided
                        # If a relative path is provided, this path is relative to the location of the metadata file:
                        flatFilename = "{metafileDir}/{imgfile}".format(metafileDir=metafileAbsDir, imgfile=flatFilename)

                log("  Flat Field(s):         {ff}".format(ff=flatFilename))

                flatNumber    = getFieldOrNone(jsonDict, "output", "projections", "flat_field", "number")
                log("  Number of Flat Fields: {n}".format(n=flatNumber))

                flatCorrected = getFieldOrNone(jsonDict, "output", "projections", "flat_field", "projections_corrected")
                log("  FF correction applied: {ffapplied}".format(ffapplied=flatCorrected))

                log("---")

                if projFilename is None:
                    raise Exception("No projection filename provided in metadata file: {filename}".format(filename=metafile))

                projections = ImageStack(
                    filePattern=projFilename,
                    width=width,
                    height=height,
                    dataType=projDataType,
                    byteOrder=projByteOrder,
                    rawFileHeaderSize=projHeaderSizeFile,
                    rawImageHeaderSize=projHeaderSizeImage,
                    slices=projNumber,
                    flipByteOrder=False
                    )

                pipeline.setInputFileStack(projections)
                pipeline.setOutputFileStack(None)  # set later, if corrections will be applied
                pipeline.step(0).setOffsetCorrection(False)
                pipeline.step(0).setGainCorrection(False)

                correctedFilename = os.path.dirname(projFilename) + "/corrected/" + os.path.basename(projFilename)
                outputFiles = ImageStack(
                    filePattern=correctedFilename,
                    width=width,
                    height=height,
                    dataType=projDataType,
                    byteOrder=projByteOrder,
                    rawFileHeaderSize=projHeaderSizeFile,
                    rawImageHeaderSize=projHeaderSizeImage,
                    slices=projNumber,
                    flipByteOrder=False
                    )
                #log("Corrected projections are stored in: {correctedPath}".format(correctedPath=correctedFilename))

                if not darkCorrected:
                    if darkNumber > 0:
                        if darkFilename != None:
                            darkStack = ImageStack(
                                            filePattern=darkFilename,
                                            width=width,
                                            height=height,
                                            dataType=projDataType,
                                            byteOrder=projByteOrder,
                                            rawFileHeaderSize=projHeaderSizeFile,
                                            rawImageHeaderSize=projHeaderSizeImage,
                                            slices=darkNumber,
                                            flipByteOrder=False
                                        )
                            pipeline.step(0).setDarkFileStack(darkStack)
                            pipeline.setOutputFileStack(outputFiles)

                if not flatCorrected:
                    if flatNumber > 0:
                        if flatFilename != None:
                            flatStack = ImageStack(
                                            filePattern=flatFilename,
                                            width=width,
                                            height=height,
                                            dataType=projDataType,
                                            byteOrder=projByteOrder,
                                            rawFileHeaderSize=projHeaderSizeFile,
                                            rawImageHeaderSize=projHeaderSizeImage,
                                            slices=flatNumber,
                                            flipByteOrder=False
                                        )
                            pipeline.step(0).setFlatFileStack(flatStack)
                            pipeline.setOutputFileStack(outputFiles)

                pipeline.run()

            else:
                raise Exception("Cannot access metadata file: {filename}".format(filename=metafile))

        if evaluationStep != None:
            # Run the follow-up procedure to output the results.
            evaluationStep.followUp()