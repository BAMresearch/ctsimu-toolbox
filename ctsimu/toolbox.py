# -*- coding: UTF-8 -*-
import os
import json

from .helpers import *
from .geometry import Geometry
from .image import ImageStack
from .scenario import Scenario
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
    def __init__(self, op, *metadata_files, **kwargs):
        if op == "info":
            self.info(metadata_files)
        elif op == "correction":
            self.correction(*metadata_files, **kwargs)
        elif op == "post_process_recursive":
            self.post_process_recursive(kwargs)
        elif op.startswith("2D-"):
            # Possibly a 2D test... try it!
            self.test_2D(op, *metadata_files, **kwargs)

    def info(self, *scenario_files:str) -> bool:
        """Print geometry information about given scenario files.

        Parameters
        ----------
        *scenario_files : str
            Strings that define CTSimU scenario files.

        Returns
        -------
        success : bool
        """

        for scenario_file in scenario_files:
            s = Scenario(scenario_file)
            s.set_frame(0)
            geo = s.current_geometry()
            print(geo.info())

        return True

    def correction(self, *metadata_files, **kwargs) -> bool:
        """Run flat-field and dark-field correction on projections
        defined in metadata files.

        Parameters
        ----------
        *metadata_files : str
            Strings that define CTSimU metadata files for uncorrected
            projection images.

        **kwargs : float

            + `ffRescaleFactor` : float

                Rescale factor after flat-field division.

                Standard value: `60000`

            + `ffRescaleOffset` : float

                Offset to be added to all pixels after flat-field correction
                and rescaling.

                Standard value: `0`

        Returns
        -------
        success : bool
        """
        ffRescaleFactor = float(60000.0)
        ffRescaleOffset = 0

        for key, value in kwargs.items():
            if key == "rescaleFactor":
                ffRescaleFactor = value
            elif key == "offsetAfterRescale":
                ffRescaleOffset = value

        for metadata_file in metadata_files:
            # Prepare a pipeline
            pipeline = self.get_ff_pipeline_from_metadata(metadata_file, ffRescaleFactor, ffRescaleOffset)
            pipeline.run()

        return True

    def get_ff_pipeline_from_metadata(self, metadata_file:str, ff_rescale_factor:float=60000, ff_rescale_offset:float=0) -> 'Pipeline':
        """Create a pipeline with a flat-field correction step based on the given metadata file.

        Parameters
        ----------
        metadata_file : str
            Path to a metadata file.

        ff_rescale_factor : float
            Rescale factor after flat-field division.

            Standard value: `60000`

        ff_rescale_offset : float
            Offset to be added to all pixels after flat-field correction
            and rescaling.

            Standard value: `0`

        Returns
        -------
        pipeline : ctsimu.processing.pipeline.Pipeline
            Pipeline with a flat-field correction step.
        """
        if os.path.isfile(metadata_file):
            if os.path.exists(metadata_file):
                metafileAbsolute = ""
                metafileAbsDir   = ""

                log("Metadata File: {mf}".format(mf=metadata_file))

                jsonDict = read_json_file(metadata_file)

                # Try to find the working directory with the metadata files:
                metafileAbsolute = os.path.abspath(metadata_file)
                metafileAbsDir   = os.path.dirname(metafileAbsolute)

                projFilename = get_value_or_none(jsonDict, "output", "projections", "filename")
                if projFilename is not None:
                    if not os.path.isabs(projFilename): # Check if an absolute path is provided
                        # If a relative path is provided, this path is
                        # relative to the location of the metadata file:
                        projFilename = join_dir_and_filename(metafileAbsDir, projFilename)

                projNumber = get_value_or_none(jsonDict, "output", "projections", "number")
                projDataType = get_value_or_none(jsonDict, "output", "projections", "datatype")
                projByteOrder = get_value_or_none(jsonDict, "output", "projections", "byteorder")
                projHeaderSizeFile = get_value_or_none(jsonDict, "output", "projections", "headersize", "file")
                projHeaderSizeImage = get_value_or_none(jsonDict, "output", "projections", "headersize", "image")
                width  = get_value_or_none(jsonDict, "output", "projections", "dimensions", "x", "value")
                height = get_value_or_none(jsonDict, "output", "projections", "dimensions", "y", "value")

                darkFilename = get_value_or_none(jsonDict, "output", "projections", "dark_field", "filename")
                if darkFilename is not None:
                    if not os.path.isabs(darkFilename): # Check if an absolute path is provided
                        # If a relative path is provided, this path is
                        # relative to the location of the metadata file:
                        darkFilename = join_dir_and_filename(metafileAbsDir, darkFilename)

                darkNumber    = get_value(jsonDict, ["output", "projections", "dark_field", "number"], 0)
                darkCorrected = get_value(jsonDict, ["output", "projections", "dark_field", "projections_corrected"], False)

                flatFilename  = get_value_or_none(jsonDict, "output", "projections", "flat_field", "filename")
                if flatFilename != None:
                    if not os.path.isabs(flatFilename): # Check if an absolute path is provided
                        # If a relative path is provided, this path is
                        # relative to the location of the metadata file:
                        flatFilename = join_dir_and_filename(metafileAbsDir, flatFilename)

                flatNumber    = get_value(jsonDict, ["output", "projections", "flat_field", "number"], 0)
                flatCorrected = get_value(jsonDict, ["output", "projections", "flat_field", "projections_corrected"], False)

                if projFilename is None:
                    raise Exception(f"No projection filename provided in metadata file: {metadata_file}")

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

                ffCorrection = Step_FlatFieldCorrection()
                ffCorrection.setFlatFieldRescaleFactor(ff_rescale_factor, ff_rescale_offset)

                pipeline = Pipeline()
                pipeline.setInputFileStack(projections)
                pipeline.setOutputFileStack(None)  # set later, if corrections will be applied
                ffCorrection.setOffsetCorrection(False)
                ffCorrection.setGainCorrection(False)

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
                        if darkFilename is not None:
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
                            ffCorrection.setDarkFileStack(darkStack)
                            pipeline.setOutputFileStack(outputFiles)

                if not flatCorrected:
                    if flatNumber > 0:
                        if flatFilename is not None:
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
                            ffCorrection.setFlatFileStack(flatStack)
                            pipeline.setOutputFileStack(outputFiles)

                pipeline.addStep(ffCorrection)

                return pipeline
            else:
                log(f"Metadata file does not exist: {metadata_file}")
                return False
        else:
            log(f"Invalid metadata file path: {metadata_file}")
            return False

    def post_process_recursive(self, directory:str, correction:bool=True, recon_configs:bool=True, overwrite:bool=False, openct_variant:str='free', openct_abspaths:bool=False, verbose:bool=True):
        """Run post-processing recursively on a whole directory. Searches for
        metadata files and automatically runs flat-field corrections and
        creates reconstruction config files.

        Parameters
        ----------
        directory : str
            Directory in which recursive post-processing will take place.

        correction : bool
            Run flat-field correction where applicable?

            Standard value: `True`

        recon_configs : bool
            Create reconstruction config files?

            Standard value: `True`

        overwrite : bool
            Overwrite existing files?

            Standard value: `False`

        openct_variant : str
            When reconstruction config files are created,
            the variant of the OpenCT config file.

            Possible values: `"free"`, `"circular"`

            Standard value: `"free"`

        openct_abspaths : bool
            Use absolute paths in OpenCT config file?

            Standard value: `False`

        verbose : bool
            Print text output?

            Standard value: `True`
        """
        pass

    def test_2D(self, test_name:str, *metadata, **kwargs) -> bool:
        """Run a projection-based 2D test on projections
        defined in one or more metadata files.

        Parameters
        ----------
        test_name : str
            Identifier of the CTSimU test, e.g. `"2D-SW-1"`.

        *metadata : str
            One or more paths to metadata files.

        **kwargs : str
            Possible subtest identifiers (as keys) and
            associated metadata files (as values).

        Returns
        -------
        success : bool
        """

        evaluationStep = None
        if test_name=="2D-FB-2":
            evaluationStep = Test2D_FB_2()
        elif test_name=="2D-FB-1":
            evaluationStep = Test2D_FB_1()
        elif test_name=="2D-SW-1":
            evaluationStep = Test2D_SW_1()
        elif test_name=="2D-SW-2":
            evaluationStep = Test2D_SW_2()
        elif test_name=="2D-SW-3":
            evaluationStep = Test2D_SW_3()
        elif test_name=="2D-SW-4":
            evaluationStep = Test2D_SW_4()
        elif test_name=="2D-DW-1":
            evaluationStep = Test2D_DW_1()
        elif test_name=="2D-WE-1":
            evaluationStep = Test2D_WE_1()
        elif test_name=="2D-WE-2":
            evaluationStep = Test2D_WE_2()
        elif test_name=="2D-HS-1":
            evaluationStep = Test2D_HS_1()
        else:
            log("Not a valid test name: {test_name}".format(test_name=test_name))
            return False

        # The flat-field rescale factor is fixed for all tests:
        ffRescaleFactor = float(60000.0)
        ffRescaleOffset = 0

        # Evaluate keyword list and see if metadata files are provided
        # together with keywords that might identify subtest scenarios.
        metadata_list = []
        for key, metafile in kwargs.items():
            if evaluationStep is not None:
                evaluationStep.addSubtest(key)
                metadata_list.append(metafile)
            else:
                log(f"'{test_name}' is not a valid command for the toolbox.")
                return False

        # Fallback to non-keyword list. Works only for tests that require no subtest identification.
        if len(metadata_list) == 0:
            metadata_list = list(metadata)

        for metadata_file in metadata_list:
            if os.path.isfile(metadata_file):
                if os.path.exists(metadata_file):
                    metafileAbsolute = ""
                    metafileAbsDir   = ""

                    # Try to find the working directory with the metadata files:
                    if metadata_file is not None:
                        metafileAbsolute = os.path.abspath(metadata_file)
                        metafileAbsDir   = os.path.dirname(metafileAbsolute)

                    resultFileDir = "{testname}-results".format(testname=test_name)
                    if metafileAbsDir is not None and len(metafileAbsDir) > 0:
                        resultFileDir = metafileAbsDir + "/" + resultFileDir

                    evaluationStep.setResultFileDirectory(resultFileDir)

                    pipeline = self.get_ff_pipeline_from_metadata(metadata_file, ffRescaleFactor, ffRescaleOffset)
                    pipeline.addStep(evaluationStep)

                    # Run pipeline: FF-correction and CTSimU test.
                    pipeline.run()

                    # Run the follow-up procedure to output the results.
                    evaluationStep.followUp()
                else:
                    raise Exception(f"Cannot access metadata file: {metadata_file}")
            else:
                raise Exception(f"Cannot access metadata file: {metadata_file}")

        return True