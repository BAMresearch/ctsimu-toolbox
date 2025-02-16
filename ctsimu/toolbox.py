# -*- coding: UTF-8 -*-
"""
Quick metadata-based post-processing with simple commands.

.. include:: ./toolbox.md
"""
import os
import shutil
from datetime import datetime

from .helpers import *
from .geometry import Geometry
from .image import ImageStack
from .scenario import Scenario
from .processing.pipeline import Pipeline
from .processing.flat_field import Step_FlatFieldCorrection
from .evaluation.test2D_FB_2 import Test2D_FB_2
from .evaluation.test2D_FB_1 import Test2D_FB_1
from .evaluation.test2D_SW_1 import Test2D_SW_1
from .evaluation.test2D_SW_2 import Test2D_SW_2
from .evaluation.test2D_SW_3 import Test2D_SW_3
from .evaluation.test2D_SW_4 import Test2D_SW_4
from .evaluation.test2D_DW_1 import Test2D_DW_1
from .evaluation.test2D_WE_1 import Test2D_WE_1
from .evaluation.test2D_WE_2 import Test2D_WE_2
from .evaluation.test2D_HS_1 import Test2D_HS_1
from .evaluation.testDigTwin import TestDigTwin
from .responses.measurands import Measurands

class Toolbox:
    """ Manages a test run, including preliminary flat field corrections, based on metadata JSON input files. """
    def __init__(self, operation, *args, **kwargs):
        if operation == "info":
            self.info(*args)
        elif operation == "correction":
            self.correction(*args, **kwargs)
        elif operation == "standardize":
            self.standardize(*args)
        elif operation == "recon_config":
            self.recon_config(*args, **kwargs)
        elif (operation == "post-processing") or (operation == "post_processing"):
            self.post_processing(*args, **kwargs)
        elif operation.startswith("2D-"):
            # Possibly a 2D test... try it!
            self.test_2D(operation, *args, **kwargs)
        elif operation == "testDigitalTwin":
            self.testDigitalTwin(*args, **kwargs)

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

            + `rescaleFactor` : float

                Rescale factor after flat-field division. If `None`, the value
                will be imported from the key `max_intensity` in the metadata
                file, or set to `60000` if this fails.

                Standard value: `None`

            + `offsetAfterRescale` : float

                Offset to be added to all pixels after flat-field correction
                and rescaling.

                Standard value: `0`

            + `overwrite` : bool

                Overwrite existing, corrected projection images?

                Standard value: `True`

        Returns
        -------
        success : bool
        """

        # Default values for keyword arguments:
        settings = {
            "rescaleFactor": None,
            "offsetAfterRescale": 0,
            "overwrite": True
        }

        # Change default settings for keyword arguments that are set:
        for key, value in kwargs.items():
            if key in settings:
                settings[key] = value

        for metadata_file in metadata_files:
            # Prepare a pipeline
            try:
                pipeline = self.get_ff_pipeline_from_metadata(metadata_file, settings["rescaleFactor"], settings["offsetAfterRescale"])
                pipeline.run(overwrite=settings["overwrite"])
            except Exception as e:
                log(f"Error: {metadata_file}: {str(e)}")

        return True

    def get_ff_pipeline_from_metadata(self, metadata_file:str, rescaleFactor:float=None, offsetAfterRescale:float=0) -> 'Pipeline':
        """Create a pipeline with a flat-field correction step based on the given metadata file.

        Parameters
        ----------
        metadata_file : str
            Path to a metadata file.

        rescaleFactor : float
            Rescale factor after flat-field division. If `None`, the value
            will be imported from key `max_intensity` in metadata file,
            or set to `60000` if this fails.

            Standard value: `None`

        offsetAfterRescale : float
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

                if rescaleFactor is None:
                    rescaleFactor = get_value(jsonDict, ["output", "projections", "max_intensity"], 60000)

                darkFilename = get_value_or_none(jsonDict, "output", "projections", "dark_field", "filename")
                if darkFilename is not None:
                    if not os.path.isabs(darkFilename): # Check if an absolute path is provided
                        # If a relative path is provided, this path is
                        # relative to the location of the metadata file:
                        darkFilename = join_dir_and_filename(metafileAbsDir, darkFilename)

                darkNumber    = get_value(jsonDict, ["output", "projections", "dark_field", "number"], 0)
                darkCorrected = get_value(jsonDict, ["output", "projections", "dark_field", "projections_corrected"], False)

                flatFilename  = get_value_or_none(jsonDict, "output", "projections", "flat_field", "filename")
                if flatFilename is not None:
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
                ffCorrection.setFlatFieldRescaleFactor(rescaleFactor, offsetAfterRescale)

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
                raise Exception(f"Metadata file does not exist: {metadata_file}")
        else:
            raise Exception(f"Invalid metadata file path: {metadata_file}")

    def standardize(self, *scenario_files, **kwargs):
        """Standardize CTSimU scenario files to the
        file format version currently supported by the toolbox.

        The previous version of the file will be kept in a
        copy called {filename}.json.old

        Parameters
        ----------
        *scenario_files : str
            Scenario files to be standardized.
        """
        settings = {

        }

        # Change default settings for keyword arguments that are set:
        for key, value in kwargs.items():
            if key in settings:
                settings[key] = value

        now = datetime.now()

        for scenario_file in scenario_files:
            if os.path.isfile(scenario_file):
                if os.path.exists(scenario_file):
                    log(f"Standardizing {scenario_file} ...")
                    scenario_file_old = scenario_file + ".old"
                    if os.path.exists(scenario_file_old):
                        # An old version seems to already exist.
                        # Skip standardization of this file.
                        log(f"   Skipped standardization, old file exists:")
                        log(f"   {scenario_file_old}")
                        continue
                    else:
                        try:
                            # Try to read scenario:
                            s = Scenario(scenario_file)

                            # Backup previous file:
                            shutil.copy2(scenario_file, scenario_file_old)

                            # Check if backup exists:
                            if os.path.exists(scenario_file_old):
                                s.file.date_changed.set(now.strftime("%Y-%m-%d"))
                                s.file.file_format_version.major.set(ctsimu_supported_scenario_version["major"])
                                s.file.file_format_version.minor.set(ctsimu_supported_scenario_version["minor"])
                                s.write(scenario_file)
                            else:
                                raise Exception(f"Failed to create scenario backup file: {scenario_file_old}")
                        except Exception as e:
                            log(f"   Error: {str(e)}")
                else:
                    log(f"   Error: file not found: '{scenario_file}'")
            else:
                log(f"   Error: not a scenario file: '{scenario_file}'")

    def recon_config(self, *metadata_files, **kwargs):
        """Create reconstruction configuration files for
        given metadata files.

        Parameters
        ----------
        *metadata_files : str
            Metadata files that describe reconstruction data.

        **kwargs

            + `openct` : `bool`

                Create OpenCT config file?

                Standard value: `True`

            + `cera` : `bool`

                Create CERA config file?

                Standard value: `True`

            + `openct_variant` : `str`

                When reconstruction config files are created,
                the variant of the OpenCT config file.

                Possible values: `"free"`, `"circular"`

                Standard value: `"free"`

            + `openct_abspaths` : `bool`

                Use absolute paths in OpenCT config file?

                Standard value: `False`

            + `create_vgi` : `bool`

                Create VGI file for the reconstruction volume?

                Standard value: `True`

            + `overwrite` : `bool`

                Overwrite existing output files?

                Standard value: `False`
        """

        # Default values for keyword arguments:
        settings = {
            "openct": True,
            "cera": True,
            "openct_variant": "free",
            "openct_abspaths": False,
            "create_vgi": True,
            "overwrite": False
        }

        # Change default settings for keyword arguments that are set:
        for key, value in kwargs.items():
            if key in settings:
                settings[key] = value

        for metadata in metadata_files:
            try:
                log(f"Creating recon config for: {metadata}")
                s = Scenario()
                s.read_metadata(filename=metadata, import_referenced_scenario=True)

                recon_config_dir, recon_config_metafile = os.path.split(metadata)

                basename, extension = os.path.splitext(recon_config_metafile)

                # Remove '_metadata' from basename:
                basename = basename.replace("_metadata", "")

                if settings["cera"] is True:
                    cera_filename = basename + "_cera.config"
                    cera_filepath = join_dir_and_filename(recon_config_dir, cera_filename)
                    cera_write = True
                    if os.path.exists(cera_filepath):
                        if settings["overwrite"] is False:
                            cera_write = False

                    if cera_write is True:
                        log(f"  Writing CERA config files to:   '{recon_config_dir}' ...")
                        s.write_CERA_config(
                            save_dir=recon_config_dir,
                            basename=f"{basename}_cera",
                            create_vgi=settings["create_vgi"]
                            )

                if settings["openct"] is True:
                    openct_filename = basename + "_openCT.json"
                    openct_filepath = join_dir_and_filename(recon_config_dir, openct_filename)
                    openct_write = True
                    if os.path.exists(openct_filepath):
                        if settings["overwrite"] is False:
                            openct_write = False

                    if openct_write is True:
                        log(f"  Writing OpenCT config files to: '{recon_config_dir}' ...")
                        s.write_OpenCT_config(
                            save_dir=recon_config_dir,
                            basename=f"{basename}_openCT",
                            create_vgi=settings["create_vgi"],
                            variant=settings["openct_variant"],
                            abspaths=settings["openct_abspaths"]
                            )
            except Exception as e:
                log(f"   Error: {metadata}: {str(e)}")

    def post_processing(self, *directories, **kwargs):
        """Run post-processing recursively on whole directories.

        Searches for metadata files in the given directories
        (and their subdirectories) and automatically runs
        flat-field corrections and creates reconstruction config files.

        Parameters
        ----------
        *directories : str
            Directories in which recursive post-processing will take place.

        **kwargs

            + `correction` : `bool`

                Run flat-field correction where applicable?

                Standard value: `False`

            + `rescaleFactor` : `float`

                Rescale factor after flat-field division. If `None`, the value
                will be imported from the key `max_intensity` in the metadata
                file, or set to `60000` if this fails.

                Standard value: `None`

            + `offsetAfterRescale` : `float`

                Offset to be added to all pixels after flat-field correction
                and rescaling.

                Standard value: `0`

            + `recon_config` : `bool`

                Create reconstruction config files?

                Standard value: `False`

            + `openct_variant` : `str`

                When reconstruction config files are created,
                the variant of the OpenCT config file.

                Possible values: `"free"`, `"circular"`

                Standard value: `"free"`

            + `openct_abspaths` : `bool`

                Use absolute paths in OpenCT config file?

                Standard value: `False`

            + `standardize` : `bool`

                Standardize CTSimU scenario files to the
                file format version currently supported by the toolbox.

                The previous version of the file will be kept in a
                copy called {filename}.json.old

                Standard value: `False`

            + `overwrite` : `bool`

                Overwrite existing output files?

                Does not apply to standardization: the old scenario backup
                files will never be overwritten.

                Standard value: `False`
        """

        # Default values for keyword arguments:
        settings = {
            "correction": False,
            "rescaleFactor": None,
            "offsetAfterRescale": 0,
            "recon_config": False,
            "cera": True,
            "openct": True,
            "openct_variant": "free",
            "openct_abspaths": False,
            "standardize": False,
            "overwrite": False
        }

        # Change default settings for keyword arguments that are set:
        for key, value in kwargs.items():
            if key in settings:
                settings[key] = value

        # Walk through directories and collect paths of scenario files,
        # projection metadata and reconstruction metadata files.
        scenario_files = list()
        metadata_projection = list()
        metadata_reconstruction = list()
        for directory in directories:
            print(f"Searching directory: {directory}")
            walker = os.walk(directory)
            for dirpath, dirnames, filenames in walker:
                json_files = [f for f in filenames if f.endswith(".json")]
                # Check JSON files:
                for j in json_files:
                    try:
                        jsonpath = join_dir_and_filename(dirpath, j)
                        jd = read_json_file(jsonpath)
                        file_type = get_value(jd, ["file", "file_type"])

                        if file_type is not None:
                            if file_type == "CTSimU Scenario":
                                scenario_files.append(jsonpath)
                                log(f"   Found scenario file:           {j}")
                            elif file_type == "CTSimU Metadata":
                                if json_exists(jd, ["output", "tomogram", "filename"]):
                                    # Seems to be a reconstruction file.
                                    metadata_reconstruction.append(jsonpath)
                                    log(f"   Found reconstruction metadata: {j}")
                                else:
                                    # Assume a projection metadata file.
                                    metadata_projection.append(jsonpath)
                                    log(f"   Found projection metadata:     {j}")
                    except Exception as e:
                        log(f"Error: {str(e)}")

        if settings["standardize"] is True:
            # Standardize scenario files.
            for scenario_file in scenario_files:
                try:
                    self.standardize(scenario_file, **settings)
                except Exception as e:
                    log(f"Error: {str(e)}")

        if settings["recon_config"] is True:
            # Create reconstruction configs:
            for metadata in metadata_reconstruction:
                self.recon_config(metadata, **settings)

        if settings["correction"] is True:
            # Run flat-field corrections:
            for metadata in metadata_projection:
                self.correction(metadata, **settings)

        print("Done.")

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
        rescaleFactor = float(60000.0)
        offsetAfterRescale = 0

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
                    if metafileAbsDir is not None:
                        if isinstance(metafileAbsDir, str):
                            if len(metafileAbsDir) > 0:
                                resultFileDir = metafileAbsDir + "/" + resultFileDir

                    evaluationStep.setResultFileDirectory(resultFileDir)

                    pipeline = self.get_ff_pipeline_from_metadata(metadata_file, rescaleFactor, offsetAfterRescale)
                    pipeline.addStep(evaluationStep)

                    # Run pipeline: FF-correction and CTSimU test.
                    pipeline.run()
                else:
                    raise Exception(f"Cannot access metadata file: {metadata_file}")
            else:
                raise Exception(f"Cannot access metadata file: {metadata_file}")

        # Run the follow-up procedure to output the results.
        if evaluationStep is not None:
            evaluationStep.followUp()

        return True

    def testDigitalTwin(self, *metadata_files, **kwargs):
        """Run the digital twin test on dimensional measurements
        defined in a metadata file.

        Parameters
        ----------
        *metadata : str
            One or more paths to metadata files.

        **kwargs : str
            Possible subtest identifiers (as keys) and
            associated metadata files (as values).
        """

        #for key, value in kwargs.items():
        #    if key in settings:
        #        settings[key] = value

        for metadata_file in metadata_files:
            # Prepare a pipeline
            try:
                #print(metadata_file)
                dataprep = TestDigTwin(metadata_file)
                #print(bla)
                
                self.RealValues = dataprep.read_and_filter_csv_files(dataprep.real_folder_path, "real_ct")
                #print(self.RealValues)
                
                self.SimValues = dataprep.read_and_filter_csv_files(dataprep.sim_folder_path, "simulation_ct")
                #print(self.SimValues)
                #self.testDigTwin(metadata_file)
                #pipeline.run(overwrite=settings["overwrite"])
                dataprep.En_calc(self.RealValues, self.SimValues)
                
                dataprep.plotResults()
                dataprep.TwinTest_report()

            except Exception as e:
                log(f"Error: {metadata_file}: {str(e)}")
        #evaluationStep = testDigTwin(metadata)

