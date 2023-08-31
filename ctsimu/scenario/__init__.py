# -*- coding: UTF-8 -*-
"""Tools to set up, read and write [CTSimU scenarios].
[CTSimU scenarios]: https://bamresearch.github.io/ctsimu-scenarios
"""

import math
import os
from datetime import datetime

from ..helpers import *
from ..geometry import *
from . import *

from .detector import Detector
from .source import Source
from .stage import Stage
from .sample import Sample
from .file import File
from .environment import Environment
from .acquisition import Acquisition
from .material import Material
from .metadata import Metadata

class Scenario:
    def __init__(self, filename:str=None):
        self.detector = Detector(_root=self)
        self.source   = Source(_root=self)
        self.stage    = Stage(_root=self)
        self.samples  = list()

        self.file = File(_root=self)
        self.environment = Environment(_root=self)
        self.acquisition = Acquisition(_root=self)
        self.materials = list()
        self.simulation = None # simply imported as dict

        # CTSimU metadata structure:
        # Necessary for generation of reconstruction configs.
        self.metadata = Metadata(_root=self)

        # Subgroups are necessary for the 'get' and 'set' command:
        self.subgroups = [
            self.detector,
            self.source,
            self.stage,
            self.file,
            self.environment,
            self.acquisition,
            self.metadata
        ]

        self.current_frame = 0
        self.current_scenario_path = None
        self.current_scenario_file = None
        self.current_scenario_basename = None
        self.current_scenario_directory = None

        self.current_metadata_path = None
        self.current_metadata_file = None
        self.current_metadata_basename = None
        self.current_metadata_directory = None
        self.metadata_is_set = False

        if filename is not None:
            self.read(filename=filename)

    def read(self, filename:str=None, json_dict:dict=None):
        """Import a CTSimU scenario from a file or a given
        scenario dictionary.

        Parameters
        ----------
        filename : str
            Path to a CTSimU scenario file.

            Default value: `None`

        json_dict : dict
            Provide a dictionary with a scenario structure instead of a file.

            Default value: `None`
        """
        self.current_scenario_path = None

        if filename is not None:
            json_dict = read_json_file(filename=filename)
            self.current_scenario_path = filename
            self.current_scenario_directory = os.path.dirname(filename)
            self.current_scenario_file = os.path.basename(filename)
            self.current_scenario_basename, extension = os.path.splitext(self.current_scenario_file)
        elif not isinstance(json_dict, dict):
            raise Exception("Scenario: read() function expects either a filename as a string or a CTSimU JSON dictionary as a Python dict.")
            return False

        self.detector.set_from_json(json_dict)
        self.source.set_from_json(json_dict)
        self.stage.set_from_json(json_dict)

        json_samples = json_extract(json_dict, ["samples"])
        if json_samples is not None:
            for json_sample in json_samples:
                s = Sample(_root=self)
                s.set_from_json(json_sample, self.stage.coordinate_system)
                self.samples.append(s)

        self.file.set_from_json(json_extract(json_dict, [self.file._name]))
        self.environment.set_from_json(json_extract(json_dict, [self.environment._name]))
        self.acquisition.set_from_json(json_extract(json_dict, [self.acquisition._name]))
        self.simulation = json_extract(json_dict, ["simulation"])

        json_materials = json_extract(json_dict, ["materials"])
        for json_material in json_materials:
            m = Material(_root=self)
            m.set_from_json(json_material)
            self.materials.append(m)

        self.set_frame(0, reconstruction=False)

        if not self.metadata_is_set:
            self.create_default_metadata()

    def reset_metadata(self):
        """Reset scenario's metadata information."""
        self.current_metadata_path = None
        self.current_metadata_file = None
        self.current_metadata_basename = None
        self.current_metadata_directory = None
        self.metadata_is_set = False

        # Create new, empty metadata:
        self.metadata = Metadata()

    def read_metadata(self, filename:str=None, json_dict:dict=None, import_referenced_scenario:bool=False):
        """Import metadata from a CTSimU metadata file or a given
        metadata dictionary.

        Parameters
        ----------
        filename : str
            Path to a CTSimU metadata file.

            Default value: `None`

        json_dict : dict
            Provide a dictionary with a metadata structure instead of a file.

            Default value: `None`

        import_referenced_scenario : bool
            Import the scenario JSON file that's referenced in the metadata file?
            Generates a warning if this fails.

            The scenario definition will be searched at two locations in the following order:

            1. Try to read from external file defined by `acquisition_geometry.path_to_CTSimU_JSON`.

            2. Try to read embedded scenario definition from
            `simulation.ctsimu_scenario`. Note that external drift files
            specified in the scenario will likely fail to load. An error
            will be issued in this case.

            Default value: `False`
        """
        if filename is not None:
            json_dict = read_json_file(filename=filename)

            # If a file is read, we want to make sure that it is a valid
            # and supported metadata file:
            if isinstance(json_dict, dict):
                file_type = get_value(json_dict, ["file", "file_type"])
                if file_type != "CTSimU Metadata":
                    raise Exception(f"Invalid metadata structure: the string 'CTSimU Metadata' was not found in 'file.file_type' in the metadata file {file}.")
            else:
                raise Exception(f"Error when reading the metadata file: {filename}")

            self.current_metadata_path = filename
            self.current_metadata_directory = os.path.dirname(filename)
            self.current_metadata_file = os.path.basename(filename)
            self.current_metadata_basename, extension = os.path.splitext(self.current_metadata_file)

        if json_dict is not None:
            # If we get a `json_dict` as function parameter, we do not
            # test for a valid version because reduced/simplified metadata
            # structures should be supported as well.
            self.metadata.set_from_json(json_dict)
            self.metadata_is_set = True

            if import_referenced_scenario:
                # Import the scenario that's referenced in the metadata file.
                ref_file = self.metadata.get(["acquisition_geometry", "path_to_CTSimU_JSON"])

                import_success = False

                try:
                    if (ref_file is not None) and (ref_file != ""):
                        if isinstance(ref_file, str):
                            ref_file = abspath_of_referenced_file(self.current_metadata_path, ref_file)
                        else:
                            raise Exception("read_metadata: path_to_CTSimU_JSON is not a string.")

                        # Try to import scenario:
                        self.read(filename=ref_file)
                        import_success = True
                except Exception as e:
                    warnings.warn(str(e))
                    import_success = False

                if not import_success:
                    # Try to import the embedded scenario structure.
                    if json_exists_and_not_null(json_dict, ["simulation", "ctsimu_scenario"]):
                        self.read(json_dict=json_dict["simulation"]["ctsimu_scenario"])
                        import_success = True

                # Create default metadata in case the original
                # metadata file did not contain all information that's needed:
                self.create_default_metadata()
                # Re-import metadata:
                self.metadata.set_from_json(json_dict)


    def create_default_metadata(self):
        """Set default metadata from scenario information,
        to use if no metadata file is available."""
        self.reset_metadata()

        geo = self.current_geometry()
        cera_parameters = geo.get_CERA_standard_circular_parameters()

        # Basename:
        if self.current_scenario_basename is not None:
            self.current_metadata_basename = self.current_scenario_basename

        basename = self.current_metadata_basename
        n_projections = self.acquisition.get("number_of_projections")
        projection_filename = f"{basename}_{counter_format(n_projections)}.raw"

        # Prepare filename for dark fields:
        n_darks = self.acquisition.dark_field.get("number")
        dark_filename = None
        if n_darks is not None:
            if n_darks > 0:
                if n_darks > 1:
                    dark_filename = f"{basename}_dark_{counter_format(n_darks)}.raw"
                else:
                    dark_filename = f"{basename}_dark.raw"

        # Prepare filename for flat fields:
        n_flats = self.acquisition.flat_field.get("number")
        flat_filename = None
        if n_flats is not None:
            if n_flats > 0:
                if n_flats > 1:
                    flat_filename = f"{basename}_flat_{counter_format(n_flats)}.raw"
                else:
                    flat_filename = f"{basename}_flat.raw"


        n_cols = self.detector.get("columns")
        n_rows = self.detector.get("rows")
        pixel_size_u = self.detector.pixel_pitch.get("u")
        pixel_size_v = self.detector.pixel_pitch.get("v")

        voxels_x = n_cols
        voxels_y = n_cols
        voxels_z = n_rows

        voxelsize_x = cera_parameters["voxelsize"]["x"]
        voxelsize_y = cera_parameters["voxelsize"]["y"]
        voxelsize_z = cera_parameters["voxelsize"]["z"]

        now = datetime.now()

        metadata = {
            "file": {
                "name": basename,
                "description": self.file.get("description"),

                "contact": self.file.get("contact"),
                "date_created": now.strftime("%Y-%m-%d"),
                "date_changed": now.strftime("%Y-%m-%d"),

                "file_type": "CTSimU Metadata",
                "file_format_version": {
                    "major": ctsimu_supported_metadata_version["major"],
                    "minor": ctsimu_supported_metadata_version["minor"]
                }
            },

            "output": {
                "system": None,
                "date_measured": None,
                "projections": {
                    "filename": projection_filename,
                    "number": n_projections,
                    "frame_average": self.acquisition.get("frame_average"),
                    "max_intensity": self.detector.get(["gray_value", "imax"]),
                    "datatype": "uint16",
                    "byteorder": "little",
                    "headersize": {
                        "file": 0,
                        "image": 0
                    },
                    "dimensions": {
                        "x": {"value": n_cols, "unit": "px"},
                        "y": {"value": n_rows, "unit": "px"}
                    },
                    "pixelsize": {
                        "x": {"value": pixel_size_u, "unit": "mm"},
                        "y": {"value": pixel_size_v, "unit": "mm"}
                    },
                    "dark_field": {
                        "number": n_darks,
                        "frame_average": self.acquisition.dark_field.get("frame_average"),
                        "filename": dark_filename,
                        "projections_corrected": False
                    },
                    "flat_field": {
                        "number": n_flats,
                        "frame_average": self.acquisition.flat_field.get("frame_average"),
                        "filename": flat_filename,
                        "projections_corrected": False
                    },
                    "bad_pixel_map": {
                        "filename": None,
                        "projections_corrected": False
                    }
                },
                "tomogram":
                {
                    "filename":  f"{basename}_recon.raw",
                    "datatype":  "float32",
                    "byteorder": "little",
                    "headersize": {
                        "file": 0,
                        "image": 0
                    },
                    "dimensions": {
                        "x": {"value": voxels_x, "unit": "px"},
                        "y": {"value": voxels_y, "unit": "px"},
                        "z": {"value": voxels_z, "unit": "px"}
                    },
                    "voxelsize": {
                        "x": {"value": voxelsize_x, "unit": "mm"},
                        "y": {"value": voxelsize_y, "unit": "mm"},
                        "z": {"value": voxelsize_z, "unit": "mm"}
                    }
                }
            },

            "acquisition_geometry": {
                "path_to_CTSimU_JSON": self.current_scenario_path
            },

            "reconstruction": {
                "software": None,
                "settings": { }
            },

            "simulation": {
                "ctsimu_scenario": None
            }
        }

        self.read_metadata(json_dict=metadata, import_referenced_scenario=False)

    def write(self, filename:str):
        """Write a scenario JSON file.

        Parameters
        ----------
        filename : str
            Filename for the scenario file.
        """
        if filename is not None:
            self.file.file_format_version.set("major", ctsimu_supported_scenario_version["major"])
            self.file.file_format_version.set("minor", ctsimu_supported_scenario_version["minor"])

            write_json_file(filename=filename, dictionary=self.json_dict())

    def write_metadata(self, filename:str):
        """Write a metadata JSON file for the scenario.

        Parameters
        ----------
        filename : str
            Filename for the metadata file.
        """
        if filename is not None:
            metadata_dict = self.metadata.json_dict()
            # potentially add simulation.ctsimu_scenario here:
            metadata_dict["simulation"]["ctsimu_scenario"] = self.json_dict()
            write_json_file(filename=filename, dictionary=metadata_dict)


    def get(self, key:list) -> float | str | bool:
        """Get the current value of the parameter identified by a list of keys.

        Parameters
        ----------
        key : list
            List of strings that identify the key of the requested
            parameter within the CTSimU scenario structure.

        Returns
        -------
        value : float or str or bool
            Current value of the requested parameter.
        """
        if isinstance(key, list):
            # Special treatment for the source geometry extras: type or beam_divergence:
            if len(key) > 2:
                if key[0:2] == ["geometry", "source"]:
                    return self.source.source_geometry_extras.get(key[2:])

            # Standard treatment:
            if len(key) > 1:
                for s in self.subgroups:
                    if s._name == key[0]:
                        return s.get(key[1:])

        raise Exception(f"Error in get: key not found: {key}")

    def path_of_external_file(self, filename:str) -> str:
        """Get the path of an external file referenced in the currently
        imported JSON scenario.

        Parameters
        ----------
        filename : str
            Possibly relative file path from JSON scenario file.

        Returns
        -------
        abs_path : str
            Absolute path to the referred external file.
        """
        if os.path.isabs(filename):
            # Already absolute path?
            return filename

        if self.current_scenario_path is not None:
            if isinstance(self.current_scenario_path, str):
                json_dirname = os.path.dirname(self.current_scenario_path)
                filename = f"{json_dirname}/{filename}"

        # On fail, simply return the filename.
        return filename

    def json_dict(self) -> dict:
        """Create a CTSimU JSON dictionary from the scenario.

        Returns
        -------
        json_dict : dict
        """
        jd = dict()
        jd["file"]        = self.file.json_dict()
        jd["environment"] = self.environment.json_dict()

        jd["geometry"]    = dict()
        jd["geometry"]["detector"] = self.detector.geometry_dict()

        jd["geometry"]["source"] = self.source.geometry_dict()

        jd["geometry"]["stage"] = self.stage.geometry_dict()

        jd["detector"] = self.detector.json_dict()
        jd["source"]   = self.source.json_dict()
        jd["samples"]  = []
        for sample in self.samples:
            jd["samples"].append(sample.json_dict())

        jd["acquisition"] = self.acquisition.json_dict()
        jd["materials"] = []
        for material in self.materials:
            jd["materials"].append(material.json_dict())

        jd["simulation"] = self.simulation

        return jd

    def n_frames(self) -> int:
        """Number of frames in the scenario.

        Returns
        -------
        n_frames : int
            Number of frames in the scenario.
        """

        # 'Frame' is in this context a projection image.
        # However, if we assume frame averaging, the number of
        # frames could also be: n_frames = nProjection * nFrameAverages
        n_frames = self.acquisition.get("number_of_projections")
        return n_frames

    def get_current_stage_rotation_angle(self):
        """Stage rotation angle (in deg) for the current frame.

        Returns
        -------
        stage_rotation_angle : float
            Current stage rotation angle (in deg).
        """

        start_angle = float(self.acquisition.get("start_angle"))
        stop_angle  = float(self.acquisition.get("stop_angle"))
        nPositions  = float(self.n_frames())

        # If the final projection is taken at the stop angle
        # (and not one step before), the number of positions
        # has to be decreased by 1, resulting in one less
        # angular step being performed.
        if self.acquisition.get("include_final_angle") is True:
            if nPositions > 0:
                nPositions -= 1

        angular_range = 0
        if start_angle <= stop_angle:
            angular_range = stop_angle - start_angle
        else:
            raise Exception("The start angle cannot be greater than the stop angle. Scan direction must be specified by the acquisition 'direction' keyword (CCW or CW).")

        angular_position = start_angle
        if nPositions != 0:
            angular_position = start_angle + self.current_frame*angular_range/nPositions

        # Mathematically negative:
        if self.acquisition.get("direction") == "CW":
            angular_position = -angular_position

        return angular_position

    def set_frame(self, frame:float=0, reconstruction:bool=False):
        self.current_frame = frame

        # Number of frames:
        n_frames = self.n_frames()

        stage_deg = self.get_current_stage_rotation_angle()
        stage_rot = math.radians(stage_deg)

        # Update materials:
        for material in self.materials:
            material.set_frame(frame, n_frames, reconstruction)

        # Update stage, source, detector and other parameters:
        self.stage.set_frame(frame, n_frames, stage_rot, None, reconstruction)
        self.source.set_frame(frame, n_frames, 0, None, reconstruction)
        self.detector.set_frame(frame, n_frames, 0, None, reconstruction)

        self.file.set_frame(frame, n_frames, reconstruction)
        self.environment.set_frame(frame, n_frames, reconstruction)
        self.acquisition.set_frame(frame, n_frames, reconstruction)

        # Update samples:
        stage_cs = self.stage.coordinate_system
        for sample in self.samples:
            sample.set_frame(frame, n_frames, 0, stage_cs, reconstruction)

    def current_geometry(self) -> 'ctsimu.geometry.Geometry':
        """Return a `ctsimu.geometry.Geometry` object for the
        current setup of the scenario.

        Returns
        -------
        geometry : ctsimu.geometry.Geometry
            Geometry for current frame.
        """
        geo = Geometry()
        geo.detector.copy_cs(self.detector.coordinate_system)
        geo.source.copy_cs(self.source.coordinate_system)
        geo.stage.copy_cs(self.stage.coordinate_system)

        geo.detector.set_size(
            pixels_u=self.detector.get("columns"),
            pixels_v=self.detector.get("rows"),
            pitch_u=self.detector.pixel_pitch.get("u"),
            pitch_v=self.detector.pixel_pitch.get("v")
        )

        return geo

    def write_recon_VGI(self, name:str="", volume_filename:str="", vgi_filename:str=None):
        """Write a VGI file for the reconstruction volume such that it can be loaded with VGSTUDIO."""

        voxels_x = self.metadata.output.get(["tomogram", "dimensions", "x"])
        voxels_y = self.metadata.output.get(["tomogram", "dimensions", "y"])
        voxels_z = self.metadata.output.get(["tomogram", "dimensions", "z"])

        voxelsize_x = self.metadata.output.get(["tomogram", "voxelsize", "x"])
        voxelsize_y = self.metadata.output.get(["tomogram", "voxelsize", "y"])
        voxelsize_z = self.metadata.output.get(["tomogram", "voxelsize", "z"])

        output_datatype = self.metadata.output.get(["tomogram", "datatype"])
        if output_datatype == "uint16":
            dataTypeOutput = "unsigned integer"
            bits = 16
            datarangelower = 0
            datarangeupper = -1
        else:
            dataTypeOutput = "float"
            bits = 32
            datarangelower = -1
            datarangeupper = 1

        vgi_content = f"""{{volume1}}
[representation]
size = {voxels_x} {voxels_y} {voxels_z}
datatype = {dataTypeOutput}
datarange = {datarangelower} {datarangeupper}
bitsperelement = {bits}
[file1]
SkipHeader = 0
FileFormat = raw
Size = {voxels_x} {voxels_y} {voxels_z}
Name = {volume_filename}
Datatype = {dataTypeOutput}
datarange = {datarangelower} {datarangeupper}
BitsPerElement = {bits}
{{volumeprimitive12}}
[geometry]
resolution = {voxelsize_x} {voxelsize_y} {voxelsize_z}
unit = mm
[volume]
volume = volume1
[description]
text = {name}"""

        if vgi_filename is not None:
            touch_directory(vgi_filename)
            with open(vgi_filename, 'w') as f:
                f.write(vgi_content)
                f.close()

        return vgi_content

    def write_CERA_config(self, save_dir:str=None, basename:str=None, create_vgi:bool=False):
        """Write CERA reconstruction config files.

        Parameters
        ----------
        save_dir : str
            Folder where to place the CERA config files. This is meant to be the
            same directory where the reconstruction metadata file is located,
            such that relative paths will match.

            If `None` is given, a directory will be inferred:

            - If only a JSON scenario file was imported to set up the scenario,
            the config files will be stored in a subdirectory next to the
            JSON scenario file of the following pattern:

                `{json_scenario_basename}/reconstruction`

            - If a reconstruction metadata file was imported, the config files
            will be stored next to the metadata file.

            Default value: `None`

        basename : str
            Base name for the created files. If `None` is given, the base
            name will be inferred from the scenario's metadata.

            Default value: `None`

        create_vgi : bool
            Write VGI file for future reconstruction volume?
        """

        matrices = []

        if basename is None:
            # Extract base name from metadata
            metadata_basename = self.metadata.get(["file", "name"])
            if metadata_basename is not None:
                basename = f"{metadata_basename}_recon_cera"
            else:
                basename = f"recon_cera"

        # Projection files
        n_projections = self.acquisition.get("number_of_projections")
        projection_file_pattern = self.metadata.get(["output", "projections", "filename"])
        projection_datatype = self.metadata.get(["output", "projections", "datatype"])
        projection_file_byteorder = self.metadata.get(["output", "projections", "byteorder"])
        projection_headersize = self.metadata.get(["output", "projections", "headersize", "file"])

        projection_filetype = "tiff"
        if projection_file_pattern.lower().endswith(".raw"):
            projection_filetype = "raw"

        # Acquisition
        start_angle = self.acquisition.get("start_angle")
        stop_angle  = self.acquisition.get("stop_angle")
        total_angle = stop_angle - start_angle

        for p in range(n_projections):
            self.set_frame(frame=p, reconstruction=True)

            # CERA projection matrix for projection p:
            m = self.current_geometry().projection_matrix(mode="CERA")
            matrices.append(m)

        # Go back to frame zero:
        self.set_frame(frame=0, reconstruction=True)

        if create_vgi:
            vgi_filename = join_dir_and_filename(save_dir, f"{basename}.vgi")
            volume_filename = f"{basename}.raw"

            self.write_recon_VGI(vgi_filename=vgi_filename, name=basename, volume_filename=volume_filename)

        create_CERA_config(
            geo=self.current_geometry(),
            projection_file_pattern=projection_file_pattern,
            basename=basename,
            save_dir=save_dir,
            n_projections=n_projections,
            projection_datatype=projection_datatype,
            projection_filetype=projection_filetype,
            projection_byteorder=projection_file_byteorder,
            projection_headersize=projection_headersize,
            start_angle=0,  # do not compensate the scenario start angle in the reconstruction
            total_angle=total_angle,
            scan_direction=self.acquisition.get("direction"),
            voxels_x=self.metadata.output.tomogram.dimensions.x.get(),
            voxels_y=self.metadata.output.tomogram.dimensions.y.get(),
            voxels_z=self.metadata.output.tomogram.dimensions.z.get(),
            voxelsize_x=self.metadata.output.tomogram.voxelsize.x.get(),
            voxelsize_y=self.metadata.output.tomogram.voxelsize.y.get(),
            voxelsize_z=self.metadata.output.tomogram.voxelsize.z.get(),
            i0max=self.metadata.output.get(["projections", "max_intensity"]),
            output_datatype=convert(cera_converter["datatype"], self.metadata.output.get(["tomogram", "datatype"])),
            matrices=matrices
        )

    def write_OpenCT_config(self, save_dir:str=None, basename:str=None, create_vgi:bool=False, variant:str='free', abspaths:bool=False):
        """Write OpenCT reconstruction config files.

        Parameters
        ----------
        save_dir : str
            Folder where to place the CERA config files. This is meant to be the
            same directory where the reconstruction metadata file is located,
            such that relative paths will match.

            If `None` is given, a directory will be inferred:

            - If only a JSON scenario file was imported to set up the scenario,
            the config files will be stored in a subdirectory next to the
            JSON scenario file of the following pattern:

                `{json_scenario_basename}/reconstruction`

            - If a reconstruction metadata file was imported, the config files
            will be stored next to the metadata file.

            Default value: `None`

        basename : str
            Base name for the created config file. If `None` is given, the base
            name will be inferred from the scenario's metadata.

            Default value: `None`

        create_vgi : bool
            Write VGI file for future reconstruction volume?

        variant : str
            Which variant of the OpenCT file format will be created: free trajectory
            or circular trajectory.

            Possible values: `"free"`, `"circular"`

            Default value: `"free"`

        abspaths : bool
            Set to `True` if absolute paths should be used in the OpenCT
            config file.

            Default value: `False` (relative paths)

        Returns
        -------
        openct_dict : dict
            Dictionary with the OpenCT JSON structure.
        """

        matrices = []

        if basename is None:
            # Extract base name from metadata
            metadata_basename = self.metadata.get(["file", "name"])
            if metadata_basename is not None:
                basename = f"{metadata_basename}_recon_openCT"
            else:
                basename = f"recon_openCT"

        # Name of config file:
        openct_config_filename = join_dir_and_filename(save_dir, f"{basename}.json")

        def projections_from_pattern(json_dict:dict):
            n = get_value(json_dict, ["number"]) # number of images
            pattern  = None
            filedir  = None
            filename = None
            filelist = list()

            if n is not None:
                if n > 0:
                    pattern = get_value(json_dict, ["filename"])
                    if abspaths is True:
                        pattern = abspath_of_referenced_file(openct_config_filename, pattern)

                    if pattern is not None:
                        filedir, filename = os.path.split(pattern)

                    # Generate list of projection files:
                    if filename is not None:
                        # Auto-generate projection file list.
                        # List of sequentially numbered projection images,
                        # starting at 0000.
                        for p in range(n):
                            if '%' in filename:
                                filelist.append(filename % p)
                            else:
                                filelist.append(filename)

            return n, filedir, filename, filelist

        # Projection files
        n_projections, projection_filedir, projection_filename, projection_filelist = projections_from_pattern(self.metadata.output.projections.json_dict())

        projection_datatype = self.metadata.get(["output", "projections", "datatype"])
        projection_file_byteorder = self.metadata.get(["output", "projections", "byteorder"])
        projection_headersize = self.metadata.get(["output", "projections", "headersize", "file"])

        projection_filetype = "tiff"
        if projection_filename is not None:
            if projection_filename.lower().endswith(".raw"):
                projection_filetype = "raw"

        # Dark files; only is corrections need to be applied.
        openct_dark_image = None
        if not self.metadata.output.projections.dark_field.projections_corrected.get() is True:
            n_darks, dark_filedir, dark_filename, dark_filelist = projections_from_pattern(self.metadata.output.projections.dark_field.json_dict())
            if isinstance(dark_filelist, list):
                if len(dark_filelist) > 0:
                    openct_dark_image = join_dir_and_filename(dark_filedir, dark_filelist[0])

        # Bright files; only is corrections need to be applied.
        flat_filedir = None
        flat_filelist = None
        if not self.metadata.output.projections.flat_field.projections_corrected.get() is True:
            n_flats, flat_filedir, flat_filename, flat_filelist = projections_from_pattern(self.metadata.output.projections.flat_field.json_dict())

        # Acquisition
        start_angle = self.acquisition.get("start_angle")
        stop_angle  = self.acquisition.get("stop_angle")
        total_angle = stop_angle - start_angle

        for p in range(n_projections):
            self.set_frame(frame=p, reconstruction=True)

            # CERA projection matrix for projection p:
            m = self.current_geometry().projection_matrix(mode="OpenCT")
            matrices.append(m)

        # Go back to frame zero:
        self.set_frame(frame=0, reconstruction=True)

        volume_filename = f"{basename}.img"
        if create_vgi:
            vgi_filename = join_dir_and_filename(save_dir, f"{basename}.vgi")

            self.write_recon_VGI(vgi_filename=vgi_filename, name=basename, volume_filename=volume_filename)

        openct_dict = create_OpenCT_config(
            geo=self.current_geometry(),
            filename=openct_config_filename,
            variant=variant,
            projection_files=projection_filelist,
            projection_dir=projection_filedir,
            projection_datatype=projection_datatype,
            projection_filetype=projection_filetype,
            projection_headersize=projection_headersize,
            projection_byteorder=projection_file_byteorder,
            total_angle=total_angle,
            scan_direction=self.acquisition.get("direction"),
            matrices=matrices,
            volumename=volume_filename,
            voxels_x=self.metadata.output.tomogram.dimensions.x.get(),
            voxels_y=self.metadata.output.tomogram.dimensions.y.get(),
            voxels_z=self.metadata.output.tomogram.dimensions.z.get(),
            voxelsize_x=self.metadata.output.tomogram.voxelsize.x.get(),
            voxelsize_y=self.metadata.output.tomogram.voxelsize.y.get(),
            voxelsize_z=self.metadata.output.tomogram.voxelsize.z.get(),
            bright_image_dir=flat_filedir,
            bright_images=flat_filelist,
            dark_image=openct_dark_image)

        return openct_dict