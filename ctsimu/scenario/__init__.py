# -*- coding: UTF-8 -*-
"""Tools to set up, read and write [CTSimU scenarios].
[CTSimU scenarios]: https://bamresearch.github.io/ctsimu-scenarios
"""

import math

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
    def __init__(self):
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
        self.current_scenario_file = None
        self.current_scenario_directory = None
        self.current_metadata_directory = None
        self.metadata_is_set = False

    def read(self, file:str=None, json_dict:dict=None):
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
        self.current_scenario_file = None

        if file is not None:
            json_dict = read_json_file(filename=file)
            self.current_scenario_file = file
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

        self.file.set_from_json(json_extract(json_dict, [self.file.name]))
        self.environment.set_from_json(json_extract(json_dict, [self.environment.name]))
        self.acquisition.set_from_json(json_extract(json_dict, [self.acquisition.name]))
        self.simulation = json_extract(json_dict, ["simulation"])

        json_materials = json_extract(json_dict, ["materials"])
        for json_material in json_materials:
            m = Material(_root=self)
            m.set_from_json(json_material)
            self.materials.append(m)

        self.set_frame(0, reconstruction=False)

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

            Default value: `False`
        """
        self.metadata_is_set = False

        if filename is not None:
            json_dict = read_json_file(filename=filename)

            # If a file is read, we want to make sure that it is a valid
            # and supported metadata file:
            if isinstance(json_dict, dict):
                file_type = get_value(json_dict, ["file", "file_type"])
                if file_type != "CTSimU Metadata":
                    raise Exception(f"Invalid metadata structure: the string 'CTSimU Metadata' was not found in 'file.file_type' in the metadata file {filename}.")
            else:
                raise Exception(f"Error when reading the metadata file: {filename}")

        if json_dict is not None:
            # If we get a `json_dict` as function parameter, we do not
            # test for a valid version because reduced/simplified metadata
            # structures should be supported as well.
            self.metadata.set_from_json(json_dict)
            self.metadata_is_set = True

    def write(self, file:str=None):
        if file is not None:
            self.file.file_format_version.set("major", ctsimu_supported_scenario_version["major"])
            self.file.file_format_version.set("minor", ctsimu_supported_scenario_version["minor"])

            write_json_file(filename=file, dictionary=self.json_dict())

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
                    if s.name == key[0]:
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

        if self.current_scenario_file is not None:
            if isinstance(self.current_scenario_file, str):
                json_dirname = os.path.dirname(self.current_scenario_file)
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

    def create_recon_VGI(self, name:str="", volume_filename:str="", vgi_filename:str=None):
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

    def write_CERA_config(self, save_dir:str=None, basename=None, metadata:dict=None, metadata_file:str=None, create_vgi:bool=False):
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
            name will be inferred from the given metadata.

            Default value: `None`

        metadata : dict
            A CTSimU metadata dictionary that defines the tomogram output parameters.
            Only necessary if no metadata has been loaded yet. Alternatively,
            a path to a metadata file can be provided (see below).

            Default value: `None`

        metadata_file : str
            Path to a CTSimU metadata file that defines the tomogram output parameters.
            Only necessary if no metadata has been loaded yet. Alternatively,
            a path to a metadata dictionary can be provided (see above).

            Default value: `None`
        """

        matrices = []

        if metadata_file is not None:
            metadata = read_json_file(metadata_file)

        if metadata is not None:
            self.metadata.set_from_json(metadata)

        if basename is None:
            # Extract base name from metadata
            basename = self.metadata.get(["file", "name"])
            basename += "_recon_cera"

        # Projection files
        n_projections = self.acquisition.get("number_of_projections")
        projection_file_pattern = self.metadata.get(["output", "projections", "filename"])
        projection_file_datatype = self.metadata.get(["output", "projections", "datatype"])
        projection_file_byteorder = self.metadata.get(["output", "projections", "byteorder"])
        projection_file_headersize = self.metadata.get(["output", "projections", "headersize", "file"])

        projection_file_type = "tiff"
        if projection_file_pattern.lower().endswith(".raw"):
            if projection_file_datatype == "uint16":
                projection_file_type = "raw_uint16"
            elif projection_file_datatype == "float32":
                projection_file_type = "raw_float"
            else:
                raise Exception(f"Projection datatype not supported by CERA: {projection_file_datatype}. Supported datatypes: 'uint16' and 'float32'.")

        big_endian = False
        if projection_file_byteorder == "big":
            big_endian = True

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
            vgi_filename = f"{save_dir}/{basename}.vgi"
            volume_filename = f"{basename}.raw"

            self.create_recon_VGI(vgi_filename=vgi_filename, name=basename, volume_filename=volume_filename)

        create_CERA_config(
            geo=self.current_geometry(),
            projection_file_pattern=projection_file_pattern,
            basename=basename,
            save_dir=save_dir,
            n_projections=n_projections,
            projection_file_type=projection_file_type,
            start_angle=0,  # do not compensate the scenario start angle in the reconstruction
            total_angle=total_angle,
            scan_direction=self.acquisition.get("direction"),
            i0max=self.metadata.output.get(["projections", "max_intensity"]),
            big_endian=big_endian,
            raw_header_size=projection_file_headersize,
            output_datatype=convert(cera_converter["datatype"], self.metadata.output.get(["tomogram", "datatype"])),
            matrices=matrices
        )

    def write_OpenCT_config(self, save_dir:str=".", basename=None, metadata:dict=None, metadata_file:str=None, create_vgi:bool=False, variant='free'):
        """Write OpenCT reconstruction config files.

        Parameters
        ----------
        save_dir : str
            Folder where to place the OpenCT config file. This is meant to be the
            same directory where the reconstruction metadata file is located,
            such that relative paths will match.

        basename : str
            Base name for the created files. If `None` is given, the base
            name will be inferred from the given metadata.

            Default value: `None`

        metadata : dict
            A CTSimU metadata dictionary that defines the tomogram output parameters.
            Only necessary if no metadata has been loaded yet. Alternatively,
            a path to a metadata file can be provided (see below).

            Default value: `None`

        metadata_file : str
            Path to a CTSimU metadata file that defines the tomogram output parameters.
            Only necessary if no metadata has been loaded yet. Alternatively,
            a path to a metadata dictionary can be provided (see above).

            Default value: `None`

        variant : str
            Which variant of the OpenCT file format will be created: free trajectory
            or circular trajectory.

            Possible values: `"free"`, `"circular"`

            Default value: `"free"`
        """

        matrices = []

        if metadata_file is not None:
            metadata = read_json_file(metadata_file)

        if metadata is not None:
            self.metadata.set_from_json(metadata)

        if basename is None:
            # Extract base name from metadata
            basename = self.metadata.get(["file", "name"])
            basename += "_recon_openCT"

        filename = f"{save_dir}/{basename}.json"

        # Projection files
        n_projections = self.acquisition.get("number_of_projections")
        projection_file_pattern = self.metadata.get(["output", "projections", "filename"])
        projection_file_datatype = self.metadata.get(["output", "projections", "datatype"])
        projection_file_byteorder = self.metadata.get(["output", "projections", "byteorder"])
        projection_file_headersize = self.metadata.get(["output", "projections", "headersize", "file"])

        projection_file_type = "tiff"
        if projection_file_pattern.lower().endswith(".raw"):
            if projection_file_datatype == "uint16":
                projection_file_type = "raw_uint16"
            elif projection_file_datatype == "float32":
                projection_file_type = "raw_float"
            else:
                raise Exception(f"Projection datatype not supported by CERA: {projection_file_datatype}. Supported datatypes: 'uint16' and 'float32'.")

        big_endian = False
        if projection_file_byteorder == "big":
            big_endian = True

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
            vgi_filename = f"{save_dir}/{basename}.vgi"
            volume_filename = f"{basename}.raw"

            self.create_recon_VGI(vgi_filename=vgi_filename, name=basename, volume_filename=volume_filename)

        create_OpenCT_config(
            geo=self.current_geometry(),
            projection_file_pattern=projection_file_pattern,
            basename=basename,
            save_dir=save_dir,
            n_projections=n_projections,
            projection_file_type=projection_file_type,
            start_angle=0,  # do not compensate the scenario start angle in the reconstruction
            total_angle=total_angle,
            scan_direction=self.acquisition.get("direction"),
            i0max=self.metadata.output.get(["projections", "max_intensity"]),
            big_endian=big_endian,
            raw_header_size=projection_file_headersize,
            output_datatype=convert(cera_converter["datatype"], self.metadata.output.get(["tomogram", "datatype"])),
            matrices=matrices
        )
