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

		self.current_frame = 0
		self.current_json_file = None

	def read(self, file:str=None, json_dict:dict=None):
		self.current_json_file = None

		if file is not None:
			json_dict = read_json_file(filename=file)
			self.current_json_file = file
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

	def write(self, file:str=None):
		if file is not None:
			self.file.file_format_version.set("major", 1)
			self.file.file_format_version.set("minor", 2)

			write_json_file(filename=file, dictionary=self.json_dict())

	def get(self, key:list) -> float | str | bool:
		"""Get a parameter value identified by a list of key strings.

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
	    """Get the path of an external file referred to in the currently
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

	    if self.current_json_file is not None:
	        if isinstance(self.current_json_file, str):
	            json_dirname = os.path.dirname(self.current_json_file)
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

	def n_frames(self):
		"""Number of frames in the scenario.

		Returns
		-------
		nFrames : int
			Number of frames in the scenario.
		"""

		# 'Frame' is in this context a projection image.
		# However, if we assume frame averaging, the number of
		# frames could also be: nFrames = nProjection * nFrameAverages
		nFrames = self.acquisition.get("number_of_projections")
		return nFrames

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
		nFrames = self.n_frames()

		stage_deg = self.get_current_stage_rotation_angle()
		stage_rot = math.radians(stage_deg)

		# Update materials:
		for material in self.materials:
			material.set_frame(frame, nFrames, reconstruction)

		# Update stage, source, detector and other parameters:
		self.stage.set_frame(frame, nFrames, stage_rot, None, reconstruction)
		self.source.set_frame(frame, nFrames, 0, None, reconstruction)
		self.detector.set_frame(frame, nFrames, 0, None, reconstruction)

		self.file.set_frame(frame, nFrames, reconstruction)
		self.environment.set_frame(frame, nFrames, reconstruction)
		self.acquisition.set_frame(frame, nFrames, reconstruction)

		# Update samples:
		stage_cs = self.stage.coordinate_system
		for sample in self.samples:
			sample.set_frame(frame, nFrames, 0, stage_cs, reconstruction)

	def current_geometry(self) -> 'ctsimu.geometry.Geometry':
		"""Return a 'ctsimu.geometry.Geometry' object for the
		current setup of the scenario.

		Returns
		-------
		geometry : ctsimu.geometry.Geometry
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
