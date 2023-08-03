# -*- coding: UTF-8 -*-
"""Tools to set up, read and write [CTSimU scenarios].
[CTSimU scenarios]: https://bamresearch.github.io/ctsimu-scenarios
"""

from ..helpers import *
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
		self.detector = Detector()
		self.source   = Source()
		self.stage    = Stage()
		self.samples  = list()

		self.file = File()
		self.environment = Environment()
		self.acquisition = Acquisition()
		self.materials = list()
		self.simulation = None # simply imported as dict

	def read(self, file:str=None, json_dict:dict=None):
		if file is not None:
			json_dict = read_json_file(filename=file)
		elif not isinstance(json_dict, dict):
			raise Exception("Scenario: read() function expects either a filename as a string or a CTSimU JSON dictionary as a Python dict.")
			return False

		self.detector.set_from_json(json_dict)
		self.source.set_from_json(json_dict)
		self.stage.set_from_json(json_dict)

		json_samples = json_extract(json_dict, ["samples"])
		if json_samples is not None:
			for json_sample in json_samples:
				s = Sample()
				s.set_from_json(json_sample, self.stage.coordinate_system)
				self.samples.append(s)

		self.file.set_from_json(json_extract(json_dict, [self.file.name]))
		self.environment.set_from_json(json_extract(json_dict, [self.environment.name]))
		self.acquisition.set_from_json(json_extract(json_dict, [self.acquisition.name]))
		self.simulation = json_extract(json_dict, ["simulation"])

		json_materials = json_extract(json_dict, ["materials"])
		for json_material in json_materials:
			m = Material()
			m.set_from_json(json_material)
			self.materials.append(m)

	def write(self, file:str=None):
		if file is not None:
			self.file.file_format_version.set("major", 1)
			self.file.file_format_version.set("minor", 2)

			write_json_file(filename=file, dictionary=self.json_dict())

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