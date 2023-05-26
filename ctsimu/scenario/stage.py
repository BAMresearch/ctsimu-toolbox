# -*- coding: UTF-8 -*-
"""
A CTSimU sample stage: position, orientation.
"""

from ..helpers import *
from ..geometry import *
from .part import Part
from .parameter import Parameter

class Stage(Part):
	"""CTSimU sample stage."""
	def __init__(self, name:str=""):
		"""A name can be passed when initializing the stage.

		Parameters
		----------
		name : str
			Stage name.
		"""
		Part.__init__(self, name)

	def set_from_json(self, json_scenario:dict):
		"""Import the stage definition and geometry from the JSON object.
		The JSON object should contain the complete content
		of the scenario definition file
		(at least the geometry section).

		Parameters
		----------
		json_scenario : dict
			A complete CTSimU scenario object, as imported from a JSON structure.
		"""
		self.reset()

		# Extract the stage's geometry:
		geo = json_extract(json_scenario, ["geometry", "stage"])
		self.set_geometry(geo)