# -*- coding: UTF-8 -*-
"""
Generic sample: position and orientation, size and material properties.
"""

from ..helpers import *
from ..geometry import *
from .part import Part
from .parameter import Parameter
from .group import Group

class Sample(Part):
	"""Generic sample."""
	def __init__(self, name:str=""):
		"""A name can be passed when initializing the sample.

		Parameters
		----------
		name : str
			Sample name.
		"""
		Part.__init__(self, name)

		self.set(key="name", value=name, native_unit="string", simple=True)
		self.set(key="file", value=None, native_unit="string")  # surface mesh file

		# Mesh file unit of length:
		self.set(key="unit", value="mm", native_unit="string", simple=True)

		self.scaling_factor = Group("scaling_factor")
		self.scaling_factor.set(key="r", value=1.0, native_unit=None)
		self.scaling_factor.set(key="s", value=1.0, native_unit=None)
		self.scaling_factor.set(key="t", value=1.0, native_unit=None)
		self.add_subgroup(self.scaling_factor)

		self.set(key="material_id", value=None, native_unit="string", simple=True)

	def set_from_json(self, json_object:dict, stage_coordinate_system:'CoordinateSystem'=None):
		"""Import the sample geometry from the JSON sample object.
		The `stage_coordinate_system` must be given as a
		`ctsimu.geometry.CoordinateSystem` object. If this part is not attached
		to the stage, `None` can be passed instead.

		Parameters
		----------
		json_object : dict
			A CTSimU sample object, as imported from a JSON structure.

		stage_coordinate_system : CoordinateSystem
			The stage coordinate system. `None` is also accepted if this part is
			not attached to the stage.
		"""
		self.reset()
		self.set_name(get_value(json_object, ["name"], fail_value="Sample"))

		# Extract the sample's geometry:
		geo = json_extract(json_object, ["position"])
		self.set_geometry(geo, stage_coordinate_system)

		Group.set_from_json(self, json_object)

	def json_dict(self) -> dict:
		"""Create a dictionary of this sample for a CTSimU JSON file.

		Returns
		-------
		json_dict : dict
			The sample's JSON dictionary.
		"""

		jd = Part.json_dict(self)
		jd["position"] = Part.geometry_dict(self)

		return jd