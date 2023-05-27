# -*- coding: UTF-8 -*-
"""
Generic sample: position and orientation, size and material properties.
"""

from ..helpers import *
from ..geometry import *
from .part import Part
from .parameter import Parameter

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

		self.set(key="surface_mesh_file", value="", native_unit="string")

		# Mesh file unit of length:
		self.set(key="unit", value="mm", native_unit="string")
		self.set(key="scaling_factor_r", value=1.0, native_unit="")
		self.set(key="scaling_factor_s", value=1.0, native_unit="")
		self.set(key="scaling_factor_t", value=1.0, native_unit="")
		self.set(key="material_id", value=None, native_unit="string")

	def reset(self):
		Part.reset()

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

		# Surface mesh file:
		self.set_parameter_from_key("surface_mesh_file", json_object, ["file"], fail_value=None)
		self.set_parameter_from_key("unit", json_object, ["unit"], fail_value="mm")
		self.set_parameter_from_key("scaling_factor_r", json_object, ["scaling_factor", "r"], fail_value=1.0)
		self.set_parameter_from_key("scaling_factor_s", json_object, ["scaling_factor", "s"], fail_value=1.0)
		self.set_parameter_from_key("scaling_factor_t", json_object, ["scaling_factor", "t"], fail_value=1.0)
		self.set_parameter_from_key("material_id", json_object, ["material_id"], fail_value=None)