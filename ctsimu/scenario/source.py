# -*- coding: UTF-8 -*-
"""
A CTSimU X-ray source: position, orientation, size, and other parameters.
"""

from ..helpers import *
from ..geometry import *
from .part import Part
from .parameter import Parameter
from .group import Group, Array

class Source(Part):
	"""CTSimU X-ray source."""
	def __init__(self):
		Part.__init__(self, "source")

		self.source_geometry_extras = Group()
		self.source_geometry_extras.set(key="type", value="cone", native_unit="string",
			valid_values=[None, "cone", "parallel"])
		self.beam_divergence = Group("beam_divergence")
		self.beam_divergence.set(key="u", value=0, native_unit="deg")
		self.beam_divergence.set(key="v", value=0, native_unit="deg")
		self.source_geometry_extras.add_subgroup(self.beam_divergence)

		# X-ray source parameters
		self.set(key="model",         value=None, native_unit="string", simple=True)
		self.set(key="manufacturer",  value=None, native_unit="string", simple=True)
		self.set(key="voltage",       value=None, native_unit="kV")
		self.set(key="current",       value=None, native_unit="mA")

		# Target
		self.target = Group("target")
		self.target.set(key="material_id",     value=None,         native_unit="string", simple=True)
		self.target.set(key="type",            value="reflection", native_unit="string", simple=True, valid_values=[None, "reflection", "transmission"])
		self.target.set(key="thickness",       value=None,         native_unit="mm")
		self.target_angle = Group("angle")
		self.target_angle.set(key="incidence", value=None,         native_unit="deg")
		self.target_angle.set(key="emission",  value=None,         native_unit="deg")
		self.target.add_subgroup(self.target_angle)
		self.add_subgroup(self.target)

		# Spot
		self.spot = Group("spot")
		self.spot_size = Group("size")
		self.spot_size.set(key="u", value=0, native_unit="mm")
		self.spot_size.set(key="v", value=0, native_unit="mm")
		self.spot_size.set(key="w", value=0, native_unit="mm")
		self.spot.add_subgroup(self.spot_size)

		self.spot_sigma = Group("sigma")
		self.spot_sigma.set(key="u", value=0, native_unit="mm")
		self.spot_sigma.set(key="v", value=0, native_unit="mm")
		self.spot_sigma.set(key="w", value=0, native_unit="mm")
		self.spot.add_subgroup(self.spot_sigma)

		self.intensity_map = Group("intensity_map")
		self.intensity_map.set(key="file",       value=None, native_unit="string")
		self.intensity_map.set(key="dim_x",      value=None, simple=True)
		self.intensity_map.set(key="dim_y",      value=None, simple=True)
		self.intensity_map.set(key="dim_z",      value=None, simple=True)
		self.intensity_map.set(key="type",       value="uint16", native_unit="string", simple=True)
		self.intensity_map.set(key="endian",     value="little", native_unit="string", simple=True)
		self.intensity_map.set(key="headersize", value=0, simple=True)
		self.spot.add_subgroup(self.intensity_map)

		self.add_subgroup(self.spot)

		# Spectrum
		self.spectrum = Group("spectrum")
		self.spectrum.set(key="monochromatic", value=False, native_unit="bool", simple=True)
		self.spectrum.set(key="file",          value=None,  native_unit="string")
		self.add_subgroup(self.spectrum)

		# Window
		self.window = Array("window")
		self.window.set(key="material_id", value=None, native_unit="string", simple=True)
		self.window.set(key="thickness",   value=None, native_unit="mm")
		self.add_subgroup(self.window)

		# Filters
		self.filters = Array("filters")
		self.filters.set(key="material_id", value=None, native_unit="string", simple=True)
		self.filters.set(key="thickness",   value=None, native_unit="mm")
		self.add_subgroup(self.filters)

	def check(self):
		# Check if the target type is valid:
		if not (self.target.get("type") in valid_xray_target_types):
			raise ValueError(f"Not a valid X-ray source target type: \'{self.target.get('type')}\'. Should be any of {valid_xray_target_types}.")
			return False

		return True

	def set_from_json(self, json_scenario:dict):
		"""Import the X-ray source definition and geometry from the JSON object.
		The JSON object should contain the complete content
		of the scenario definition file
		(at least the geometry and detector sections).

		Parameters
		----------
		json_scenario : dict
			A complete CTSimU scenario object, as imported from a JSON structure.
		"""
		self.reset()

		# Extract the X-ray source's geometry:
		geo = json_extract(json_scenario, ["geometry", "source"])
		self.set_geometry(json_geometry_object=geo, proper_cs="local")

		self.source_geometry_extras.set_from_json(geo)

		# Source properties:
		sourceprops = json_extract(json_scenario, ["source"])
		Group.set_from_json(self, sourceprops)

	def geometry_dict(self) -> dict:
		"""Create a dictionary of the source geometry for a CTSimU scenario file.

		Returns
		-------
		source_geometry_dict : dict
			Dictionary with the geometry of this part.
		"""
		jd = self.source_geometry_extras.json_dict()
		jd.update(Part.geometry_dict(self))
		return jd