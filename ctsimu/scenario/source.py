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
	def __init__(self, _root=None):
		Part.__init__(self, name="source", _root=_root)

		# Extra geometry parameters:
		self.new_subgroup("geometry")
		self.geometry.set(key="type", value="cone", native_unit="string",
			valid_values=[None, "cone", "parallel"])
		self.geometry.new_subgroup("beam_divergence")
		self.geometry.beam_divergence.set(key="u", value=0, native_unit="deg")
		self.geometry.beam_divergence.set(key="v", value=0, native_unit="deg")

		# X-ray source parameters
		self.set(key="model",         value=None, native_unit="string", simple=True)
		self.set(key="manufacturer",  value=None, native_unit="string", simple=True)
		self.set(key="voltage",       value=None, native_unit="kV")
		self.set(key="current",       value=None, native_unit="mA")

		# Target
		self.new_subgroup("target")
		self.target.set(key="material_id",     value=None,         native_unit="string", simple=True)
		self.target.set(key="type",            value="reflection", native_unit="string", simple=True, valid_values=[None, "reflection", "transmission"])
		self.target.set(key="thickness",       value=None,         native_unit="mm")
		self.target.new_subgroup("angle")
		self.target.angle.set(key="incidence", value=None,         native_unit="deg")
		self.target.angle.set(key="emission",  value=None,         native_unit="deg")

		# Spot
		self.new_subgroup("spot")
		self.spot.new_subgroup("size")
		self.spot.size.set(key="u", value=0, native_unit="mm")
		self.spot.size.set(key="v", value=0, native_unit="mm")
		self.spot.size.set(key="w", value=0, native_unit="mm")

		self.spot.new_subgroup("sigma")
		self.spot.sigma.set(key="u", value=0, native_unit="mm")
		self.spot.sigma.set(key="v", value=0, native_unit="mm")
		self.spot.sigma.set(key="w", value=0, native_unit="mm")

		self.spot.new_subgroup("intensity_map")
		self.spot.intensity_map.set(key="file",       value=None, native_unit="string")
		self.spot.intensity_map.set(key="dim_x",      value=None, simple=True)
		self.spot.intensity_map.set(key="dim_y",      value=None, simple=True)
		self.spot.intensity_map.set(key="dim_z",      value=None, simple=True)
		self.spot.intensity_map.set(key="type",       value="uint16", native_unit="string", simple=True)
		self.spot.intensity_map.set(key="endian",     value="little", native_unit="string", simple=True)
		self.spot.intensity_map.set(key="headersize", value=0, simple=True)

		# Spectrum
		self.new_subgroup("spectrum")
		self.spectrum.set(key="monochromatic", value=False, native_unit="bool", simple=True)
		self.spectrum.set(key="file",          value=None,  native_unit="string")

		# Window
		self.new_subgroup("window", array=True)
		self.window.set(key="material_id", value=None, native_unit="string", simple=True)
		self.window.set(key="thickness",   value=None, native_unit="mm")

		# Filters
		self.new_subgroup("filters", array=True)
		self.filters.set(key="material_id", value=None, native_unit="string", simple=True)
		self.filters.set(key="thickness",   value=None, native_unit="mm")

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

		self.geometry.set_from_json(geo)

		# Source properties:
		sourceprops = json_extract(json_scenario, ["source"])
		Group.set_from_json(self, sourceprops)

	def json_dict(self) -> dict:
		"""Create a CTSimU JSON dictionary for this source object and its subgroups.

		Returns
		-------
		json_dict : dict
		"""
		jd = Group.json_dict(self)

		# Remove the geometry section, it was only used to store
		# the type and beam_divergence. It is not part of the
		# actual 'source' definition in a CTSimU scenario, but of
		# the source geometry, which is handled differently by the toolbox.
		del jd["geometry"]

		return jd

	def geometry_dict(self) -> dict:
		"""Create a dictionary of the source geometry for a CTSimU scenario file.

		Returns
		-------
		source_geometry_dict : dict
			Dictionary with the geometry of this part.
		"""
		jd = self.geometry.json_dict()
		jd.update(Part.geometry_dict(self))
		return jd