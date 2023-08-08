# -*- coding: UTF-8 -*-
"""
A CTSimU detector: position, orientation, size, and other parameters.
"""

from ..helpers import *
from ..geometry import *
from .part import Part
from .parameter import Parameter
from .group import Group, Array

class Detector(Part):
	"""CTSimU detector."""
	def __init__(self):
		Part.__init__(self, "detector")

		# Detector parameters
		self.set(key="model",            value=None,    native_unit="string", simple=True)
		self.set(key="manufacturer",     value=None,    native_unit="string", simple=True)
		self.set(key="type",             value="ideal", native_unit="string", simple=True, valid_values=[None, "ideal", "real"])
		self.set(key="columns",          value=1000,    native_unit="px")
		self.set(key="rows",             value=1000,    native_unit="px")

		self.pixel_pitch = Group("pixel_pitch")
		self.pixel_pitch.set(key="u",    value=0.1,     native_unit="mm")
		self.pixel_pitch.set(key="v",    value=0.1,     native_unit="mm")
		self.add_subgroup(self.pixel_pitch)

		self.set(key="bit_depth",        value=16,      native_unit=None)
		self.set(key="integration_time", value=1.0,     native_unit="s")
		self.set(key="dead_time",        value=1.0,     native_unit="s")
		self.set(key="image_lag",        value=0.0,     native_unit=None)
		self.set(key="gain",             value=None,    native_unit=None)

		# Properties for gray value reproduction:
		self.gray_value = Group("gray_value")
		self.gray_value.add_alternative_name("grey_value")
		self.gray_value.set(key="imin",   value=0,     native_unit=None)
		self.gray_value.set(key="imax",   value=60000, native_unit=None)
		self.gray_value.set(key="factor", value=None,  native_unit="1/J")
		self.gray_value.set(key="offset", value=None,  native_unit=None)
		self.gray_value.set(key="intensity_characteristics_file",  value=None, native_unit="string")
		self.gray_value.set(key="efficiency_characteristics_file", value=None, native_unit="string")
		self.add_subgroup(self.gray_value)

		# Noise:
		self.noise = Group("noise")
		self.noise.set(key="snr_at_imax",                value=None, native_unit=None)
		self.noise.set(key="noise_characteristics_file", value=None, native_unit="string")
		self.add_subgroup(self.noise)

		# Unsharpness:
		self.unsharpness = Group("unsharpness")
		self.unsharpness.add_alternative_name("sharpness")
		self.unsharpness.set(key="basic_spatial_resolution", value=None, native_unit="mm")
		self.unsharpness.set(key="mtf", value=None, native_unit="string")
		self.add_subgroup(self.unsharpness)

		# Bad pixel map:
		self.bad_pixel_map = Group("bad_pixel_map")
		self.bad_pixel_map.set(key="file",       value=None, native_unit="string")
		# dim_x and dim_y not necessary. Must match detector columns/rows.
		#self.bad_pixel_map.set(key="dim_x",      value=None, simple=True)
		#self.bad_pixel_map.set(key="dim_y",      value=None, simple=True)
		self.bad_pixel_map.set(key="type",       value="int16",  native_unit="string", simple=True)
		self.bad_pixel_map.set(key="endian",     value="little", native_unit="string", simple=True)
		self.bad_pixel_map.set(key="headersize", value=0, simple=True)
		self.add_subgroup(self.bad_pixel_map)

		# Scintillator
		self.scintillator = Group("scintillator")
		self.scintillator.set(key="material_id", value=None, native_unit="string", simple=True)
		self.scintillator.set(key="thickness",   value=None, native_unit="mm")
		self.add_subgroup(self.scintillator)

		# Window
		self.window = Group("window")
		self.window_front = Array("front")
		self.window_front.set(key="material_id", value=None, native_unit="string", simple=True)
		self.window_front.set(key="thickness",   value=None, native_unit="mm")
		self.window.add_subgroup(self.window_front)
		self.window_rear = Array("rear")
		self.window_rear.set(key="material_id", value=None, native_unit="string", simple=True)
		self.window_rear.set(key="thickness",   value=None, native_unit="mm")
		self.window.add_subgroup(self.window_rear)
		self.add_subgroup(self.window)

		# Filters
		self.filters = Group("filters")
		self.filters_front = Array("front")
		self.filters_front.set(key="material_id", value=None, native_unit="string", simple=True)
		self.filters_front.set(key="thickness",   value=None, native_unit="mm")
		self.filters.add_subgroup(self.filters_front)
		self.filters_rear = Array("rear")
		self.filters_rear.set(key="material_id", value=None, native_unit="string", simple=True)
		self.filters_rear.set(key="thickness",   value=None, native_unit="mm")
		self.filters.add_subgroup(self.filters_rear)
		self.add_subgroup(self.filters)

	def check(self):
		# Check if the detector type is valid:
		if not (self.get("type") in valid_detector_types):
			raise ValueError(f"Not a valid detector type: \'{self.get('type')}\'. Should be any of {valid_detector_types}.")
			return False

		return True

	def physical_width(self) -> float:
		"""Get detector's physical width in mm.

		Returns
		-------
		physical_width : float
		"""
		return float(self.get("pitch_u") * int(self.get("columns")))

	def physical_height(self) -> float:
		"""Get detector's physical height in mm.

		Returns
		-------
		physical_height : float
		"""
		return float(self.get("pitch_v") * int(self.get("rows")))

	def pixel_area(self) -> float:
		"""Get area of a single pixel in mm².

		Returns
		-------
		area : float
			Pixel area in mm².
		"""
		return self.get("pitch_u") * self.get("pitch_v")

	def pixel_area_m2(self) -> float:
		"""Get area of a single pixel in m².

		Returns
		-------
		area : float
			Pixel area in m².
		"""
		return self.pixel_area() * 1.0e-6

	def max_gray_value(self) -> int:
		"""The maximum gray value that can be stored using
		the detector bit depth, assuming gray values are stored as
		unsigned integers.

		Returns
		-------
		max_gray_value : int
			2^bitdepth - 1
		"""
		return (2**self.get("bit_depth")) - 1

	def set_from_json(self, json_scenario:dict):
		"""Import the detector definition and geometry from the JSON object.
		The JSON object should contain the complete content
		of the scenario definition file
		(at least the geometry and detector sections).

		Parameters
		----------
		json_scenario : dict
			A complete CTSimU scenario object, as imported from a JSON structure.
		"""
		self.reset()

		# Extract the detector's geometry:
		geo = json_extract(json_scenario, ["geometry", "detector"])
		self.set_geometry(json_geometry_object=geo, proper_cs="local")

		# Detector properties:
		detprops = json_extract(json_scenario, ["detector"])
		Group.set_from_json(self, detprops)