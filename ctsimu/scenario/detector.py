# -*- coding: UTF-8 -*-
"""
A CTSimU detector: position, orientation, size, and other parameters.
"""

from ..helpers import *
from ..geometry import *
from .part import Part
from .parameter import Parameter

class Detector(Part):
	"""CTSimU Detector."""
	def __init__(self, name:str=""):
		"""A name can be passed when initializing the detector.

		Parameters
		----------
		name : str
			Detector name.
		"""
		Part.__init__(self, name)

		self.windows_front = list()
		self.windows_rear  = list()
		self.filters_front = list()
		self.filters_rear  = list()

		# Detector parameters
		self.set(key="model",            value="",      native_unit="string")
		self.set(key="manufacturer",     value="",      native_unit="string")
		self.set(key="type",             value="ideal", native_unit="string")
		self.set(key="columns",          value=1000,    native_unit=None)
		self.set(key="rows",             value=1000,    native_unit=None)
		self.set(key="pitch_u",          value=0.1,     native_unit="mm")
		self.set(key="pitch_v",          value=0.1,     native_unit="mm")
		self.set(key="bit_depth",        value=16,      native_unit=None)
		self.set(key="integration_time", value=1.0,     native_unit="s")
		self.set(key="dead_time",        value=1.0,     native_unit="s")
		self.set(key="image_lag",        value=0.0,     native_unit=None)
		self.set(key="frame_average",    value=1,       native_unit=None)
		self.set(key="multisampling",    value="3x3",   native_unit="string")

		# Properties for gray value reproduction:
		self.set(key="primary_energy_mode",             value=False, native_unit="bool")
		self.set(key="primary_intensity_mode",          value=False, native_unit="bool")
		self.set(key="imin",                            value=0,     native_unit=None)
		self.set(key="imax",                            value=60000, native_unit=None)
		self.set(key="factor",                          value=None,  native_unit=None)
		self.set(key="offset",                          value=None,  native_unit=None)
		self.set(key="gv_characteristics_file",         value="",    native_unit="string")
		self.set(key="efficiency_characteristics_file", value="",    native_unit="string")

		# Noise:
		self.set(key="snr_at_imax",                value=100,  native_unit=None)
		self.set(key="noise_characteristics_file", value=None, native_unit="string")

		# Unsharpness:
		self.set(key="basic_spatial_resolution", value=None, native_unit="mm")
		self.set(key="mtf_file",                 value=None, native_unit="string")
		self.set(key="long_range_unsharpness",   value=None, native_unit="mm")
		self.set(key="long_range_ratio",         value=None, native_unit=None)

		# Bad pixel map:
		self.set(key="bad_pixel_map",            value=None, native_unit="string")
		self.set(key="bad_pixel_map_dim_x",      value=None, native_unit=None)
		self.set(key="bad_pixel_map_dim_y",      value=None, native_unit=None)
		self.set(key="bad_pixel_map_type",       value=None, native_unit="string")
		self.set(key="bad_pixel_map_endian",     value=None, native_unit="string")
		self.set(key="bad_pixel_map_headersize", value=None, native_unit="string")

		# Scintillator
		self.set(key="scintillator_material_id", value=None, native_unit="string")
		self.set(key="scintillator_thickness",   value=None, native_unit="mm")

	def reset(self):
		Part.reset()
		self.windows_front = list()
		self.windows_rear = list()
		self.filters_front = list()
		self.filters_rear = list()

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

	def set_frame(self, frame:float, nFrames:int):
		"""Set all properties of the detector to match
		the given `frame` number, given a total of `nFrames`.
		"""

		# Update window and filter lists:
		for window in self.windows_front:
			window.set_frame(frame, nFrames, only_known_to_reconstruction=False)

		for window in self.windows_rear:
			window.set_frame(frame, nFrames, only_known_to_reconstruction=False)

		for filt in self.filters_front:
			filt.set_frame(frame, nFrames, only_known_to_reconstruction=False)

		for filt in self.filters_rear:
			filt.set_frame(frame, nFrames, only_known_to_reconstruction=False)

		Part.set_frame(frame, nFrames)

	def set_frame_for_reconstruction(self, frame:float, nFrames:int):
		"""Set all properties of the detector to match
		the given `frame` number, given a total of `nFrames`.
		"""

		# Update window and filter lists:
		for window in self.windows_front:
			window.set_frame(frame, nFrames, only_known_to_reconstruction=True)

		for window in self.windows_rear:
			window.set_frame(frame, nFrames, only_known_to_reconstruction=True)

		for filt in self.filters_front:
			filt.set_frame(frame, nFrames, only_known_to_reconstruction=True)

		for filt in self.filters_rear:
			filt.set_frame(frame, nFrames, only_known_to_reconstruction=True)

		Part.set_frame_for_reconstruction(frame, nFrames)

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

		# Extract the sample's geometry:
		geo = json_extract(json_scenario, ["geometry", "detector"])
		self.set_geometry(geo)

		# Detector properties:
		detprops = json_extract(json_scenario, ["detector"])

		self.set_parameter_from_key("model",        detprops, ["model"])
		self.set_parameter_from_key("manufacturer", detprops, ["manufacturer"])

		# Detector type:
		if self.set_parameter_from_key("type", detprops, ["type"], fail_value=None):
			# Check if the detector type is valid:
			if not (self.get("type") in valid_detector_types):
				raise ValueError(f"Not a valid detector type: \'{self.get('type')}\'. Should be any of {valid_detector_types}.")

		# Detector size:
		self.set_parameter_from_key("columns", detprops, ["columns"])
		self.set_parameter_from_key("row",     detprops, ["rows"])
		self.set_parameter_from_key("pitch_u", detprops, ["pitch_u"])
		self.set_parameter_from_key("pitch_v", detprops, ["pitch_v"])

		# More parameters:
		self.set_parameter_from_key("bit_depth",        detprops, ["bit_depth"])
		self.set_parameter_from_key("integration_time", detprops, ["integration_time"])
		self.set_parameter_from_key("dead_time",        detprops, ["dead_time"], fail_value=None)
		self.set_parameter_from_key("image_lag",        detprops, ["image_lag"], fail_value=None)

		# Gray value characteristics. (Support old British spelling, file format version < 1.0.)
		self.set_parameter_from_possible_keys("imin", detprops, [["gray_value", "imin"], ["grey_value", "imin"]], fail_value=None)
		self.set_parameter_from_possible_keys("imax", detprops, [["gray_value", "imax"], ["grey_value", "imax"]], fail_value=None)
		self.set_parameter_from_possible_keys("factor", detprops, [["gray_value", "factor"], ["grey_value", "factor"]], fail_value=None)
		self.set_parameter_from_possible_keys("offset", detprops, [["gray_value", "offset"], ["grey_value", "offset"]], fail_value=None)
		self.set_parameter_from_possible_keys("gv_characteristics_file", detprops, [["gray_value", "intensity_characteristics_file"], ["grey_value", "intensity_characteristics_file"]], fail_value=None)
		self.set_parameter_from_possible_keys("efficiency_characteristics_file", detprops, [["gray_value", "efficiency_characteristics_file"], ["grey_value", "efficiency_characteristics_file"]], fail_value=None)

		# Noise:
		self.set_parameter_from_key("snr_at_imax", detprops, ["noise", "snr_at_imax"], fail_value=None)
		self.set_parameter_from_key("noise_characteristics_file", detprops, ["noise", "noise_characteristics_file"], fail_value=None)

		# Unsharpness:
		self.set_parameter_from_possible_keys("basic_spatial_resolution", detprops, [["unsharpness", "basic_spatial_resolution"], ["sharpness", "basic_spatial_resolution"]], fail_value=None)
		self.set_parameter_from_possible_keys("mtf_file", detprops, [["unsharpness", "mtf"], ["sharpness", "mtf"]], fail_value=None)

		# Long range unsharpness is software-specific JSON parameter:
		self.set_parameter_from_key("long_range_unsharpness", detprops, ["simulation", "aRTist", "long_range_unsharpness", "extension"], fail_value=None)
		self.set_parameter_from_key("long_range_ratio",       detprops, ["simulation", "aRTist", "long_range_unsharpness", "ratio"], fail_value=None)

		# Bad pixel map:
		self.set_parameter_from_key("bad_pixel_map",    detprops, ["bad_pixel_map", "file"], fail_value=None)
		self.set_parameter_from_key("bad_pixel_map_dim_x", detprops, ["bad_pixel_map", "dim_x"], fail_value=None)
		self.set_parameter_from_key("bad_pixel_map_dim_y", detprops, ["bad_pixel_map", "dim_y"], fail_value=None)
		self.set_parameter_from_key("bad_pixel_map_type",  detprops, ["bad_pixel_map", "type"], fail_value=None)
		self.set_parameter_from_key("bad_pixel_map_type",  detprops, ["bad_pixel_map", "endian"], fail_value=None)
		self.set_parameter_from_key("bad_pixel_map_type",  detprops, ["bad_pixel_map", "headersize"], fail_value=None)

		# Scintillator:
		self.set_parameter_from_key("scintillator_material_id", detprops, ["scintillator", "material_id"], fail_value=None)
		self.set_parameter_from_key("scintillator_thickness",   detprops, ["scintillator", "thickness"], fail_value=None)

		# Windows and filters:
		self.windows_front = add_filters_to_list(self.windows_front, detprops, ["window", "front"])
		self.windows_rear  = add_filters_to_list(self.windows_rear,  detprops, ["window", "rear"])
		self.filters_front = add_filters_to_list(self.filters_front, detprops, ["filters", "front"])
		self.filters_rear  = add_filters_to_list(self.filters_rear,  detprops, ["filters", "rear"])

		# Frame averaging:
		self.set_parameter_from_key("frame_average", json_scenario, ["acquisition", "frame_average"], fail_value=None)

		# Software-specific: primary energy/intensity mode.
		self.set_parameter_from_key("primary_energy_mode", json_scenario, ["simulation", "aRTist", "primary_energies"], fail_value=None)
		self.set_parameter_from_key("primary_intensity_mode", json_scenario, ["simulation", "aRTist", "primary_intensities"], fail_value=None)