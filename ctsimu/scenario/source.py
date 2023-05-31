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

		# X-ray source parameters
		self.set(key="model",        value=None, native_unit="string", simple=True)
		self.set(key="manufacturer", value=None, native_unit="string", simple=True)
		self.set(key="voltage",      value=None, native_unit="kV")
		self.set(key="current",      value=None, native_unit="mA")

		# Target
		self.target = Group("target")
		self.target.set(key="material_id",     value=None,         native_unit="string", simple=True)
		self.target.set(key="type",            value="reflection", native_unit="string", simple=True)
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
		self.spectrum.set(key="monochromatic", value=False, native_unit="bool")
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

	def set_frame_legacy(self, frame:float, nFrames:int, w_rotation:float=0):
		"""Set all properties of the X-ray source to match
		the given `frame` number, given a total of `nFrames`.

		All drifts and deviations are obeyed.

		Parameters
		----------
		frame : float
			Current frame number.

		nFrames : int
			Total number of frames in scan.

		w_rotation : float
			An additional rotation (in rad) around the stage's w axis
			for this frame. Used for the sample stage, which
			rotates during a CT scan.
		"""

		# Update window and filter lists:
		for window in self.windows:
			window.set_frame(frame, nFrames, only_known_to_reconstruction=False)

		for filt in self.filters:
			filt.set_frame(frame, nFrames, only_known_to_reconstruction=False)

		Part.set_frame(self, frame, nFrames, w_rotation)

	def set_frame_for_reconstruction_legacy(self, frame:float, nFrames:int, w_rotation:float=0):
		"""Set all properties of the X-ray source to match
		the given `frame` number, given a total of `nFrames`.

		Only those drifts and deviations are obeyed which are known to the
		reconstruction software.

		Parameters
		----------
		frame : float
			Current frame number.

		nFrames : int
			Total number of frames in scan.

		w_rotation : float
			An additional rotation (in rad) around the stage's w axis
			for this frame. Used for the sample stage, which
			rotates during a CT scan.
		"""

		# Update window and filter lists:
		for window in self.windows:
			window.set_frame(frame, nFrames, only_known_to_reconstruction=True)

		for filt in self.filters:
			filt.set_frame(frame, nFrames, only_known_to_reconstruction=True)

		Part.set_frame_for_reconstruction(self, frame, nFrames, w_rotation)

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
		self.set_geometry(geo)

		# Source properties:
		sourceprops = json_extract(json_scenario, ["source"])
		Group.set_from_json(self, sourceprops)

		"""
		# Source properties:
		sourceprops = json_extract(json_scenario, ["source"])

		self.set_parameter_from_key("model",        sourceprops, ["model"])
		self.set_parameter_from_key("manufacturer", sourceprops, ["manufacturer"])

		self.set_parameter_from_key("voltage", sourceprops, ["voltage"])
		self.set_parameter_from_key("current", sourceprops, ["current"])

		# Target:
		self.set_parameter_from_key("target_material_id",     sourceprops, ["target", "material_id"])
		self.set_parameter_from_key("target_thickness",       sourceprops, ["target", "thickness"])
		self.set_parameter_from_key("target_angle_incidence", sourceprops, ["target", "angle", "incidence"])
		self.set_parameter_from_key("target_angle_emission",  sourceprops, ["target", "angle", "emission"])
		if self.set_parameter_from_key("target_type", sourceprops, ["target", "type"], fail_value=None):
			# Check if the target type is valid:
			if not (self.get("target_type") in valid_xray_target_types):
				raise ValueError(f"Not a valid X-ray source target type: \'{self.get('target_type')}\'. Should be any of {valid_xray_target_types}.")

		# Spot:
		self.set_parameter_from_key("spot_size_u",  sourceprops, ["spot", "size", "u"])
		self.set_parameter_from_key("spot_size_v",  sourceprops, ["spot", "size", "v"])
		self.set_parameter_from_key("spot_size_w",  sourceprops, ["spot", "size", "w"])
		self.set_parameter_from_key("spot_sigma_u", sourceprops, ["spot", "sigma", "u"])
		self.set_parameter_from_key("spot_sigma_v", sourceprops, ["spot", "sigma", "v"])
		self.set_parameter_from_key("spot_sigma_w", sourceprops, ["spot", "sigma", "w"])

		# Intensity map:
		self.set_parameter_from_key("intensity_map",            sourceprops, ["spot", "intensity_map", "file"])
		self.set_parameter_from_key("intensity_map_type",       sourceprops, ["spot", "intensity_map", "type"])
		self.set_parameter_from_key("intensity_map_dim_x",      sourceprops, ["spot", "intensity_map", "dim_x"])
		self.set_parameter_from_key("intensity_map_dim_y",      sourceprops, ["spot", "intensity_map", "dim_y"])
		self.set_parameter_from_key("intensity_map_dim_z",      sourceprops, ["spot", "intensity_map", "dim_z"])
		self.set_parameter_from_key("intensity_map_endian",     sourceprops, ["spot", "intensity_map", "endian"])
		self.set_parameter_from_key("intensity_map_headersize", sourceprops, ["spot", "intensity_map", "headersize"])

		# Spectrum:
		self.set_parameter_from_key("spectrum_monochromatic", sourceprops, ["spectrum", "monochromatic"])
		self.set_parameter_from_key("spectrum_file",          sourceprops, ["spectrum", "file"])
		self.set_parameter_from_key("spectrum_resolution",    json_scenario, ["simulation", "aRTist", "spectral_resolution"])

		# Windows and filters:
		self.windows = add_filters_to_list(self.windows, detprops, ["window"])
		self.filters = add_filters_to_list(self.filters, detprops, ["filters"])
		"""