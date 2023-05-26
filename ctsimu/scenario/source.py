# -*- coding: UTF-8 -*-
"""
A CTSimU X-ray source: position, orientation, size, and other parameters.
"""

from ..helpers import *
from ..geometry import *
from .part import Part
from .parameter import Parameter

class Source(Part):
	"""CTSimU X-ray source."""
	def __init__(self, name:str=""):
		"""A name can be passed when initializing the X-ray source.

		Parameters
		----------
		name : str
			Detector name.
		"""
		Part.__init__(self, name)

		self.windows = list()
		self.filters  = list()

		# X-ray source parameters
		self.set(key="model",        value="",  native_unit="string")
		self.set(key="manufacturer", value="",  native_unit="string")
		self.set(key="voltage",      value=130, native_unit="kV")
		self.set(key="current",      value=0.1, native_unit="mA")

		# Target
		self.set(key="target_material_id",     value="W",          native_unit="string")
		self.set(key="target_type",            value="reflection", native_unit="string")
		self.set(key="target_thickness",       value=10,           native_unit="mm")
		self.set(key="target_angle_incidence", value=45,           native_unit="deg")
		self.set(key="target_angle_emission",  value=45,           native_unit="deg")

		# Spot
		self.set(key="spot_size_u",        value=0,    native_unit="mm")
		self.set(key="spot_size_v",        value=0,    native_unit="mm")
		self.set(key="spot_size_w",        value=0,    native_unit="mm")
		self.set(key="spot_sigma_u",       value=0,    native_unit="mm")
		self.set(key="spot_sigma_v",       value=0,    native_unit="mm")
		self.set(key="spot_sigma_w",       value=0,    native_unit="mm")
		self.set(key="spot_multisampling", value="20", native_unit="string")

		# Intensity map
		self.set(key="intensity_map",            value=None, native_unit="string")
		self.set(key="intensity_map_dim_x",      value=None, native_unit=None)
		self.set(key="intensity_map_dim_y",      value=None, native_unit=None)
		self.set(key="intensity_map_dim_z",      value=None, native_unit=None)
		self.set(key="intensity_map_type",       value=None, native_unit="string")
		self.set(key="intensity_map_endian",     value=None, native_unit="string")
		self.set(key="intensity_map_headersize", value=None, native_unit="string")

		# Spectrum
		self.set(key="spectrum_monochromatic", value=False, native_unit="bool")
		self.set(key="spectrum_file",          value="",    native_unit="string")
		self.set(key="spectrum_resolution",    value=1.0,   native_unit="keV")


	def reset(self):
		Part.reset()
		self.windows = list()
		self.filters = list()

	def set_frame(self, frame:float, nFrames:int, w_rotation:float=0):
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

		Part.set_frame(frame, nFrames, w_rotation)

	def set_frame_for_reconstruction(self, frame:float, nFrames:int, w_rotation:float=0):
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

		Part.set_frame_for_reconstruction(frame, nFrames, w_rotation)

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
