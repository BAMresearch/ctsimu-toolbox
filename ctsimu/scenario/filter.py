# -*- coding: UTF-8 -*-
"""
A filter or window for an X-ray tube or detector.
"""
import numbers

from ..helpers import *
from .parameter import Parameter
from .drift import Drift

class Filter:
	"""A filter or window for an X-ray tube or detector.

	Attributes
	----------
	material_id : str
		Filter material, referencing the ID of a material
		from the scenario description.

	thickness : ctsimu.scenario.parameter.Parameter
		Filter thickness (in mm).
	"""
	def __init__(self, material_id:str="Fe", thickness:float=1.0):
		"""Both attributes can be set when the filter object
		is initialized.

		Parameters
		----------
		material_id : str
			Filter material, referencing the ID of a material
			from the scenario description.

		thickness : float or ctsimu.scenario.parameter.Parameter
			Filter thickness (in mm).
		"""
		self.material_id = material_id
		self.thickness = None
		self.set_thickness(thickness)

	def set_thickness(self, thickness):
		"""Set the filter thickness.

		Parameters
		----------
		thickness : float or ctsimu.scenario.parameter.Parameter
			The filter thickness (in mm). Given as a simple numerical
			value (`float`), it will be converted to a Parameter
			with this as the standard value.
			It can also be given directly as a `Parameter` object.

		Raises
		------
		TypeError
			If the `thickness` is neither of type `float` nor an instance
			of `ctsimu.scenario.parameter.Parameter`.
		"""
		if isinstance(thickness, numbers.Number):
			self.thickness = Parameter(native_unit="mm", standard_value=thickness)
		elif isinstance(thickness, Parameter):
			self.thickness = thickness
		else:
			raise TypeError(f"CTSimU Filter: set_thickness: invalid type of given thickness: '{type(density)}'. Must be 'float' or 'Parameter'.")

	def set_frame(self, frame:float, nFrames:int, reconstruction:bool=False) -> bool:
		"""Prepare the thickness parameter for the given `frame`
		number (out of a total of `nFrames`). Handles possible drifts.

		Parameters
		----------
		frame : float
			Current frame number.

		nFrames : int
			Total number of frames in scan.

		reconstruction : bool
			`True` if only those drifts known to the reconstruction software
			must be taken into account, `False` if all drifts are taken
			into account.

		Returns
		-------
		changed : bool
			`True` if the thickness has changed since the
			previous state, `False` if not.
		"""

		changed = self.thickness.set_frame(frame, nFrames, reconstruction)
		return changed

	def set_from_json(self, json_object:dict):
		"""Sets the filter properties from a given CTSimU JSON filter object.

		Parameters
		----------
		json_object : dict
			CTSimU filter definition, as imported from a JSON file.

		Raises
		------
		ValueError
			If `material_id` or `thickness` cannot be read from the `json_object`.
		"""

		if json_exists_and_not_null(json_object, ["material_id"]):
			self.material_id = get_value(json_object, ["material_id"])
		else:
			raise ValueError(f"CTSimU Filter: set_from_json: Error reading a filter's 'material_id'. Value wrong or not specified.")

		if not self.thickness.set_parameter_from_key(json_object, ["thickness"]):
			raise ValueError(f"CTSimU Filter: set_from_json: Error reading a filter's 'thickness'. Value wrong or not specified.")

		self.set_frame(0, 1, False)