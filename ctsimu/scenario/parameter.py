# -*- coding: UTF-8 -*-

import numbers
from ..helpers import *
from .drift import Drift

class Parameter:
	"""Parameter value, includes handling of parameter drifts.

	Attributes
	----------
	standard_value : float or str or bool
		The parameter's standard value, i.e. its value without any drifts.

	current_value : float or str or bool
		The parameter's current value for the given frame number,
		as set in `Parameter.set_frame()`.

	native_unit : str
		The parameter's native unit.
		Possible values: `None`, `"mm"`, `"rad"`, `"deg"`, `"s"`, `"mA"`, `"kV"`, `"g/cm^3"`, `"lp/mm"`, `"bool"`, `"string"`.

	drifts : list
		A list of subsequent parameter drifts, all of them
		are objects of class `ctsimu.scenario.drift.Drift`.

	value_has_changed : bool
		`True` if the parameter's value has changed since the last value
		has been acknowledged. A parameter value can be acknowledged by
		calling the method `Parameter.acknowledge_change()`.
	"""

	def __init__(self, native_unit:str=None, standard_value=0):
		"""A new parameter object must be assigned a valid native unit
		to enable the JSON parser to convert the drift values from the
		JSON file, if necessary.

		Optionally, a standard value can be passed to the constructor.
		The standard value is the "actual" value defined for this
		parameter in the JSON file. If a JSON object is used to set up
		this parameter, the standard value provided in the constructor
		is overwritten by the value given in the JSON file.

		Parameters
		----------
		standard_value : float or str or bool
			The parameter's standard value, i.e. its value without any drifts.

		native_unit : str
			The parameter's native unit.
			Possible values: `None`, `"mm"`, `"rad"`, `"deg"`, `"s"`, `"mA"`, `"kV"`, `"g/cm^3"`, `"lp/mm"`, `"bool"`, `"string"`.
		"""
		self.standard_value = standard_value
		self.native_unit = native_unit
		self.drifts = []
		self.current_value = standard_value
		self.value_has_changed = True

		self.reset()

	def __str__(self):
		s  = f"Standard value: {self.standard_value}, "
		s += f"Current value: {self.current_value}, "
		s += f"Native unit: {self.native_unit}, "
		s += f"nDrifts: {len(self.drifts)}"

		return s

	def reset(self):
		"""Delete all drifts and set the parameter's current value
		to the standard value.
		"""

		self.drifts = []
		self.current_value = self.standard_value
		self.value_has_changed = True

	def has_changed(self) -> bool:
		"""Has the parameter changed since the last acknowledged change?
		(See function `acknowledge_change`).

		Returns
		-------
		value_has_changed : bool
			`True` if the value has changed since last acknowledgment,
			otherwise `False`.
		"""
		return value_has_changed

	def has_drifts(self) -> bool:
		"""Does the parameter specify any drift?

		Returns
		-------
		has_drifts : bool
			`True` if the parameter specifies any drift, otherwise `False`.
		"""
		return (len(self.drifts) > 0)

	def maximum_value(self, nFrames:int, only_drifts_known_to_reconstruction:bool=False) -> float:
		"""Get the parameter's maximum value during the evolution of `nFrames`, given drifts.

		Parameters
		----------
		nFrames : int
			Number of frames in the CT scan.

		only_drifts_known_to_reconstruction : bool
			Obey only those drifts that are labeled as known to
			the reconstruction software?

		Returns
		-------
		maximum_value : float
			Maximum value of this parameter, considering the parameter
			drifts throughout the CT scan.
		"""
		if self.has_drifts():
			if (self.native_unit != "string") and (self.native_unit != "bool") and (self.standard_value != None):
				total_drift_max = self.get_total_drift_value_for_frame(
					0,
					nFrames,
					only_drifts_known_to_reconstruction
				)

				for f in range(nFrames):
					total_drift_for_frame = self.get_total_drift_value_for_frame(
						f,
						nFrames,
						only_drifts_known_to_reconstruction
					)

					if total_drift_for_frame > total_drift_max:
						total_drift_max = total_drift_for_frame

				return self.standard_value + total_drift_max

		return standard_value

	def minimum_value(self, nFrames:int, only_drifts_known_to_reconstruction:bool=False) -> float:
		"""Get the parameter's minimum value during the evolution of `nFrames`, given drifts.

		Parameters
		----------
		nFrames : int
			Number of frames in the CT scan.

		only_drifts_known_to_reconstruction : bool
			Obey only those drifts that are labeled as known to
			the reconstruction software?

		Returns
		-------
		minimum_value : float
			Minimum value of this parameter, considering the parameter
			drifts throughout the CT scan.
		"""
		if self.has_drifts():
			if (self.native_unit != "string") and (self.native_unit != "bool") and (self.standard_value != None):
				total_drift_min = self.get_total_drift_value_for_frame(
					0,
					nFrames,
					only_drifts_known_to_reconstruction
				)

				for f in range(nFrames):
					total_drift_for_frame = self.get_total_drift_value_for_frame(
						f,
						nFrames,
						only_drifts_known_to_reconstruction
					)

					if total_drift_for_frame < total_drift_min:
						total_drift_min = total_drift_for_frame

				return self.standard_value + total_drift_min

		return standard_value

	def set_native_unit(self, native_unit:str):
		"""Set the parameter's native unit.

		Parameters
		----------
		native_unit : str
			New native unit for the parameter.
			Possible values: `None`, `"mm"`, `"rad"`, `"deg"`, `"s"`, `"mA"`, `"kV"`, `"g/cm^3"`, `"lp/mm"`, `"bool"`, `"string"`.
		"""
		if is_valid_native_unit(native_unit):
			self.native_unit = native_unit

	def set_standard_value(self, value):
		"""Set the parameter's standard value.
		Automatically sets the current value to the standard value.

		Parameters
		----------
		value : float or str or bool
			New standard value for the parameter.
		"""
		self.standard_value = value
		self.current_value = value

	def acknowledge_change(self, new_change_state:bool=False):
		"""Acknowledge a change of the parameter due to a drift.
		After the acknowledgment, the function `has_changed`
		will return the `new_change_state` value (standard: `False`).

		Parameters
		----------
		new_change_state : bool
			The new state to be returned by the method `has_changed`.
			Standard value: `False`.
		"""
		self.value_has_changed = new_change_state

	def set_frame_and_get_value(self, frame:float, nFrames:int=1, only_drifts_known_to_reconstruction:bool=False) -> float | str | bool:
		"""Set the new frame number, return the new current value.

		Parameters
		----------
		frame : float
			Current frame number.

		nFrames : int
			Total number of frames in scan.

		only_drifts_known_to_reconstruction : bool
			Obey only those drifts that are labeled as known to the
			reconstruction software?

		Returns
		-------
		current_value : float or str or bool
			The parameter's current value for the given `frame` number.
		"""
		self.set_frame(frame, nFrames, only_drifts_known_to_reconstruction)
		return self.current_value

	def get_total_drift_value_for_frame(self, frame:float, nFrames:int, only_drifts_known_to_reconstruction:bool=False) -> float | str | bool:
		"""Calculate the total drift value from all drift components,
		for the given `frame` out of a total of `nFrames`,
		depending on whether all drifts are applied or only
		drifts known to the reconstruction software.

		Parameters
		----------
		frame : float
			Current frame number.

		nFrames : int
			Total number of frames in CT scan.

		only_drifts_known_to_reconstruction : bool
			Obey only those drifts that are labeled as known to the
			reconstruction software?

		Returns
		-------
		total_drift_value : float or str or bool
			The parameter's total (accumulated) drift value
			for the given `frame` number.
		"""
		total_drift = 0

		for d in self.drifts:
			if only_drifts_known_to_reconstruction is True:
				if d.known_to_reconstruction is False:
					# Skip this drift if it is unknown to the reconstruction,
					# but we only want to obey drifts that are actually known
					# to the reconstruction...
					continue

			if (self.native_unit == "string") or (self.native_unit == "bool"):
				# A string-type parameter can only be one string,
				# nothing is added, and the drifts array should only
				# contain one element. Otherwise, the last drift is the
				# one that has precedence. Same for bool.
				total_drift = d.get_value_for_frame(frame=frame, nFrames=nFrames)
			else:
				# Add up all drift values for requested frame:
				total_drift += d.get_value_for_frame(frame=frame, nFrames=nFrames)

		return total_drift

	def add_drift_json(self, json_drift_object:dict):
		"""Generate a `ctsimu.scenario.drift.Drift` object from a JSON structure
		that defines a drift and add it to the parameter's internal list of
		drifts to handle.

		Parameters
		----------
		json_drift_object : dict
			A CTSimU drift structure, imported from a JSON structure.
		"""
		d = Drift(native_unit=self.native_unit)
		d.set_from_json(json_drift_object)
		self.drifts.append(d)

	def set_from_json(self, json_parameter_object:dict) -> bool:
		"""Set up this parameter from a CTSimU parameter structure.
		The proper `native_unit` must be set up correctly before
		running this function.

		Parameters
		----------
		json_parameter_object : dict
			A CTSimU parameter structure, imported from a JSON structure.

		Returns
		-------
		success : bool
			`True` on success, `False` if an error occurred.
		"""
		self.reset()
		success = False

		if isinstance(json_parameter_object, numbers.Number):
			# Parameter is given as a single number, not as a
			# parameter structure (with value, unit, drift, uncertainty)
			if self.native_unit != "string":
				self.set_standard_value(get_value(json_parameter_object))
				success = True
		elif isinstance(json_parameter_object, str):
			# Parameter is given as a single string, not as a
			# parameter structure (with value, unit, drift, uncertainty)
			if self.native_unit == "string":
				self.set_standard_value(get_value(json_parameter_object))
				success = True
		elif isinstance(json_parameter_object, bool):
			# Parameter is given as a boolean.
			if self.native_unit == "bool":
				self.set_standard_value(json_convert_to_native_unit(self.native_unit, json_parameter_object))
				success = True
		elif isinstance(json_parameter_object, dict):
			# Hopefully a valid parameter structure...

			# Value, automatically converted to parameter's native unit:
			if json_exists_and_not_null(json_parameter_object, ["value"]):
				if not object_value_is_null(json_parameter_object):
					self.set_standard_value(json_convert_to_native_unit(self.native_unit, json_parameter_object))
					success = True

			self.current_value = self.standard_value

			# Drifts:
			if json_exists_and_not_null(json_parameter_object, ["drifts"]):
				drifts = json_extract(json_parameter_object, ["drifts"])

				if isinstance(drifts, list):
					# An array of drift structures.
					for d in drifts:
						self.add_drift(d)
				elif isinstance(drifts, dict):
					# A single drift structure (apparently).
					# This is not standard behavior, but let's support it.
					self.add_drift(drifts)

		return success

	def set_parameter_from_key(self, json_object:dict, key_sequence:list) -> bool:
		"""Tries to find a valid parameter object at the given
		`key_sequence` in the given `json_object`. Sets the parameter
		if possible and returns `True` on success, `False` otherwise.

		Parameters
		----------
		json_object : dict
			Piece of a CTSimU JSON structure.

		key_sequence : list
			List of strings that identify the key sequence where the
			parameter is found in the given JSON structure.

		Returns
		-------
		success : bool
			`True` on success, `False` if an error occurred.
		"""
		if json_exists_and_not_null(json_object, key_sequence):
			if self.set_from_json(json_extract(json_object, key_sequence)):
				return True

		return False

	def set_parameter_from_possible_keys(self, json_object:dict, key_sequences:list) -> bool:
		"""Searches the JSON object for each key sequence in the given list
		of `key_sequences`.	The first sequence that can be found is taken to
		set up the parameter.

		Parameters
		----------
		json_object : dict
			Piece of a CTSimU JSON structure.

		key_sequences : list
			A list of possible key sequences, i.e., a list of lists of strings.

		Returns
		-------
		success : bool
			`True` on success (i.e., at least one of the key sequences
			identified a valid parameter in the given JSON structure),
			`False` if an error occurred.
		"""
		for keyseq in key_sequences:
			if self.set_parameter_from_key(json_object, keyseq):
				# Returned successfully, so we are done here.
				return True

		return False

	def set_frame(self, frame:float, nFrames:int, only_drifts_known_to_reconstruction:bool=False) -> bool:
		"""Prepares the `current_value` for the given `frame` number
		(assuming a total of `nFrames`). This takes into account all drifts
		(or only the ones known to reconstruction).

		Parameters
		----------
		frame : float
			Current frame number.

		nFrames : int
			Total number of frames in scan.

		only_drifts_known_to_reconstruction : bool
			Obey only those drifts that are labeled as known to the
			reconstruction software?

		Returns
		-------
		value_has_changed : bool
			If the value has changed from the previous value (e.g. due to drifts.)
		"""
		new_value = self.standard_value
		total_drift = self.get_total_drift_value_for_frame(frame, nFrames, only_drifts_known_to_reconstruction)

		if self.native_unit == "string":
			if isinstance(total_drift, str):
				# Only assign new value if `total_drift` is actually a string:
				new_value = total_drift
		else:
			if self.standard_value is not None:
				new_value = self.standard_value + total_drift

		# Check if the value has changed when compared to the previous value:
		value_has_changed = False
		if self.current_value != new_value:
			value_has_changed = True
			self.current_value = new_value

		return value_has_changed