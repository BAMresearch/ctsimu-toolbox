# -*- coding: UTF-8 -*-
"""
Groups are collections of parameters.
"""

import copy
from ..helpers import *
from ..geometry import *
from .parameter import Parameter

class Group:
	"""A collection of parameters.

	A group is used to manage a collection of parameters (i.e., the group's
	properties) and set the frame number for all of those parameters at once.

	Attributes
	----------
	properties : dict
		A dictionary of part-specific properties.
		The dictionary elements are usually
		`ctsimu.scenario.parameter.Parameter` objects.
	"""
	def __init__(self, name:str="", _root=None):
		self._root = _root  # root scenario object

		self.name = name
		self.properties = dict() # Properties of class Parameter.
		self.subgroups = list()
		self.alternative_names = list()

	def reset(self):
		"""Reset all parameters in this group to their standard values."""
		for key in self.properties:
			self.properties[key].reset()

		for subgroup in self.subgroups:
			subgroup.reset()

	def json_dict(self) -> dict:
		"""Create a CTSimU JSON dictionary for this group and its subgroups.

		Returns
		-------
		json_dict : dict
		"""
		jd = dict()
		for key in self.properties:
			jd[key] = self.parameter(key).json_dict()

		for subgroup in self.subgroups:
			jd[subgroup.get_name()] = subgroup.json_dict()

		return jd

	def get_name(self):
		return self.name

	def set_name(self, name:str):
		self.name = name

	def add_alternative_name(self, alternative_name:str):
		if isinstance(alternative_name, str):
			self.alternative_names.append(alternative_name)
		else:
			raise TypeError("add_alternative_name: given alternative is not a string.")

	def add_subgroup(self, subgroup:'Group'):
		if isinstance(subgroup, Group):
			self.subgroups.append(subgroup)
		else:
			raise TypeError("add_subgroup: given subgroup is not an object of class 'Group'.")

	def parameter(self, key:str) -> 'Parameter':
		"""The parameter object behind a given property `key`.

		Parameters
		----------
		key : str
			Identifier for the requested property.

		Returns
		-------
		p : ctsimu.scenario.parameter.Parameter
			Parameter for the requested property.
		"""
		if key in self.properties:
			return self.properties[key]
		else:
			raise Exception(f"Part '{self.name}' does not have a property called '{key}'.")

	def get(self, key:str) -> float | str | bool:
		"""The current value for a given property `key`.

		Parameters
		----------
		key : str or list
			If key is a string: identifier for the local requested property.
			If key is a list of strings: key sequence to identify an
			item in a subgroup.

		Returns
		-------
		current_value : float or str or bool
			Current value of the requested property.
		"""
		if isinstance(key, str):
			return self.parameter(key).current_value
		elif isinstance(key, list):
			if len(key) == 1:
				# Only one item in the list: return parameter value.
				return self.parameter(key[0]).current_value
			elif len(key) == 0:
				raise Exception("Get: empty list does not identify a valid key.")
			else:
				next_in_line = key[0]
				# Identify the subgroup that is meant by the next keyword in line:
				for s in self.subgroups:
					if s.name == next_in_line:
						# Get the key element from the subgroup, remove
						# this 'next_in_line' from the front of the list:
						return s.get(key[1:])

		raise Exception("Get: please give a valid key string or list of strings")


	def standard_value(self, key:str) -> float | str | bool:
		"""The standard value for a given property `key`
		(i.e., the value unaffected by drifts).

		Parameters
		----------
		key : str
			Identifier for the requested property.

		Returns
		-------
		standard_value : float or str or bool
			Standard value of the requested property.
		"""
		return self.parameter(key).standard_value

	def changed(self, key:str) -> bool:
		"""Has the property changed its value since the last acknowledgment?
		The method `acknowledge_change` should be used to acknowledge
		a change of the value.

		Parameters
		----------
		key : str
			Identifier for the requested property.

		Returns
		-------
		changed : bool
			`True` if the requested property has changed since the
			last acknowledgment, `False` if not.
		"""
		return self.parameter(key).changed()

	def set(self, key:str, value, native_unit:str="undefined", simple:bool=False, valid_values:list=None):
		"""Set a property value in the group's properties dictionary.
		This function sets the respective parameter's standard value.
		If the parameter already exists in the internal properties dictionary,
		the parameter is reset (i.e., all its drifts are deleted) before the new
		standard value is set.

		Parameters
		----------
		key : str
			Identifier for the property to be set.

		value : float or str or bool
			New standard value for the property.

		native_unit : str, optional
			The native unit of the property.
			Possible values: `None`, `"mm"`, `"rad"`, `"deg"`, `"s"`, `"mA"`, `"kV"`, `"g/cm^3"`, `"lp/mm"`, `"bool"`, `"string"`.

		simple : bool
			Set to `True` if the parameter is represented
			as a single value in the JSON file instead of a
			parameter object.

		valid_values : list
			A list of acceptable standard values for the parameter.
		"""

		# Check if the property already exists:
		if key in self.properties:
			# Already exists in dict.
			# Reset parameter and set its standard value:
			self.properties[key].reset()
			if native_unit != "undefined":
				self.properties[key].set_native_unit(native_unit)

			self.properties[key].set_standard_value(value=value)
		else:
			# Property does not yet exist.
			if native_unit == "undefined":
				# If no native unit is specified, set it to `None`.
				native_unit = None

			# Create new parameter:
			new_parameter = Parameter(native_unit=native_unit, standard_value=value, simple=simple, valid_values=valid_values, _root=self._root)
			self.properties[key] = new_parameter

	def check(self):
		# Virtual function for consistency checks.
		# Needs specific implementations in children.
		return True

	def acknowledge_change(self, key:str, new_change_state:bool=False):
		"""Acknowledge a change of a property value due to a drift.
		After the acknowledgment, the method `changed`
		will return the `new_change_state` value (standard: `False`).

		Parameters
		----------
		new_change_state : bool, optional
			The new state to be returned by the method `changed`.
			Standard value: `False`.
		"""
		self.parameter(key).acknowledge_change(new_change_state=new_change_state)

	def set_parameter(self, key:str, parameter:'Parameter'):
		"""Set the internal property that is identified by the
		given property `key` to the given `parameter`.

		Parameters
		----------
		key : str
			Identifier for the property to be set.

		parameter : ctsimu.scenario.parameter.Parameter
			Parameter object for the given property `key`.
		"""

		self.properties[key] = parameter

	def set_parameter_value(self, key:str, dictionary:dict, key_sequence:list, native_unit:str=None, fail_value=None) -> bool:
		"""Set the value for the parameter that is identified
		by the property `key` in the internal properties dictionary.
		The new value is taken from the given JSON `dictionary`
		and located by the given `key_sequence`.

		If no value can be found at the given `key_sequence`,
		the `fail_value` is used instead.

		A `native_unit` can be provided in case the `property` does
		not yet exist in the internal properties dictionary.
		In this case, a new `ctsimu.scenario.parameter.Parameter`
		is created for the property and given the `native_unit`.

		If the parameter already exists in the internal properties dictionary,
		it is reset (i.e., all drifts are deleted).

		Parameters
		----------
		key : str
			Identifier for the property to be set.

		dictionary : dict
			CTSimU JSON structure for value, unit, drifts and uncertainty.

		key_sequence : list
			List of strings that identify the key sequence where the
			value is found in the given JSON structure.

		native_unit : str
			Native unit for the new parameter. Only necessary if the parameter
			for the given `key` does not yet exist.
			Possible values: `None`, `"mm"`, `"rad"`, `"deg"`, `"s"`, `"mA"`, `"kV"`, `"g/cm^3"`, `"lp/mm"`, `"bool"`, `"string"`.

		fail_value : float or str or bool, optional
			Value to be used if no value can be found in the given `dictionary`
			at the given `key_sequences`.

		Returns
		-------
		success : bool
			`True` if a valid value has been found to set up the parameter
			from the dictionary. `False` if no valid value has been
			found and possibly the `fail_value` had to be used instead.
		"""

		# Check if the property already exists:
		if key not in self.properties:
			# If not, create a new parameter object
			# and insert it into the properties dictionary:
			new_param = Parameter(native_unit=native_unit, standard_value=None, _root=self._root)
			self.properties[key] = new_param

		# Extract the value and set it as standard value:
		value = get_value(dictionary, key_sequence, fail_value="undefined")
		if value != "undefined":
			# Getting the value succeeded.
			self.set(key, value)
			return True
		else:
			# Getting the value failed or it is
			# not defined in the JSON file. Use
			# the fail_value instead and return `False`.
			self.set(key, fail_value)
			return False

		return False

	def set_parameter_from_key(self, key:str, dictionary:dict, key_sequence:list, fail_value="undefined", native_unit:str=None) -> bool:
		"""Set up a parameter object for the given
		property `key` from the `key_sequence` in the given `dictionary`.
		The object located at the key sequence must at least
		have a `value` property.

		If no value can be found at the given `key_sequence`,
		the `fail_value` is used instead.

		A `native_unit` can be provided in case the `property` does
		not yet exist in the internal properties dictionary.
		In this case, a new `ctsimu.scenario.parameter.Parameter`
		is created for the property and given the `native_unit`.

		If the parameter already exists in the internal properties dictionary,
		it is reset (i.e., all drifts are deleted).

		Parameters
		----------
		key : str
			Identifier for the property to be set.

		dictionary : dict
			CTSimU JSON structure for value, unit, drifts and uncertainty.

		key_sequence : list
			List of strings that identify the key sequence where the
			parameter is found in the given JSON structure.

		fail_value : float or str or bool, optional
			Value to be used if no value can be found in the given `dictionary`
			at the given `key_sequences`.

		native_unit : str
			Native unit for the new parameter. Only necessary if the parameter
			for the given `key` does not yet exist.
			Possible values: `None`, `"mm"`, `"rad"`, `"deg"`, `"s"`, `"mA"`, `"kV"`, `"g/cm^3"`, `"lp/mm"`, `"bool"`, `"string"`.

		Returns
		-------
		success : bool
			`True` if a valid value has been found to set up the parameter
			from the dictionary. `False` if no valid value has been
			found and possibly the `fail_value` had to be used instead.
		"""

		# Check if the property already exists:
		if key not in self.properties:
			# If not, create a new parameter object
			# and insert it into the properties dictionary:
			new_param = Parameter(native_unit=native_unit, standard_value=fail_value, _root=self._root)
			self.properties[key] = new_param

		# Try to set up the parameter:
		if not self.properties[key].set_parameter_from_key(dictionary, key_sequence):
			# Use fail value if this failed:
			if fail_value != "undefined":
				self.set(key, fail_value)

			return False

		return True

	def set_parameter_from_possible_keys(self, key:str, dictionary:dict, key_sequences:list, native_unit:str=None, fail_value="undefined") -> bool:
		"""Set up a parameter object for the given
		property `key` from the first matching key sequence in
		the `key_sequences` list. The object located at a
		key sequence must at least have a `value` property.

		If no value can be found at the given `key_sequences`,
		or no key sequence matches, the `fail_value` is used instead.

		A `native_unit` can be provided in case the `property` does
		not yet exist in the internal properties dictionary.
		In this case, a new `ctsimu.scenario.parameter.Parameter`
		is created for the property and given the `native_unit`.

		If the parameter already exists in the internal properties
		dictionary, it is reset (i.e., all drifts are deleted).

		Parameters
		----------
		key : str
			Identifier for the property to be set.

		dictionary : dict
			CTSimU JSON structure for value, unit, drifts and uncertainty.

		key_sequences : list
			List of lists of strings that identify the key sequences
			where the parameter may be found in the given JSON structure.

		native_unit : str
			Native unit for the new parameter. Only necessary if the parameter
			for the given `key` does not yet exist.
			Possible values: `None`, `"mm"`, `"rad"`, `"deg"`, `"s"`, `"mA"`, `"kV"`, `"g/cm^3"`, `"lp/mm"`, `"bool"`, `"string"`.

		fail_value : float or str or bool, optional
			Value to be used if no value can be found in the given `dictionary`
			at the given `key_sequences`.

		Returns
		-------
		success : bool
			`True` if a valid value has been found to set up the parameter
			from the dictionary. `False` if no valid value has been
			found and possibly the `fail_value` had to be used instead.
		"""

		# Check if the property already exists:
		if key not in self.properties:
			# If not, create a new parameter object
			# and insert it into the properties dictionary:
			new_param = Parameter(native_unit=native_unit, standard_value=fail_value, _root=self._root)
			self.properties[key] = new_param

		if not self.properties[key].set_parameter_from_possible_keys(dictionary, key_sequences):
			# Use fail value if this failed:
			if fail_value != "undefined":
				self.set(key, fail_value)

			return False

		return True

	def set_from_json(self, json_obj:dict):
		"""Set all properties of this group from the
		properties in the given `json_obj`. The group's
		property keys must match the keys in the `json_obj`.

		Very simple implementation that only works for
		trivial, linear (i.e., non-nested) parameter groups.
		Usually, this function is meant to be replaced by
		the child which inherits from the `Group` class.

		Parameters
		----------
		json_obj : dict
			A group of properties from a CTSimU JSON scenario
			that matches this group's properties.
		"""

		if json_obj is None:
			self.reset()
			return

		for key in self.properties:
			try:
				self.set_parameter_from_key(key, json_obj, [key], fail_value=self.parameter(key).fail_value)
			except Exception as e:
				raise Exception(f"Error setting key '{key}' of '{self.name}': {e}")

		for subgroup in self.subgroups:
			json_subobj = None

			# Check if subgroup's name is in json_obj,
			# extract if found:
			if subgroup.get_name() in json_obj:
				json_subobj = json_extract(json_obj, [subgroup.get_name()])
			else:
				# If subgroup is not found, try alternative names:
				for alternative_name in subgroup.alternative_names:
					if alternative_name in json_obj:
						json_subobj = json_extract(json_obj, [alternative_name])
						break

			subgroup.set_from_json(json_subobj)

	def set_frame(self, frame:float, nFrames:int, reconstruction:bool=False):
		"""
		Set up the group for the given `frame` number, obeying all drifts.

		Parameters
		----------
		frame : float
			Current frame number.

		nFrames : int
			Total number of frames in scan.

		reconstruction : bool
			If `True`, set frame as seen by reconstruction software.
			Default: `False`.
		"""

		# Set the frame for all elements of the properties dictionary:
		for key in self.properties:
			self.properties[key].set_frame(
				frame=frame,
				nFrames=nFrames,
				reconstruction=reconstruction
			)

		for subgroup in self.subgroups:
			subgroup.set_frame(frame, nFrames, reconstruction)

class Array(Group):
	"""Array of objects that are all of the same kind, such as
	the CTSimU materials list."""
	def __init__(self, name:str="", _root=None):
		self._root = _root  # root scenario object

		Group.__init__(self, name, _root)
		self.elements = []

	def reset(self):
		for element in self.elements:
			element.reset()

	def json_dict(self) -> list:
		jd = list()
		for element in self.elements:
			jd.append(element.json_dict())

		return jd

	def get(self, i:int):
		if isinstance(i, int):
			if i < len(self.elements):
				return self.elements[i]
		elif isinstance(i, list):
			pos = i[0]
			if isinstance(pos, int):
				return self.elements[pos].get(i[1:])
			else:
				raise Exception(f"'{pos}' does not identify a position in the array '{self.name}'. Please provide an integer.")

		raise Exception(f"Array '{self.name}': cannot get element at key sequence {i}.")

	def n_elements(self):
		return len(self.elements)

	def add_element(self, json_obj:dict):
		# Create new group:
		new_group = Group()
		new_group.properties = copy.deepcopy(self.properties)
		new_group.subgroups  = copy.deepcopy(self.subgroups)
		new_group.set_from_json(json_obj)

		self.elements.append(new_group)

	def set_from_json(self, json_obj:list):
		if json_obj is None:
			self.elements = list()
			self.reset()
			return

		if isinstance(json_obj, list):
			for element in json_obj:
				self.add_element(element)
		else:
			# Given object is not a list. Try to interpret as single object.
			if isinstance(json_obj, dict):
				self.add_element(json_obj)
			else:
				raise TypeError("Array: set_from_json: given object is neither a list nor a single JSON dictionary.")

	def set_frame(self, frame:float, nFrames:int, reconstruction:bool=False):
		"""
		Set up the group for the given `frame` number.

		Parameters
		----------
		frame : float
			Current frame number.

		nFrames : int
			Total number of frames in scan.

		reconstruction : bool
			If `True`, set frame as seen by reconstruction software.
			Default: `False`.
		"""

		for element in self.elements:
			element.set_frame(frame, nFrames, reconstruction)