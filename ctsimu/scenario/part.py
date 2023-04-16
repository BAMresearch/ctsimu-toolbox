# -*- coding: UTF-8 -*-
"""
Parts are objects in the scene: detector, source, stage and samples.
"""

from ..helpers import *
from ..geometry import *
from .deviation import Deviation
from .parameter import Parameter
from .scenevector import Scenevector

class Part:
	"""
	Parts are objects in the scene: detector, source, stage and samples.

	They have a coordinate system and can define deviations from their
	standard geometries (translations and rotations around given axes).
	The center, vectors and deviations can all have drifts,
	allowing for an evolution in time.

	Attributes
	----------
	name : str
		Name of the object.

	properties : dict
		A dictionary of part-specific properties.
		The dictionary elements are usually
		`ctsimu.scenario.parameter.Parameter` objects.

	attached_to_stage : bool
		Is the part attached to the sample stage?
		If so, it will follow the stage's movements throughout the CT scan.

	coordinate_system : ctsimu.geometry.CoordinateSystem
		The current coordinate system. This variable is automatically
		set up whenever the `set_frame` method is called.

	center : ctsimu.scenario.scenevector.Scenevector
		Coordinates of the center point.

	u : ctsimu.scenario.scenevector.Scenevector
		u vector of the part's coordinate system.

	w : ctsimu.scenario.scenevector.Scenevector
		w vector of the part's coordinate system.

	deviations : list
		All the geometric deviations for this part as
		`ctsimu.scenario.deviation.Deviation` objects.

	legacy_deviations : list
		All the legacy deviations for this part as
		`ctsimu.scenario.deviation.Deviation` objects.
		Legacy deviations are geometric deviations from
		CTSimU scenario descriptions prior to
		file format version 1.0. In general, they are incompatible
		with the new geometric deviation format and must be
		treated differently.
	"""
	def __init__(self, name:str="Part"):
		"""The part's name can be set upon initialization.

		Parameters
		----------
		name : str
			Name for the part.
		"""
		self.name = name
		self.properties = dict() # Properties of class Parameter.

		self.attached_to_stage = False
		self.coordinate_system = CoordinateSystem()
		self.center = SceneVector(native_unit="mm")
		self.u = SceneVector(native_unit=None)
		self.w = SceneVector(native_unit=None)
		self.deviations = list()
		self.legacy_deviations = list() # for deviations prior to file format 1.0

		# Internal parameters used for efficiency and
		# to track current state.
		self._static = False # non-drifting object if `True`.

		# If certain geometry deviations or drifts are unknown to the
		# reconstruction software, the part may be set up in two different ways:
		self._cs_initialized_real = False  # initialized to real coordinate system?
		self._cs_initialized_recon = False # initialized to recon coordinate system?

	def reset(self):
		self.attached_to_stage = False
		self._static = False
		self._cs_initialized_real = False
		self._cs_initialized_recon = False

		self.coordinate_system.reset()

		self.deviations = list()
		self.legacy_deviations = list()

		for key in self.properties:
			self.properties[key].reset()

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
			raise Exception(f"Part '{self.name}'' does not have a property called '{key}'.")

	def get(self, key:str) -> float | str | bool:
		"""The current value for a given property `key`.

		Parameters
		----------
		key : str
			Identifier for the requested property.

		Returns
		-------
		current_value : float or str or bool
			Current value of the requested property.
		"""
		return self.parameter(key).current_value()

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
		return self.parameter(key).standard_value()

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

	def is_attached_to_stage(self) -> bool:
		"""Is the part attached to the sample stage?

		Returns
		-------
		attached_to_stage : bool
			`True` is part is attached to sample stage, i.e., its reference
			coordinate system is the stage coordinate system.
			`False` if the part is not attached to the sample stage, i.e.,
			its reference coordinate system is the world coordinate system.
		"""
		return self.attached_to_stage

	def set(self, key:str, value, native_unit:str="undefined"):
		"""Set a property value in the part's properties dictionary.
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
			new_parameter = Parameter(native_unit=native_unit, standard_value=value)
			self.properties[key] = new_parameter

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
			parameter is found in the given JSON structure.

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
			new_param = Parameter(native_unit=native_unit, standard_value=None)
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

	def set_parameter_from_key(self, key:str, dictionary:dict, key_sequence:list, native_unit:str=None, fail_value="undefined") -> bool:
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
			new_param = Parameter(native_unit=native_unit, standard_value=fail_value)
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
			new_param = Parameter(native_unit=native_unit, standard_value=fail_value)
			self.properties[key] = new_param

		if not self.properties[key].set_parameter_from_possible_keys(dictionary, key_sequences):
			# Use fail value if this failed:
			if fail_value != "undefined":
				self.set(key, fail_value)

			return False

		return True

	def set_name(self, name:str):
		"""Set the part's name.

		Parameters
		----------
		name : str
			New name for the part.
		"""
		self.name = name

	def set_center(self, center:'Scenevector'):
		"""Set the part's center.

		Parameters
		----------
		center : Scenevector
			New center for the part.
		"""
		self.center = center

	def set_u(self, u:'Scenevector'):
		"""Set the part's `u` vector.

		Parameters
		----------
		u : Scenevector
			New `u` vector for the part.
		"""
		self.u = u

	def set_w(self, w:'Scenevector'):
		"""Set the part's `w` vector.

		Parameters
		----------
		w : Scenevector
			New `w` vector for the part.
		"""
		self.w = w

	def attach_to_stage(self, attached:bool=True):
		"""Set the part's `attached_to_stage` property.

		Parameters
		----------
		attached : bool
			`True` if part is attached to the sample stage,
			`False` if not.
		"""
		self.attached_to_stage = attached

	def _set_static_if_no_drifts(self):
		"""Set the object to 'static' if it does not drift,
		i.e., it is not moving.
		In this case, its coordinate system does not need to
		be re-assembled for each frame.
		"""

		self._static = False

		if self.attached_to_stage is False:
			# Count drifts:
			if self.center.has_drifts() or self.u.has_drifts() or self.w.has_drifts():
				return

			for dev in self.deviations:
				if dev.has_drifts():
					return

			for ldev in self.legacy_deviations:
				if ldev.has_drifts():
					return

			self._static = True

	def set_geometry(self, json_geometry_object:dict, stage_coordinate_system:'CoordinateSystem'=None) -> bool:
		"""
		Set up the part from a CTSimU JSON geometry definition.
		The `stage_coordinate_system` must only be provided if this
		part is attached to the stage.

		Parameters
		----------
		json_geometry_object : dict
			A CTSimU geometry object, as imported from a JSON structure.

		Returns
		-------
		success : bool
			`True` on success, `False` if an error occurred.
		"""
		if stage_coordinate_system is None:
			stage_coordinate_system = ctsimu_world

		self.reset()

		geo = json_geometry_object

		# Try to set up the part from world coordinate notation (x, y, z).
		# We also have to support legacy spelling of "centre" ;-)
		if (json_exists_and_not_null(geo, ["center", "x"]) or \
			json_exists_and_not_null(geo, ["centre", "x"])) and \
		   (json_exists_and_not_null(geo, ["center", "y"]) or \
			json_exists_and_not_null(geo, ["centre", "y"])) and \
		   (json_exists_and_not_null(geo, ["center", "z"]) or \
			json_exists_and_not_null(geo, ["centre", "z"])):
			# *******************************
			#           Part is in
			#     WORLD COORDINATE SYSTEM
			# *******************************
			self.attach_to_stage(attached=False)

			# Center
			# ------
			if self.center.set_from_json(json_extract_from_possible_keys(geo, [["center"], ["centre"]])):
				# success
				pass
			else:
				raise Exception(f"Part '{self.name}': failed setting the object center from the JSON dictionary.")
				return False

			# Orientation
			# -----------
			# Vectors can be either u, w (for source, stage, detector)
			# or r, t (for samples).
			if self.u.set_from_json(json_extract_from_possible_keys(geo, [["vector_u"], ["vector_r"]])) and \
			   self.w.set_from_json(json_extract_from_possible_keys(geo, [["vector_w"], ["vector_t"]])):
				# success
				pass
			else:
				raise Exception(f"Part {self.name} is placed in world coordinate system, but its vectors u and w (or r and t, for samples) are not properly defined (each with an x, y and z component).")
				return False

		elif (json_exists_and_not_null(geo, ["center", "u"]) or \
		      json_exists_and_not_null(geo, ["centre", "u"])) and \
		     (json_exists_and_not_null(geo, ["center", "v"]) or \
		      json_exists_and_not_null(geo, ["centre", "v"])) and \
		     (json_exists_and_not_null(geo, ["center", "w"]) or \
		      json_exists_and_not_null(geo, ["centre", "w"])):
			# *******************************
			#           Part is in
			#     STAGE COORDINATE SYSTEM
			# *******************************
			self.attach_to_stage(attached=True)

			# Center
			# ------
			if self.center.set_from_json(json_extract_from_possible_keys(geo, [["center"], ["centre"]])):
				# success
				pass
			else:
				raise Exception(f"Part '{self.name}': failed setting the object center from the JSON dictionary.")
				return False

			# Orientation
			# -----------
			# Vectors can only be r, t
			# (because only samples can be attached to the stage).
			if self.u.set_from_json(json_extract(geo, ["vector_r"])) and \
			   self.w.set_from_json(json_extract(geo, ["vector_t"])):
				# success
				pass
			else:
				raise Exception(f"Part {self.name} is placed in stage system, but its vectors r and t are not properly defined (each with a u, v and w component).")
				return False

		else:
			raise Exception(f"Failed to set geometry for part '{self.name}'. Found no valid center definition in JSON file.")
			return False

		# *******************************
		#     DEVIATIONS
		# *******************************
		if json_exists_and_not_null(geo, ["deviations"]):
			devs = json_extract(geo, ["deviations"])
			if isinstance(devs, list):
				# Go through all elements in the deviations array
				# and add them to this part's list of deviations.
				for dev in devs:
					new_deviation = Deviation()
					if new_deviation.set_from_json(dev):
						self.deviations.append(new_deviation)
					else:
						raise Exception(f"An error occurred when setting a geometrical deviation (translation) for part '{self.name}'.")
						return False
			elif isinstance(devs, dict):
				# Only one drift defined directly as a JSON object?
				# Actually not supported by file format,
				# but let's be generous and try...
				new_deviation = Deviation()
				if new_deviation.set_from_json(devs):
					self.deviations.append(new_deviation)
				else:
					raise Exception(f"An error occurred when setting a geometrical deviation (rotation) for part '{self.name}'.")
					return False
			else:
				raise Exception(f"Error reading geometrical deviations for part '{self.name}'")

		# Support for legacy deviations, prior to
		# file format version 0.9:
		# ------------------------------------------
		if json_exists_and_not_null(geo, ["deviation"]):
			known_to_recon = True
			if json_exists_and_not_null(geo, ["deviation", "known_to_reconstruction"]):
				known_to_recon = get_value_in_native_unit("bool", geo, ["deviation", "known_to_reconstruction"])

			for axis in ctsimu_valid_axes:
				# Deviations in position
				# -------------------------------------
				# Positional deviations along sample axes r, s, t
				# have not been part of the legacy file formats
				# prior to version 0.9, but we still add them here
				# because now we easily can... ;-)
				if json_exists_and_not_null(geo, ["deviation", "position", axis]):
					pos_dev = Deviation()
					pos_dev.set_type("translation")
					pos_dev.set_axis(axis)
					pos_dev.set_known_to_reconstruction(known_to_recon)
					if pos_dev.amount.set_from_json(json_extract(geo, ["deviation", "position", axis])):
						# Legacy_deviations not necessary here
						# because positional translations are fully
						# compatible with the new file format:
						self.deviations.append(pos_dev)
					else:
						raise Exception(f"An error occurred when setting a geometrical deviation (translation) for part '{self.name}'.")
						return False

			for axis in ctsimu_valid_axes:
				# Deviations in rotation
				# -------------------------------------
				# File formats prior to version 0.9 only supported
				# rotations around u, v and w, in the order wv'u'',
				# and ts'r'' for samples. We need to take care
				# to keep this order here; it is ensured by the order
				# of elements in the `valid_axes` list. This means we
				# also add support for x, y, z (zy'x''),
				# just because we can.
				if json_exists_and_not_null(geo, ["deviation", "rotation", axis]):
					rot_dev = Deviation()
					rot_dev.set_type("rotation")

					# Prior to 0.9, all deviations were meant to take place
					# before the stage rotation. This means they need to be
					# stored as scene vectors to designate a constant deviation axis.
					rot_dev.set_axis(axis)
					rot_dev.set_known_to_reconstruction(known_to_recon)
					if rot_dev.amount.set_from_json(json_extract(geo, ["deviation", "rotation", axis])):
						# success
						self.legacy_deviations.append(rot_dev)
					else:
						raise Exception(f"An error occurred when setting a geometrical deviation (rotation) for part '{self.name}'.")
						return False

		self.set_frame(stage_coordinate_system, frame=0, nFrames=1, w_rotation_in_rad=0)
		return True
