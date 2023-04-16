# -*- coding: UTF-8 -*-
"""
Geometrical deviation of a coordinate system, i.e. a translation or a rotation
with respect to an arbitrary axis.
"""

from ..helpers import *
from ..primitives import Vector, Matrix
from .parameter import Parameter
from .scenevector import Scenevector

class Deviation:
	""" Geometrical deviation of a coordinate system,
	i.e. a translation or a rotation with respect to one of the
	axes x, y, z (world), u, v, w (local), or r, s, t (sample)
	or any other arbitrary vector.

	Like any parameter, they can have drifts, which means they
	can change over time.

	Attributes
	----------
	type : str
		Transformation type. Can be `"translation"` or `"rotation"`.

	amount : ctsimu.scenario.parameter.Parameter
		Transformation amount: length for a translation, angle for a rotation.

	axis : ctsimu.scenario.scenevector.Scenevector
		Transformation axis: direction for a translation, rotation axis for
		a rotation.

	pivot : ctsimu.scenario.scenevector.Scenevector
		For rotational deviations: position of the rotation's pivot point.

	known_to_reconstruction : bool
		Should this geometrical deviation be considered by the reconstruction software?
		This attribute is obeyed when calculating projection matrices.
	"""

	def __init__(self):
		# The axis and pivot point are `Scenevector`
		# objects that can handle vector drifts and
		# conversion between coordinate systems:
		self.axis  = Scenevector(native_unit=None)
		self.pivot = Scenevector(native_unit="mm")

		# Set a standard pivot which refers to the object's center:
		self.pivot.set_simple(0, 0, 0)
		# Set "sample" as the "most local" pivot reference.
		# Will be converted to "local" if necessary:
		self.pivot.set_reference("sample")

		# The transformation amount is a
		# `Parameter` that can handle drifts.
		self.amount = Parameter(native_unit=None, standard_value=0)

		self.set_type("translation")
		self.known_to_reconstruction = False

	def has_drifts(self) -> bool:
		"""Does the deviation specify any drift?

		Returns
		-------
		has_drifts : bool
			`True` if the deviation specifies any drift, otherwise `False`.
		"""
		return self.amount.has_drifts() or self.axis.has_drifts() or self.pivot.has_drifts()

	def set_type(self, transformation_type):
		"""Set the transformation type.

		Parameters
		----------
		transformation_type : str
			Transformation type. `"translation"` or `"rotation"`.
		"""
		if (transformation_type == "rotation") or (transformation_type == "translation"):
			self.type = transformation_type

			# Set the correct native unit for the amount:
			if self.transformation_type == "rotation":
				self.amount.set_native_unit("rad")
			elif self.transformation_type == "translation":
				self.amount.set_native_unit("mm")
		else:
			raise Exception(f"CTSimU Deviation: set_type: Not a valid transformation type: '{transformation_type}'. Valid types are 'translation' and 'rotation'.")

	def set_axis(self, axis):
		"""Set the transformation axis.
		Can be a string identifier: `"x"`, `"y"`, `"z"` (for world axes),
		`"u"`, `"v"`, `"w"` (for local axes), `"r"`, `"s"`, `"t"` (for sample axes),
		or an arbitrary `ctsimu.scenario.scenevector.Scenevector` object.

		Parameters
		----------
		axis : str or ctsimu.scenario.scenevector.Scenevector
			Transformation axis.
		"""
		if isinstance(axis, str):
			if axis in ctsimu_world_axis_designations:
				# Given axis is "x", "y" or "z"
				# -> vector in world coordinate system
				self.axis.set_reference("world")
				if axis == "x": self.axis.set_simple(1, 0, 0)
				if axis == "y": self.axis.set_simple(0, 1, 0)
				if axis == "z": self.axis.set_simple(0, 0, 1)
			elif axis in ctsimu_local_axis_designations:
				# Given axis is "u", "v" or "w"
				# -> vector in local coordinate system
				self.axis.set_reference("local")
				if axis == "u": self.axis.set_simple(1, 0, 0)
				if axis == "v": self.axis.set_simple(0, 1, 0)
				if axis == "w": self.axis.set_simple(0, 0, 1)
			elif axis in ctsimu_sample_axis_designations:
				# Given axis is "r", "s" or "t"
				# -> vector in sample coordinate system
				self.axis.set_reference("sample")
				if axis == "r": self.axis.set_simple(1, 0, 0)
				if axis == "s": self.axis.set_simple(0, 1, 0)
				if axis == "t": self.axis.set_simple(0, 0, 1)
			else:
				raise Exception(f"CTSimU Deviation: set_axis: Not a valid axis designation string: '{axis}'. Valid options are: {ctsimu_axis_strings}.")
		elif isinstance(axis, Scenevector):
			self.axis = axis
		else:
			raise Exception(f"CTSimU Deviation: set_axis: Not a valid deviation axis type: '{type(axis)}'. Valid options are: str, Scenevector.")

	def set_pivot(self, pivot:'Scenevector'):
		"""Set the pivot point for rotations.

		Parameters
		----------
		pivot : Scenevector
		"""
		if isinstance(pivot, Scenevector):
			self.pivot = pivot
		else:
			raise Exception(f"CTSimU Deviation: set_pivot: expects an instance of 'Scenevector', but given pivot is of type {type(pivot)}.")

	def set_known_to_reconstruction(self, known_to_reconstruction:bool):
		"""Sets the `known_to_reconstruction` attribute.

		Parameters
		----------
		known_to_reconstruction : bool
			Should this geometrical deviation be considered by the reconstruction software?
			This attribute is obeyed when calculating projection matrices.
		"""
		if known_to_reconstruction == True:
			self.known_to_reconstruction = True
		else:
			self.known_to_reconstruction = False

	def set_amount_from_json(self, json_parameter_object:dict) -> bool:
		"""Set the deviation's amount from a JSON object, which
		is a parameter with a value and potentially a drift.

		This function is usually not called from the outside,
		but used by `set_from_json`.

		Parameters
		----------
		json_parameter_object : dict
			A CTSimU parameter structure, as imported from a JSON scenario description.

		Returns
		-------
		success : bool
			`True` on success, `False` if an error occurred.
		"""
		return self.amount.set_from_json(json_parameter_object)

	def set_from_json(self, json_deviation_object:dict) -> bool:
		"""Set the deviation from a CTSimU JSON deviation structure.

		Parameters
		----------
		json_deviation_object : dict
			A CTSimU deviation structure, as imported from a JSON scenario description.

		Returns
		-------
		success : bool
			`True` on success, `False` if an error occurred.

		Raises
		------
		Exception
			When an error occurred.
		"""
		if json_exists_and_not_null(json_deviation_object, ["type"]):
			self.set_type(get_value(json_deviation_object, ["type"]))
		else:
			raise Exception('CTSimU Deviation: set_from_json: a geometrical deviation must provide a "type": either "rotation" or "translation".')
			return False

		# Transformation axis:
		if json_exists_and_not_null(json_deviation_object, ["axis"]):
			axis_obj = json_extract(json_deviation_object, ["axis"])
			if isinstance(axis_obj, str):
				# Axis is given as a single string.
				if axis_obj in ctsimu_axis_strings:
					self.set_axis(axis_obj)
				else:
					raise Exception(f'CTSimU Deviation: set_from_json: The deviation "axis" string is incorrect: must be any of {ctsimu_axis_strings} or a free vector definition.')
					return False
			elif isinstance(axis_obj, dict):
				# Free vector definition.
				if self.axis.set_from_json(axis_obj):
					# Success.
					pass
				else:
					raise Exception(f'CTSimU Deviation: set_from_json: Failed to set up deviation axis from JSON file. Vector definition seems to be incorrect.')
					return False
			else:
				raise Exception("CTSimU Deviation: set_from_json: Failed to set up deviation axis from JSON structure.")

		# Pivot point for rotations:
		# Set a standard pivot which refers to the object's center:
		self.pivot.set_simple(0, 0, 0)

		# Set "sample" as the "most local" pivot reference.
		# Will be converted to "local" if necessary:
		self.pivot.set_reference("sample")

		if json_exists_and_not_null(json_deviation_object, ["pivot"]):
			# If another pivot is defined in the
			# JSON file, take this one instead...
			if self.pivot.set_from_json(json_extract(json_deviation_object, ["pivot"])):
				# Success
				pass
			else:
				raise Exception("CTSimU Deviation: set_from_json: Failed to set up deviation's pivot point from JSON structure. Vector definition seems to be incorrect.")
				return False

		self.set_amount_from_json(json_extract(json_deviation_object, ["amount"]))
		self.set_known_to_reconstruction(get_value_in_native_unit(
			native_unit="bool",
			dictionary=json_deviation_object,
			keys=["known_to_reconstruction"],
			fail_value=True
		))

		return True