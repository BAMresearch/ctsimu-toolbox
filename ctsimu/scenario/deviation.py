# -*- coding: UTF-8 -*-
"""
Geometrical deviation of a coordinate system, i.e. a translation or a rotation
with respect to an arbitrary axis.
"""

from ..helpers import *
from ..primitives import Vector, Matrix
from ..geometry import CoordinateSystem, ctsimu_world
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

	def __init__(self, pivot_reference:str="sample", _root=None):
		self._root = _root  # root scenario object

		# The axis and pivot point are `Scenevector`
		# objects that can handle vector drifts and
		# conversion between coordinate systems:
		self.axis  = Scenevector(native_unit=None, _root=self._root)
		self.pivot = Scenevector(native_unit="mm", _root=self._root)

		# Set a standard pivot which refers to the object's center:
		self.pivot.set_simple(0, 0, 0)
		# Set "sample" as the "most local" pivot reference.
		# Will be converted to "local" if necessary:
		self.pivot.set_reference(pivot_reference)

		# The transformation amount is a
		# `Parameter` that can handle drifts.
		self.amount = Parameter(native_unit=None, standard_value=0, _root=self._root)

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
			if self.type == "rotation":
				self.amount.set_native_unit("rad")
			elif self.type == "translation":
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

	def json_dict(self) -> dict:
		"""Create a dictionary of this deviation for a CTSimU JSON file."""
		jd = dict()

		jd["type"] = self.type
		jd["axis"] = self.axis.json_dict()
		if self.type == "rotation":
			jd["pivot"] = self.pivot.json_dict()

		jd["amount"] = self.amount.json_dict()
		jd["known_to_reconstruction"] = self.known_to_reconstruction

		return jd

	def deviate(self, coordinate_system:'CoordinateSystem', frame:float, n_frames:int, reconstruction=False, attached_to_stage:bool=False, stage_coordinate_system:'CoordinateSystem'=None) -> 'CoordinateSystem':
		"""
		Apply this deviation to the given coordinate system.

		This function checks if the deviation's `known_to_reconstruction`
		property applies to this case,

		Parameters
		----------
		coordinate_system : ctsimu.geometry.CoordinateSystem
			The coordinate system that will be deviated.

		frame : float
			Frame number.

		n_frames : int
			Total number of frames in scan.

		reconstruction : bool
			If `True`, the deviation is only applied if it is
			known to the reconstruction software. The same applies
			to whether its drifts are considered or not.
			If `False`, the full deviation and all of its drifts will
			be applied to the coordinate system.

		attached_to_stage : bool
			`True` if the coordinate system is attached to the sample stage,
			`False` if its reference coordinate system is a world coordinate system.

		stage_coordinate_system : ctsimu.geometry.CoordinateSystem, optional
			The stage coordinate system. Only necessary if the part is attached
			to the sample stage. Otherwise, `None` can be passed.

		Returns
		-------
		deviated_cs : ctsimu.geometry.CoordinateSystem
			Coordinate system with the applied deviation. Note that the
			original `coordinate_system` passed as an argument will be
			altered equally.
		"""

		if (self.known_to_reconstruction is True) or (reconstruction is False):
			amount = self.amount.set_frame_and_get_value(
				frame=frame,
				n_frames=n_frames,
				reconstruction=reconstruction
			)

			if self.type == "translation":
				if self.amount.native_unit == "mm":
					if attached_to_stage is False:
						# Object in world coordinate system.
						# --------------------------------------
						# The deviation axis can undergo drifts and could
						# be expressed in any coordinate system (world, local, sample).
						# Therefore, the axis is a `Scenevector`, which can
						# give us the translation vector for the current frame:
						translation_axis = self.axis.direction_in_world(
							local=coordinate_system,
							sample=ctsimu_world,
							frame=frame,
							n_frames=n_frames,
							reconstruction=reconstruction
						)

						coordinate_system.translate_in_direction(direction=translation_axis, distance=amount)
					else:
						# Object is in stage coordinate system.
						# --------------------------------------
						translation_axis = self.axis.direction_in_local(
							local=stage_coordinate_system,
							sample=coordinate_system,
							frame=frame,
							n_frames=n_frames,
							reconstruction=reconstruction
						)

						coordinate_system.translate_in_direction(direction=translation_axis, distance=amount)
				else:
					raise Exception(f"CTSimU Deviation: deviate: Translational deviation must be given in a unit of length. The given unit '{self.amount.native_unit}' cannot be used for a translation.")
			elif self.type == "rotation":
				if self.amount.native_unit == "rad":
					if attached_to_stage is False:
						# Object in world coordinate system.
						# --------------------------------------
						rotation_axis = self.axis.direction_in_world(
							local=coordinate_system,
							sample=ctsimu_world,
							frame=frame,
							n_frames=n_frames,
							reconstruction=reconstruction
						)
						pivot_point = self.pivot.point_in_world(
							local=coordinate_system,
							sample=ctsimu_world,
							frame=frame,
							n_frames=n_frames,
							reconstruction=reconstruction
						)

						coordinate_system.rotate_around_pivot_point(
							axis=rotation_axis, angle=amount, pivot=pivot_point
						)
					else:
						# Object is in stage coordinate system.
						# --------------------------------------
						rotation_axis = self.axis.direction_in_local(
							local=stage_coordinate_system,
							sample=coordinate_system,
							frame=frame,
							n_frames=n_frames,
							reconstruction=reconstruction
						)
						pivot_point = self.pivot.point_in_local(
							local=stage_coordinate_system,
							sample=coordinate_system,
							frame=frame,
							n_frames=n_frames,
							reconstruction=reconstruction
						)

						coordinate_system.rotate_around_pivot_point(
							axis=rotation_axis, angle=amount, pivot=pivot_point
						)
				else:
					raise Exception(f"CTSimU Deviation: deviate: Rotational deviation must be given in an angular unit. The given unit '{self.amount.native_unit}' cannot be used for a rotation.")

		return coordinate_system