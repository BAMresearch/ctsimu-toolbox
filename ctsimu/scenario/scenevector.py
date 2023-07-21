# -*- coding: UTF-8 -*-
"""
A scene vector is a 3D vector that knows its reference coordinate system.
It can be converted between these coordinate systems and handles drifts.
"""

from ..helpers import *
from ..primitives import Vector, Matrix
from ..geometry import *
from .parameter import Parameter

class Scenevector:
	"""A scene vector is a 3D vector that knows the type of its
	reference coordinate system, given as world, local or sample.
	It provides functions to convert between these coordinate systems
	and it can handle drifts. Therefore, all three vector components
	are stored as `ctsimu.scenario.parameter.Parameter` objects.

	Useful for vectors that can change due to drifts,
	such as rotation axis and pivot point of a deviation,
	or, in general, the coordinate system vectors.

	Attributes
	----------
	reference : str
		Identifies the reference coordinate system.
		Can be: `"world"`, `local"`, or `"sample"`.

	c0 : ctsimu.scenario.parameter.Parameter
		Value of the first vector component.

	c1 : ctsimu.scenario.parameter.Parameter
		Value of the second vector component.

	c2 : ctsimu.scenario.parameter.Parameter
		Value of the third vector component.
	"""

	def __init__(self, native_unit:str=None):
		"""For initialization, the native unit can be given.

		Parameters
		----------
		native_unit : str, optional
			The native unit of the vector's components. Standard is `None`,
			but `"mm"` would make sense for positions in space.
			Possible values: `None`, `"mm"`, `"rad"`, `"deg"`, `"s"`, `"mA"`, `"kV"`, `"g/cm^3"`, `"lp/mm"`, `"bool"`, `"string"`.
		"""
		self.c0 = Parameter(native_unit=native_unit) # 1st vector component
		self.c1 = Parameter(native_unit=native_unit) # 2nd vector component
		self.c2 = Parameter(native_unit=native_unit) # 3rd vector component

		# Reference coordinate system:
		self.reference = "world" # "world", local", "sample"

	def __str__(self):
		s = "({c0}, {c1}, {c2}) in {ref}".format(
			c0 = self.c0.current_value,
			c1 = self.c1.current_value,
			c2 = self.c2.current_value,
			ref = self.reference
		)

		return s

	def has_drifts(self) -> bool:
		"""Does the vector drift during the CT scan?

		Returns
		-------
		has_drifts : bool
			`True` if any vector component defines a drift,
			`False` if none of the vector components defines any drift.
		"""
		return (self.c0.has_drifts()) or (self.c1.has_drifts()) or (self.c2.has_drifts())

	def set_reference(self, reference:str):
		"""Set the reference coordinate system.
		Can be `"world"`, `"local"` or `"sample"`.

		Parameters
		----------
		reference : str
			String that identifies the reference coordinate system
			(see above for valid strings).
		"""
		valid_references = ["world", "local", "sample"]
		if reference in valid_references:
			self.reference = reference
		else:
			raise Exception(f"Scenevector: set_reference(): not a valid reference coordinate system: '{reference}'. Valid reference coordinate systems are: {valid_references}.")

	def set_native_unit(self, native_unit:str):
		"""Set native unit of vector components.
		Necessary for the location of points such as
		the center points of coordinate systems,
		usually given in `"mm"` as native unit.

		Parameters
		----------
		native_unit : str
			The native unit of the vector's components. Standard is `None`,
			but `"mm"` would make sense for positions in space.
			Possible values: `None`, `"mm"`, `"rad"`, `"deg"`, `"s"`, `"mA"`, `"kV"`, `"g/cm^3"`, `"lp/mm"`, `"bool"`, `"string"`.
		"""
		if is_valid_native_unit(native_unit):
			self.c0.set_native_unit(native_unit)
			self.c1.set_native_unit(native_unit)
			self.c2.set_native_unit(native_unit)

	def set_simple(self, c0:float, c1:float, c2:float):
		"""Set a simple scene vector from three numbers,
		results in a scene vector without drifts.

		Parameters
		----------
		c0 : float
			Value of the first vector component.

		c1 : float
			Value of the second vector component.

		c2 : float
			Value of the third vector component.
		"""
		self.c0.set_standard_value(c0)
		self.c1.set_standard_value(c1)
		self.c2.set_standard_value(c2)

		# Delete all drifts and set current value to standard value:
		self.c0.reset()
		self.c1.reset()
		self.c2.reset()

	def set_component(self, i:int, p:'Parameter'):
		"""Set the `i`th vector component to the given parameter `p`
		(which must be a `ctsimu.scenario.parameter.Parameter` object).

		Parameters
		----------
		i : int
			Index of the vector component (`0`, `1` or `2`).

		p : ctsimu.scenario.parameter.Parameter
			Parameter object for the given vector component.
		"""
		if isinstance(p, Parameter):
			if i == 0:
				self.c0 = p
			elif i == 1:
				self.c1 = p
			elif i == 2:
				self.c2 = p
		else:
			raise Exception(f"CTSimU Scenevector: set_component(): the given parameter is not an instance of class 'ctsimu.scenario.parameter.Parameter', but of type {type(p)}.")

	def standard_vector(self) -> 'Vector':
		"""A `ctsimu.primitives.Vector` that represents this
		scene vector in its standard orientation
		(i.e., without any drifts applied).

		Returns
		-------
		v : ctsimu.primitives.Vector
			Standard vector (without any drifts applied).
		"""
		return Vector(
			x=self.c0.get_standard_value(),
			y=self.c1.get_standard_value(),
			z=self.c2.get_standard_value()
		)

	def drift_vector(self, frame:float, nFrames:int, only_known_to_reconstruction:bool=False) -> 'Vector':
		"""A `ctsimu.primitives.Vector` that represents
		only the drift values for the given `frame` number.

		Can later be added to the standard value to get
		the resulting vector respecting all drifts.

		Parameters
		----------
		frame : float
			Current frame number.

		nFrames : int
			Total number of frames in CT scan.

		only_known_to_reconstruction : bool
			Obey only those drifts that are labeled as known to the
			reconstruction software?

		Returns
		-------
		drift_vector : ctsimu.primitives.Vector
			Vector that contains the drift values for the given `frame` number.
		"""

		return Vector(
			x=self.c0.get_total_drift_value_for_frame(frame, nFrames, only_known_to_reconstruction),
			y=self.c1.get_total_drift_value_for_frame(frame, nFrames, only_known_to_reconstruction),
			z=self.c2.get_total_drift_value_for_frame(frame, nFrames, only_known_to_reconstruction)
		)

	def vector_for_frame(self, frame:float, nFrames:int, only_known_to_reconstruction:bool=False) -> 'Vector':
		"""A `ctsimu.primitives.Vector` for the given `frame` number,
		respecting all drifts.

		Parameters
		----------
		frame : float
			Current frame number.

		nFrames : int
			Total number of frames in CT scan.

		only_known_to_reconstruction : bool
			Obey only those drifts that are labeled as known to the
			reconstruction software?

		Returns
		-------
		vec : ctsimu.primitives.Vector
			Vector for the given `frame` number.
		"""
		return Vector(
			x=self.c0.set_frame_and_get_value(frame, nFrames, only_known_to_reconstruction),
			y=self.c1.set_frame_and_get_value(frame, nFrames, only_known_to_reconstruction),
			z=self.c2.set_frame_and_get_value(frame, nFrames, only_known_to_reconstruction)
		)

	def point_in_world(self, local:'CoordinateSystem', sample:'CoordinateSystem', frame:float, nFrames:int, only_known_to_reconstruction:bool=False) -> 'Vector':
		"""A `ctsimu.primitives.Vector` for point coordinates
		in terms of the world coordinate system for the given `frame` number,
		respecting all relevant drifts.

		Parameters
		----------
		local : ctsimu.geometry.CoordinateSystem
			A coordinate system that represents the object's local coordinate
			system in terms of world coordinates.

		sample : ctsimu.geometry.CoordinateSystem
			A coordinate system that represents the sample in terms of the
			stage coordinate system. If you don't want to convert from a
			sample vector, pass `None`.

		frame : float
			The number of the current frame.

		nFrames : int
			The total number of frames in the CT scan.

		only_known_to_reconstruction : bool
			Only handle drifts that are labeled as known to the reconstruction software?

		Returns
		-------
		vec_world : ctsimu.primitives.Vector
			A vector in terms of the world coordinate system.
		"""

		v = self.vector_for_frame(frame, nFrames, only_known_to_reconstruction)

		if self.reference == "world":
			# Already in world.
			return v
		elif self.reference == "local":
			# Convert from local to world.
			v_in_world = change_reference_frame_of_point(v, local, ctsimu_world)
			return v_in_world
		elif self.reference == "sample":
			# The sample's "world" is the stage (here: local).
			# To get the sample coordinates in stage coordinates,
			# we therefore transform to the world a first time...
			# ...and a second time to transform it
			# from the stage to the world:
			v_in_stage = change_reference_frame_of_point(v, sample, ctsimu_world)
			v_in_world = change_reference_frame_of_point(v_in_stage, local, ctsimu_world)
			return v_in_world

	def point_in_local(self, local:'CoordinateSystem', sample:'CoordinateSystem', frame:float, nFrames:int, only_known_to_reconstruction:bool=False) -> 'Vector':
		"""A `ctsimu.primitives.Vector` for point coordinates
		in terms of the local coordinate system for the given `frame` number,
		respecting all relevant drifts.

		Parameters
		----------
		local : ctsimu.geometry.CoordinateSystem
			A coordinate system that represents the object's local coordinate
			system in terms of world coordinates.

		sample : ctsimu.geometry.CoordinateSystem
			A coordinate system that represents the sample in terms of the
			stage coordinate system. If you don't want to convert from a
			sample vector, pass `None`.

		frame : float
			The number of the current frame.

		nFrames : int
			The total number of frames in the CT scan.

		only_known_to_reconstruction : bool
			Only handle drifts that are labeled as known to the reconstruction software?

		Returns
		-------
		vec_local : ctsimu.primitives.Vector
			A vector in terms of the local coordinate system.
		"""

		v = self.vector_for_frame(frame, nFrames, only_known_to_reconstruction)

		if self.reference == "world":
			# Convert from world to local.
			v_in_local = change_reference_frame_of_point(v, ctsimu_world, local)
			return v_in_local
		elif self.reference == "local":
			# Already in local.
			return v
		elif self.reference == "sample":
			# The sample's "world" is the stage (here: local).
			# To get the sample coordinates in stage coordinates,
			# we therefore transform to the world.
			v_in_stage = change_reference_frame_of_point(v, sample, ctsimu_world)
			return v_in_local

	def point_in_sample(self, stage:'CoordinateSystem', sample:'CoordinateSystem', frame:float, nFrames:int, only_known_to_reconstruction:bool=False) -> 'Vector':
		"""A `ctsimu.primitives.Vector` for point coordinates
		in terms of the sample coordinate system for the given `frame` number,
		respecting all relevant drifts.

		The sample must be attached to the stage, otherwise use `in_local`
		instead. Note that a direct conversion from one sample CS to another is
		not possible with this function.

		Parameters
		----------
		stage : ctsimu.geometry.CoordinateSystem
			A coordinate system that represents the stage's local coordinate
			system in terms of world coordinates.

		sample : ctsimu.geometry.CoordinateSystem
			A coordinate system that represents the sample in terms of the
			stage coordinate system. If you don't want to convert from a
			sample vector, pass `None`.
			This must be the same sample coordinate system to which this
			scene vector refers.

		frame : float
			The number of the current frame.

		nFrames : int
			The total number of frames in the CT scan.

		only_known_to_reconstruction : bool
			Only handle drifts that are labeled as known to the reconstruction software?

		Returns
		-------
		vec_sample : ctsimu.primitives.Vector
			A vector in terms of the sample coordinate system.
		"""

		v = self.vector_for_frame(frame, nFrames, only_known_to_reconstruction)

		if self.reference == "world":
			# From world to stage...
			# ...and a second time from world to sample
			# because the stage is the sample's "world":
			v_in_stage  = change_reference_frame_of_point(v, ctsimu_world, stage)
			v_in_sample = change_reference_frame_of_point(v_in_stage, ctsimu_world, sample)
			return v_in_sample
		elif self.reference == "local":
			# Convert from stage to sample. Because
			# the stage is the sample's world, we actually
			# convert from world to sample.
			v_in_sample = change_reference_frame_of_point(v, ctsimu_world, sample)
			return v_in_sample
		elif self.reference == "sample":
			# Already in sample coordinate system:
			return v

	def direction_in_world(self, local:'CoordinateSystem', sample:'CoordinateSystem', frame:float, nFrames:int, only_known_to_reconstruction:bool=False) -> 'Vector':
		"""A `ctsimu.primitives.Vector` for a direction in terms of the world
		coordinate system for the given `frame` number, respecting all relevant drifts.

		Parameters
		----------
		local : ctsimu.geometry.CoordinateSystem
			A coordinate system that represents the object's local coordinate
			system in terms of world coordinates.

		sample : ctsimu.geometry.CoordinateSystem
			A coordinate system that represents the sample in terms of the
			stage coordinate system. If you don't want to convert from a
			sample vector, pass `None`.

		frame : float
			The number of the current frame.

		nFrames : int
			The total number of frames in the CT scan.

		only_known_to_reconstruction : bool
			Only handle drifts that are labeled as known to the reconstruction software?

		Returns
		-------
		vec_world : ctsimu.primitives.Vector
			A vector in terms of the world coordinate system.
		"""

		v = self.vector_for_frame(frame, nFrames, only_known_to_reconstruction)

		if self.reference == "world":
			# Already in world.
			return v
		elif self.reference == "local":
			# Convert from local to world.
			v_in_world = change_reference_frame_of_direction(v, local, ctsimu_world)
			return v_in_world
		elif self.reference == "sample":
			# The sample's "world" is the stage (here: local).
			# To get the sample coordinates in stage coordinates,
			# we therefore transform to the world a first time...
			# ...and a second time to transform it
			# from the stage to the world:
			v_in_stage = change_reference_frame_of_direction(v, sample, ctsimu_world)
			v_in_world = change_reference_frame_of_direction(v_in_stage, local, ctsimu_world)
			return v_in_world

	def direction_in_local(self, local:'CoordinateSystem', sample:'CoordinateSystem', frame:float, nFrames:int, only_known_to_reconstruction:bool=False) -> 'Vector':
		"""A `ctsimu.primitives.Vector` for a direction in terms of the local
		coordinate system for the given `frame` number, respecting all relevant drifts.

		Parameters
		----------
		local : ctsimu.geometry.CoordinateSystem
			A coordinate system that represents the object's local coordinate
			system in terms of world coordinates.

		sample : ctsimu.geometry.CoordinateSystem
			A coordinate system that represents the sample in terms of the
			stage coordinate system. If you don't want to convert from a
			sample vector, pass `None`.

		frame : float
			The number of the current frame.

		nFrames : int
			The total number of frames in the CT scan.

		only_known_to_reconstruction : bool
			Only handle drifts that are labeled as known to the reconstruction software?

		Returns
		-------
		vec_local : ctsimu.primitives.Vector
			A vector in terms of the local coordinate system.
		"""

		v = self.vector_for_frame(frame, nFrames, only_known_to_reconstruction)

		if self.reference == "world":
			# Convert from world to local.
			v_in_local = change_reference_frame_of_direction(v, ctsimu_world, local)
			return v_in_local
		elif self.reference == "local":
			# Already in local.
			return v
		elif self.reference == "sample":
			# The sample's "world" is the stage (here: local).
			# To get the sample coordinates in stage coordinates,
			# we therefore transform to the world.
			v_in_stage = change_reference_frame_of_direction(v, sample, ctsimu_world)
			return v_in_local

	def direction_in_sample(self, stage:'CoordinateSystem', sample:'CoordinateSystem', frame:float, nFrames:int, only_known_to_reconstruction:bool=False) -> 'Vector':
		"""A `ctsimu.primitives.Vector` for a direction in terms of the sample
		coordinate system for the given `frame` number, respecting all relevant drifts.

		The sample must be attached to the stage, otherwise use `in_local`
		instead. Note that a direct conversion from one sample CS to another is
		not possible with this function.

		Parameters
		----------
		stage : ctsimu.geometry.CoordinateSystem
			A coordinate system that represents the stage's local coordinate
			system in terms of world coordinates.

		sample : ctsimu.geometry.CoordinateSystem
			A coordinate system that represents the sample in terms of the
			stage coordinate system. If you don't want to convert from a
			sample vector, pass `None`.
			This must be the same sample coordinate system to which this
			scene vector refers.

		frame : float
			The number of the current frame.

		nFrames : int
			The total number of frames in the CT scan.

		only_known_to_reconstruction : bool
			Only handle drifts that are labeled as known to the reconstruction software?

		Returns
		-------
		vec_sample : ctsimu.primitives.Vector
			A vector in terms of the sample coordinate system.
		"""

		v = self.vector_for_frame(frame, nFrames, only_known_to_reconstruction)

		if self.reference == "world":
			# From world to stage...
			# ...and a second time from world to sample
			# because the stage is the sample's "world":
			v_in_stage  = change_reference_frame_of_direction(v, ctsimu_world, stage)
			v_in_sample = change_reference_frame_of_direction(v_in_stage, ctsimu_world, sample)
			return v_in_sample
		elif self.reference == "local":
			# Convert from stage to sample. Because
			# the stage is the sample's world, we actually
			# convert from world to sample.
			v_in_sample = change_reference_frame_of_direction(v, ctsimu_world, sample)
			return v_in_sample
		elif self.reference == "sample":
			# Already in sample coordinate system:
			return v

	def set_from_json(self, json_vector_object:dict) -> bool:
		"""Set up the scene vector from a CTSimU JSON object
		that describes a three-component vector.

		Parameters
		----------
		json_vector_object : dict
			CTSimU vector object, as imported from a JSON scenario description.

		Returns
		-------
		success : bool
			`True` on success, `False` when an error occurred.
		"""

		if json_exists_and_not_null(json_vector_object, keys=["x"]) and \
		   json_exists_and_not_null(json_vector_object, keys=["y"]) and \
		   json_exists_and_not_null(json_vector_object, keys=["z"]):
			self.set_reference("world")
			if self.c0.set_parameter_from_key(json_vector_object, ["x"]) and \
			   self.c1.set_parameter_from_key(json_vector_object, ["y"]) and \
			   self.c2.set_parameter_from_key(json_vector_object, ["z"]):
				return True
		elif json_exists_and_not_null(json_vector_object, keys=["u"]) and \
		     json_exists_and_not_null(json_vector_object, keys=["v"]) and \
		     json_exists_and_not_null(json_vector_object, keys=["w"]):
			self.set_reference("local")
			if self.c0.set_parameter_from_key(json_vector_object, ["u"]) and \
			   self.c1.set_parameter_from_key(json_vector_object, ["v"]) and \
			   self.c2.set_parameter_from_key(json_vector_object, ["w"]):
				return True
		elif json_exists_and_not_null(json_vector_object, keys=["r"]) and \
		     json_exists_and_not_null(json_vector_object, keys=["s"]) and \
		     json_exists_and_not_null(json_vector_object, keys=["t"]):
			self.set_reference("sample")
			if self.c0.set_parameter_from_key(json_vector_object, ["r"]) and \
			   self.c1.set_parameter_from_key(json_vector_object, ["s"]) and \
			   self.c2.set_parameter_from_key(json_vector_object, ["t"]):
				return True

		raise Exception("Unable to set scene vector from JSON file. A vector must be specified by the three components (x,y,z), (u,v,w) or (r,s,t).")
		return False

	def json_dict(self) -> dict:
		"""Create a dictionary of this scenevector for a CTSimU JSON file."""
		jd = dict()

		if self.reference == "world":
			jd["x"] = self.c0.json_dict()
			jd["y"] = self.c1.json_dict()
			jd["z"] = self.c2.json_dict()
		elif self.reference == "local":
			jd["u"] = self.c0.json_dict()
			jd["v"] = self.c1.json_dict()
			jd["w"] = self.c2.json_dict()
		elif self.reference == "sample":
			jd["r"] = self.c0.json_dict()
			jd["s"] = self.c1.json_dict()
			jd["t"] = self.c2.json_dict()

		return jd