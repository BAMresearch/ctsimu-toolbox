# -*- coding: UTF-8 -*-
"""
Parts are objects in the scene: detector, source, stage and samples.
"""

from ..helpers import *
from ..geometry import *
from .group import Group
from .deviation import Deviation
from .parameter import Parameter
from .scenevector import Scenevector

class Part(Group):
	"""Parts are objects in the scene: detector, source, stage and samples.

	They have a coordinate system and can define deviations from their
	standard geometries (translations and rotations around given axes).
	The center, vectors and deviations can all have drifts,
	allowing for an evolution in time.

	Attributes
	----------
	name : str
		Name of the object.

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
	def __init__(self, name:str="Part", _root=None):
		"""The part's name can be set upon initialization.

		Parameters
		----------
		name : str
			Name for the part.
		"""
		Group.__init__(self, _root=_root)
		self.name = name

		self.attached_to_stage = False
		self.coordinate_system = CoordinateSystem()
		self.center = Scenevector(native_unit="mm", _root=self._root)
		self.u = Scenevector(native_unit=None, _root=self._root)
		self.w = Scenevector(native_unit=None, _root=self._root)
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
		"""Reset to standard conditions, delete all deviations, reset
		coordinate system to standard world coordinates,
		and reset all parameters to their standard values."""
		self.attached_to_stage = False
		self._static = False
		self._cs_initialized_real = False
		self._cs_initialized_recon = False

		self.coordinate_system.reset()

		self.deviations = list()
		self.legacy_deviations = list()

		Group.reset(self)

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

	def set_geometry(self, json_geometry_object:dict, stage_coordinate_system:'CoordinateSystem'=None, proper_cs:str="local") -> bool:
		"""
		Set up the part from a CTSimU JSON geometry definition.
		The `stage_coordinate_system` must only be provided if this
		part is attached to the stage.

		Parameters
		----------
		json_geometry_object : dict
			A CTSimU geometry object, as imported from a JSON structure.

		stage_coordinate_system : CoordinateSystem
			Stage coordinate system. Only necessary for samples attached
			to the stage.

		proper_cs : str
			Which is the proper coordinate system of this object?
			Either `"world"` (x, y, z), `"local"` (u, v, w) or `"sample"` (r, s, t).

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
					new_deviation = Deviation(pivot_reference=proper_cs, _root=self._root)
					if new_deviation.set_from_json(dev):
						self.deviations.append(new_deviation)
					else:
						raise Exception(f"An error occurred when setting a geometrical deviation (translation) for part '{self.name}'.")
						return False
			elif isinstance(devs, dict):
				# Only one drift defined directly as a JSON object?
				# Actually not supported by file format,
				# but let's be generous and try...
				new_deviation = Deviation(pivot_reference=proper_cs, _root=self._root)
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
					pos_dev = Deviation(pivot_reference=proper_cs, _root=self._root)
					pos_dev.set_type("translation")
					pos_dev.set_axis(axis)
					pos_dev.set_known_to_reconstruction(known_to_recon)
					if pos_dev.amount.set_from_json(json_extract(geo, ["deviation", "position", axis])):
						# Legacy_deviations not necessary here
						# because positional translations are fully
						# compatible with the new file format:
						if not pos_dev.amount.is_zero():
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
					rot_dev = Deviation(pivot_reference=proper_cs, _root=self._root)
					rot_dev.set_type("rotation")

					# Prior to 0.9, all deviations were meant to take place
					# before the stage rotation. This means they need to be
					# stored as scene vectors to designate a constant deviation axis.
					rot_dev.set_axis(axis)
					rot_dev.set_known_to_reconstruction(known_to_recon)
					if rot_dev.amount.set_from_json(json_extract(geo, ["deviation", "rotation", axis])):
						# success
						if not rot_dev.amount.is_zero():
							self.legacy_deviations.append(rot_dev)
					else:
						raise Exception(f"An error occurred when setting a geometrical deviation (rotation) for part '{self.name}'.")
						return False

		self.set_frame(frame=0, nFrames=1, w_rotation=0, stage_coordinate_system=stage_coordinate_system)
		return True

	def geometry_dict(self) -> dict:
		"""Create a dictionary of the geometry for a CTSimU scenario file.

		Returns
		-------
		geometry_dict : dict
			Dictionary with the geometry of this part.
		"""
		jd = dict()
		jd["center"] = self.center.json_dict()

		if self.u.reference == "world":
			jd["vector_u"] = self.u.json_dict()
		elif self.u.reference == "local":
			jd["vector_r"] = self.u.json_dict()

		if self.w.reference == "world":
			jd["vector_w"] = self.w.json_dict()
		elif self.w.reference == "local":
			jd["vector_t"] = self.w.json_dict()

		# Assemble deviations and legacy deviations:
		devs = None
		if isinstance(self.deviations, list) and isinstance(self.legacy_deviations, list):
			devs = self.deviations + self.legacy_deviations
		elif isinstance(self.deviations, list):
			devs = self.deviations
		elif isinstance(self.legacy_deviations, list):
			devs = self.legacy_deviations

		if devs is not None:
			if len(devs) > 0:
				jd["deviations"] = list()
				for dev in devs:
					jd["deviations"].append(dev.json_dict())

		return jd

	def _set_frame_coordinate_system(self, frame:float, nFrames:int, reconstruction:bool=False, w_rotation:float=0, stage_coordinate_system:'CoordinateSystem'=None):
		"""
		Set up the part's current coordinate system such that
		it complies with the `frame` number and all necessary
		drifts and deviations (assuming a total number of `nFrames`).

		This function is used by `set_frame` and `set_frame_for_reconstruction`
		and is usually not called from outside the object.

		Parameters
		----------
		frame : float
			Current frame number.

		nFrames : int
			Total number of frames in scan.

		reconstruction : bool
			Set up the coordinate system as it
			would be presented to the reconstruction software.
			Ignores deviations and drifts that should not
			be considered during the reconstruction.

		w_rotation : float
			An additional rotation around the part's w axis
			for this frame. Used for the sample stage, which
			rotates during a CT scan.

		stage_coordinate_system : ctsimu.geometry.CoordinateSystem, optional
			If this part is attached to the sample stage,
			the stage coordinate system for the given `frame`
			must be passed.
		"""

		# Set up standard coordinate system at frame zero:
		center = self.center.standard_vector()
		u      = self.u.standard_vector()
		w      = self.w.standard_vector()

		self.coordinate_system.make_from_vectors(center, u, w)
		self.coordinate_system.make_unit_coordinate_system()

		# Legacy rotational deviations (prior to file format 1.0)
		# all took place before any stage rotation:
		# ----------------------------------------------------------
		for legacy_dev in self.legacy_deviations:
			self.coordinate_system = legacy_dev.deviate(
				coordinate_system=self.coordinate_system,
				frame=frame,
				nFrames=nFrames,
				reconstruction=reconstruction,
				attached_to_stage=self.attached_to_stage,
				stage_coordinate_system=stage_coordinate_system,
			)

		# Potential stage rotation:
		# ------------------------------------
		# Potential rotation around the w axis (in rad).
		if w_rotation != 0:
			self.coordinate_system.rotate_around_w(angle=w_rotation)

		# Deviations:
		# ------------------------------------
		for dev in self.deviations:
			self.coordinate_system = dev.deviate(
				coordinate_system=self.coordinate_system,
				frame=frame,
				nFrames=nFrames,
				reconstruction=reconstruction,
				attached_to_stage=self.attached_to_stage,
				stage_coordinate_system=stage_coordinate_system,
			)

		# Drifts (center and vector components):
		# -----------------------------------------------
		# Build a translation vector for the center point
		# from the total drift for this frame and apply
		# the translation:
		if self.center.has_drifts():
			 center_drift = self.center.drift_vector(
			 	frame=frame,
			 	nFrames=nFrames,
			 	reconstruction=reconstruction
			 )
			 self.coordinate_system.translate(translation_vector=center_drift)

		if self.u.has_drifts() or self.w.has_drifts():
			new_u = self.u.vector_for_frame(frame, nFrames, reconstruction)
			new_w = self.w.vector_for_frame(frame, nFrames, reconstruction)

			self.coordinate_system.make_from_vectors(
				center=coordinate_system.center,
				u=new_u,
				w=new_w
			)
			self.coordinate_system.make_unit_coordinate_system()

	def set_frame(self, frame:float, nFrames:int, w_rotation:float=0, stage_coordinate_system:'CoordinateSystem'=None, reconstruction:bool=False):
		"""
		Set up the part for the given `frame` number, obeying all
		deviations and drifts.

		Parameters
		----------
		frame : float
			Current frame number.

		nFrames : int
			Total number of frames in scan.

		w_rotation : float
			An additional rotation (in rad) around the part's w axis
			for this frame. Used for the sample stage, which
			rotates during a CT scan.

		stage_coordinate_system : ctsimu.geometry.CoordinateSystem, optional
			If this part is attached to the sample stage,
			the stage coordinate system for the given `frame`
			must be passed.

		reconstruction : bool
			If `True`, set frame as seen by reconstruction software.
			Default: `False`.
		"""

		# Set up the current coordinate system obeying all drifts:
		if (self._cs_initialized_real is False) or (self._static is False):
			self._set_frame_coordinate_system(
				frame=frame,
				nFrames=nFrames,
				reconstruction=reconstruction,
				w_rotation=w_rotation,
				stage_coordinate_system=stage_coordinate_system
			)
			self._cs_initialized_real  = not reconstruction
			self._cs_initialized_recon = reconstruction

		Group.set_frame(self, frame, nFrames, reconstruction)