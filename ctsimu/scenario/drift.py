# -*- coding: UTF-8 -*-
"""
Drift structure for a `ctsimu.scenario.parameter.Parameter`.
"""

import numbers
import math
import numpy
from ..helpers import *

class Drift:
	"""Drift structure for a parameter.
	Drifts typically belong to a `ctsimu.scenario.parameter.Parameter` object,
	which keeps a list of possible parameter drifts.

	Attributes
	----------
	trajectory : list
		List of drift values. If the number of drift values does not
		match the number of frames in the CT scan, a linear interpolation
		between the drift values can be used to calculate the drift value
		for a specific frame.

	native_unit : str
		The drift's native unit. Should match the native unit of the
		`ctsimu.scenario.parameter.Parameter` object to which this drift belongs.
		Possible values: `None`, `"mm"`, `"rad"`, `"deg"`, `"s"`, `"mA"`, `"kV"`, `"g/cm^3"`, `"lp/mm"`, `"deg/s"`, `"C"`, `"bool"`, `"string"`.

	known_to_reconstruction : bool
		Should this parameter drift be considered by the reconstruction software?
		This attribute is obeyed when calculating projection matrices.

	interpolation : bool
		Does a linear interpolation take place between the components of the
		`trajectory` list if the number of trajectory components does not match
		the number of frames in the CT scan?
	"""

	def __init__(self, native_unit:str, preferred_unit:str=None, _root=None):
		"""A new drift object must be assigned a valid native unit
		to enable the JSON parser to convert the drift values from the
		JSON file, if necessary.

		Parameters
		----------
		native_unit : str
			The drift's native unit. Should match the native unit of the
			`ctsimu.scenario.parameter.Parameter` object to which this drift belongs.
			Possible values: `None`, `"mm"`, `"rad"`, `"deg"`, `"s"`, `"mA"`, `"kV"`, `"g/cm^3"`, `"lp/mm"`, `"deg/s"`, `"C"`, `"bool"`, `"string"`.

		preferred_unit : str
			Preferred unit to represent these drift values in the JSON file.
			If set to `None`, the native unit will be used.
		"""
		self._root = _root  # root scenario object

		self.value = []
		self.native_unit = native_unit

		if preferred_unit is not None:
			self.preferred_unit = preferred_unit
		else:
			self.preferred_unit = native_unit

		self.known_to_reconstruction = True
		self.interpolation = True

		self.reset()

	def reset(self):
		"""Clear the list of drift values, set `known_to_reconstruction` to `True`."""
		self.known_to_reconstruction = True
		self.value = []

	def json_dict(self) -> dict:
		"""Create a CTSimU JSON dictionary for this drift.

		Returns
		-------
		json_dict : dict
		"""
		jd = dict()
		if len(self.value) > 0:
			jd["value"] = self.value
		else:
			jd["value"] = None

		if self.native_unit not in native_units_to_omit_in_json_file:
			unit = self.native_unit
			if self.preferred_unit is not None:
				unit = self.preferred_unit

				# Convert value to preferred unit:
				conversion_factor = convert_to_native_unit(given_unit=self.preferred_unit, native_unit=self.native_unit, value=1)

				if conversion_factor != 0:
					jd["value"] = list()
					for val in self.value:
						jd["value"].append(val/conversion_factor)

			jd["unit"] = unit

		jd["known_to_reconstruction"] = self.known_to_reconstruction
		return jd

	def set_known_to_reconstruction(self, known:bool):
		"""Set the `known_to_reconstruction` property.

		Parameters
		----------
		known : bool
			Should the reconstruction software take this drift into account?
		"""
		self.known_to_reconstruction = known

	def set_interpolation(self, intpol:bool):
		"""Set the `interpolation` property. Only necessary for
		number-type parameters. String parameters cannot run an interpolation.

		Parameters
		----------
		intpol : bool
			`True` if the drift value should be interpolated if the number of
			drift components does not match the number of frames in the CT scan.
			`False` if no interpolation should take place: the last drift value
			that would match the given frame number is taken.
		"""
		self.interpolation = intpol

	def set_native_unit(self, native_unit:str):
		"""Set the drifts's native unit.

		Parameters
		----------
		native_unit : str
			New native unit for the drift.
			Possible values: `None`, `"mm"`, `"rad"`, `"deg"`, `"s"`, `"mA"`, `"kV"`, `"g/cm^3"`, `"lp/mm"`, `"deg/s"`, `"C"`, `"bool"`, `"string"`.
		"""
		if is_valid_native_unit(native_unit):
			self.native_unit = native_unit

			if native_unit == "string":
				self.interpolation = False

	def add_drift_component(self, value):
		"""Add a drift component to the list of drift values.

		Parameters
		----------
		value : float or str or bool
			Drift value. For number-type drifts (`float`), this is given as an
			absolute deviation from the parameter's standard value.
			For string-type drifts (`str`), this is the string that should replace
			the parameter's standard value.
		"""
		self.value.append(value)

	def set_from_json(self, json_drift_object:dict) -> bool:
		"""Set the drift from a given CTSimU drift structure.
		The proper `native_unit` must be set up correctly before
		running this function.

		Parameters
		----------
		json_drift_object : dict
			A CTSimU drift structure, imported from a JSON structure.

		Returns
		-------
		success : bool
			`True` on success, `False` if an error occurred.
		"""
		self.reset()

		# Get JSON unit:
		if json_exists_and_not_null(json_drift_object, ["unit"]):
			self.preferred_unit = get_value(
				dictionary=json_drift_object,
				keys=["unit"],
				fail_value=self.native_unit
			)

		# Known to reconstruction
		if json_exists_and_not_null(json_drift_object, ["known_to_reconstruction"]):
			self.set_known_to_reconstruction(
				get_value_in_native_unit(
					native_unit="bool",
					dictionary=json_drift_object,
					keys=["known_to_reconstruction"],
					fail_value=True
				)
			)

		# Get drift values:
		if json_exists_and_not_null(json_drift_object, ["value"]):
			value = json_extract(json_drift_object, ["value"])
			if isinstance(value, list):
				# Drift value is an array.
				for v in value:
					self.add_drift_component(
						convert_to_native_unit(self.preferred_unit, self.native_unit, v)
					)

				return True
			elif isinstance(value, numbers.Number):
				# Interpret value as a single number:
				self.add_drift_component(
					convert_to_native_unit(self.preferred_unit, self.native_unit, value)
				)
				return True
		elif json_exists_and_not_null(json_drift_object, ["file"]):
			filename = json_extract(json_drift_object, ["file"])
			if isinstance(filename, str):
				# Read drift values from a file
				if self._root is not None:
					# Get path relative to JSON file:
					rel_filename = self._root.path_of_external_file(filename)
				else:
					rel_filename = filename

				values = read_csv_file(rel_filename)

				firstColumn = values[0]
				for v in firstColumn:
					self.add_drift_component(float(v))

				return True

		return False

	def get_value_for_frame(self, frame:float, n_frames:int) -> float | str | bool:
		"""Return a drift value for the given `frame` number,
		assuming a total number of `n_frames` in the CT scan.
		If interpolation is activated, linear interpolation will
		take place between drift values, but also for frame numbers
		outside of the expected range: (frame < 0) and (frame >= n_frames).
		Note that the frame number starts at 0.

		Parameters
		----------
		frame : float
			Current frame number.

		n_frames : int
			Total number of frames in scan.

		Returns
		-------
		value : float or str or bool
			The drift's value for the given `frame` number.
		"""

		nTrajectoryPoints = len(self.value)
		if nTrajectoryPoints is not None:
			if nTrajectoryPoints > 1:
				if n_frames > 1:
					# We know that we have at least two drift values, so we are
					# on the safe side for linear interpolations.
					# We also know that we have at least two frames in the scan,
					# and can therefore map our "scan progress" (the current
					# frame number) to the array of drift values.

					lastFrameNr = n_frames - 1 # frames start counting at 0

					# Frame progress on a scale between 0 and 1.
					# 0: first frame (start), 1: last frame (finish).
					progress = float(frame) / float(lastFrameNr)

					# Calculate the theoretical array index to get the drift
					# value for the current frame from the array of drift values:
					lastTrajectoryIndex = int(nTrajectoryPoints - 1)
					trajectoryIndex = progress * float(lastTrajectoryIndex)

					if (progress >= 0) and (progress <= 1):
						# We are inside the array of drift values.
						leftIndex = int(math.floor(trajectoryIndex))

						if float(leftIndex) == float(trajectoryIndex):
							# We are exactly at one trajectory point;
							# no need for interpolation.
							return self.value[leftIndex]

						if self.interpolation is True:
							# Linear interpolation:
							rightIndex = int(leftIndex + 1)

							# We return a weighted average of the two drift
							# values where the current frame is "in between".
							# Weight for the right bin is frac(trajectoryIndex).
							rightWeight = float(trajectoryIndex) - float(math.floor(trajectoryIndex))

							# Weight for the left bin is 1-rightWeight.
							leftWeight = 1.0 - rightWeight

							return leftWeight*self.value[leftIndex] \
								+ rightWeight*self.value[rightIndex]
						else:
							# No interpolation.
							# Return the value at the last drift value index
							# that would apply to this frame position.
							return self.value[leftIndex]
					else:
						# Linear interpolation beyond provided trajectory data.
						if progress > 1.0:
							# We are beyond the expected last frame.
							if self.interpolation is True:
								# Linear interpolation beyond last two drift values:
								trajectoryValue0 = self.value[int(lastTrajectoryIndex-1)]
								trajectoryValue1 = self.value[int(lastTrajectoryIndex)]
								offsetValue = trajectoryValue1

								# We assume a linear interpolation function beyond the two
								# last drift values. Taking the last frame as the zero point
								# (i.e., the starting point) of this linear interpolation,
								# the frame's position on the x axis would be:
								xTraj = trajectoryIndex - float(lastTrajectoryIndex)
							else:
								# No interpolation. Return last trajectory value:
								return self.value[lastTrajectoryIndex]
						else:
							# We are before the first frame (i.e., before frame 0).
							if self.interpolation is True:
								trajectoryValue0 = self.value[0]
								trajectoryValue1 = self.value[1]
								offsetValue = trajectoryValue0

								# We assume a linear interpolation function beyond the two
								# last drift values. Taking the last frame as the zero point
								# (i.e., the starting point) of this linear interpolation,
								# the frame's position on the x axis would be:
								xTraj = trajectoryIndex # is negative in this case
							else:
								# No interpolation. Return first trajectory value:
								return self.value[0]

						# m = the slope of our linear interpolation function.
						# We are on the axis of trajectory drift values:
						# Our step size in x direction is 1 (trajectoryValue1 and
						# trajectoryValue0 are one trajectory step apart). Not
						# to be confused with the frame number.
						m = float(trajectoryValue1 - trajectoryValue0)
						driftValue = m * xTraj + offsetValue

						return driftValue
				else:
					# n_frames <= 1
					# If "scan" only has 1 or 0 frames,
					# simply return the first trajectory value.
					return self.value[0]
			else:
				# trajectory points <= 1:
				if nTrajectoryPoints > 0:
					# Simply check if the trajectory consists of at least
					# 1 point and return this value:
					return self.value[0]

		# Drifts are per-frame absolute deviations from the standard value,
		# so 0 is a sane default value for no drift components:
		return 0