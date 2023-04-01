# -*- coding: UTF-8 -*-

import numbers
from ..helpers import *

class Drift:
	def __init__(self, native_unit:str):
		self.known_to_reconstruction = True
		self.interpolation = True
		self.trajectory = []
		self.native_unit = native_unit

		self.reset()

	def reset(self):
		self.known_to_reconstruction = True
		self.trajectory = []

	def set_known_to_reconstruction(self, known:bool):
		self.known_to_reconstruction = known

	def set_interpolation(self, intpol:bool):
		self.interpolation = intpol

	def set_native_unit(self, native_unit:str):
		self.native_unit = native_unit

		if native_unit == "string":
			self.interpolation = False

	def set_from_json(self, json_object:dict):
		"""Set the drift from a given JSON drift object."""
		self.reset()

		# Get JSON unit:
		jsonUnit = self.native_unit
		if json_exists_and_not_null(json_object, ["unit"]):
			jsonUnit = get_value(dictionary=json_object, keys=["unit"], fail_value=self.native_unit)

		# Known to reconstruction
		if json_exists_and_not_null(json_object, ["known_to_reconstruction"]):
			self.set_known_to_reconstruction(get_value_in_unit(native_unit="bool", dictionary=json_object, keys=["known_to_reconstruction"], fail_value=True))

		# Get drift values:
		if json_exists_and_not_null(json_object, ["value"]):
			value = json_extract(json_object, ["value"])
			if isinstance(value, list):
				# Drift value is an array.
				for v in value:
					self.trajectory.append(convert_to_native_unit(jsonUnit, self.native_unit, v))

				return True
			elif isinstance(value, numbers.Number)::
				# Interpret value as a single number:
				self.trajectory.append(convert_to_native_unit(jsonUnit, self.native_unit, value))
				return True
		elif json_exists_and_not_null(json_object, ["file"]):
			filename = json_extract(json_object, ["file"])
			if isinstance(filename, str):
				# Read drift values from a file
				values = read_csv_file(filename)
				firstColumn = values[0]
				for v in firstColumn:
					self.trajectory.append(float(v))

				return True

		return False

	def get_value_for_frame(self, frame:int, nFrames:int):
		...