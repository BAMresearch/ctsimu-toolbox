# -*- coding: UTF-8 -*-

class Drift:
	def __init__(self, native_unit):
		self.known_to_reconstruction = True
		self.interpolation = True
		self.trajectory = []
		self.native_unit = native_unit

		self.reset()

	def reset(self):
		self.known_to_reconstruction = True
		self.trajectory = []

	def setFromJSON(self, drift):
		self.reset()

