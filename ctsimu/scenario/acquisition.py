# -*- coding: UTF-8 -*-
"""
Collection of parameters that describe the scan acquisition.
"""

from ..helpers import *
from .group import Group
from .parameter import Parameter

class Acquisition(Group):
	"""CTSimU acquisition properties."""
	def __init__(self, _root=None):
		Group.__init__(self, name="acquisition", _root=_root)

		self.set(key="start_angle", value=None, native_unit="deg")
		self.set(key="stop_angle",  value=None, native_unit="deg")
		self.set(key="direction",   value="CCW", native_unit="string", simple=True)
		self.set(key="scan_mode",   value="stop+go", native_unit="string", simple=True)
		self.set(key="scan_speed",  value=None, native_unit="deg/s")
		self.set(key="number_of_projections", value=None, native_unit=None, simple=True)
		self.set(key="include_final_angle", value=True, native_unit="bool", simple=True)
		self.set(key="frame_average", value=None, native_unit=None, simple=True)
		self.set(key="scattering", value=False, native_unit="bool", simple=True)

		self.new_subgroup("dark_field")
		self.dark_field.set(key="number", value=None, simple=True)
		self.dark_field.set(key="frame_average", value=None, simple=True)
		self.dark_field.set(key="ideal", value=False, native_unit="bool", simple=True)
		self.dark_field.set(key="correction", value=False, native_unit="bool", simple=True)

		self.new_subgroup("flat_field")
		self.flat_field.set(key="number", value=None, simple=True)
		self.flat_field.set(key="frame_average", value=None, simple=True)
		self.flat_field.set(key="ideal", value=False, native_unit="bool", simple=True)
		self.flat_field.set(key="correction", value=False, native_unit="bool", simple=True)

		self.new_subgroup("pixel_binning")
		self.pixel_binning.set(key="u", value=1, simple=True)
		self.pixel_binning.set(key="v", value=1, simple=True)