# -*- coding: UTF-8 -*-
"""
Collection of parameters that describe the scenario file.
"""

from ..helpers import *
from ..geometry import *
from .group import Group
from .parameter import Parameter

class File(Group):
	"""CTSimU file properties."""
	def __init__(self):
		Group.__init__(self, "file")

		self.set(key="name", value=None, native_unit="string", simple=True)
		self.set(key="description", value=None, native_unit="string", simple=True)
		self.set(key="contact", value=None, native_unit="string", simple=True)
		self.set(key="date_created", value=None, native_unit="string", simple=True)
		self.set(key="date_changed", value=None, native_unit="string", simple=True)

		# version
		self.version = Group(name="version")
		self.version.set(key="major", value=None, native_unit=None, simple=True)
		self.version.set(key="minor", value=None, native_unit=None, simple=True)
		self.add_subgroup(self.version)

		# file format version
		self.file_format_version = Group(name="file_format_version")
		self.file_format_version.set(key="major", value=1, native_unit=None, simple=True)
		self.file_format_version.set(key="minor", value=2, native_unit=None, simple=True)
		self.add_subgroup(self.file_format_version)
