# -*- coding: UTF-8 -*-
"""
Collection of parameters that describe the scenario file.
"""
from datetime import datetime

from ..helpers import *
from ..geometry import *
from .group import Group
from .parameter import Parameter

class File(Group):
	"""CTSimU file properties."""
	def __init__(self, _root=None):
		Group.__init__(self, name="file", _root=_root)

		now = datetime.now()

		self.set(key="file_type", value="CTSimU Scenario", native_unit="string", simple=True)

		self.set(key="name", value=None, native_unit="string", simple=True)
		self.set(key="description", value=None, native_unit="string", simple=True)
		self.set(key="contact", value=None, native_unit="string", simple=True)
		self.set(key="date_created", value=now.strftime("%Y-%m-%d"), native_unit="string", simple=True)
		self.set(key="date_changed", value=now.strftime("%Y-%m-%d"), native_unit="string", simple=True)

		# version
		self.new_subgroup("version")
		self.version.set(key="major", value=None, native_unit=None, simple=True)
		self.version.set(key="minor", value=None, native_unit=None, simple=True)

		# file format version
		self.new_subgroup("file_format_version")
		self.file_format_version.set(key="major", value=1, native_unit=None, simple=True)
		self.file_format_version.set(key="minor", value=2, native_unit=None, simple=True)
