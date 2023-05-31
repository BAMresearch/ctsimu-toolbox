# -*- coding: UTF-8 -*-
"""
Collection of parameters that describe the environment.
"""

from ..helpers import *
from ..geometry import *
from .group import Group
from .parameter import Parameter

class Environment(Group):
	"""CTSimU environment properties."""
	def __init__(self):
		Group.__init__(self, "environment")

		self.set(key="material_id", value=None, native_unit="string", simple=True)
		self.set(key="temperature", value=None, native_unit="C")