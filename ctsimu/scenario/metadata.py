# -*- coding: UTF-8 -*-
"""
CTSimU metadata of a CT scan.
"""

from ..helpers import *
from .group import Group
from .parameter import Parameter

class Metadata(Group):
	"""CTSimU metadata properties."""
	def __init__(self, _root=None):
		Group.__init__(self, name="metadata", _root=_root)

		self.file = Group("file")
		self.file.set(key="name", value="", native_unit="string", simple=True)
		self.file.set(key="description", value="", native_unit="string", simple=True)
		self.file.set(key="contact", value="", native_unit="string", simple=True)
		self.file.set(key="date_created", value="", native_unit="string", simple=True)
		self.file.set(key="date_changed", value="", native_unit="string", simple=True)

		self.file_version = Group("version")
		self.file_version.set(key="major", value=1, simple=True)
		self.file_version.set(key="minor", value=0, simple=True)
		self.file.add_subgroup(self.file_version)

		self.add_subgroup(self.file)


		self.output = Group("output")
		self.output.set(key="system", value="", native_unit="string", simple=True)
		self.output.set(key="date_measured", value="", native_unit="string", simple=True)

		# Projections
		self.projections = Group("projections")
		self.projections.set(key="number", value=0, simple=True)
		self.projections.set(key="filename", value="", native_unit="string", simple=True)
		self.projections.set(key="datatype", value="float32", native_unit="string", simple=True)
		self.projections.set(key="byteorder", value="little", native_unit="string", simple=True)

		self.projections_headersize = Group("headersize")
		self.projections_headersize.set(key="file", value=0, simple=True)
		self.projections_headersize.set(key="image", value=0, simple=True)
		self.projections.add_subgroup(self.projections_headersize)

		self.projections_dimensions = Group("dimensions")
		self.projections_dimensions.set(key="x", native_unit="px", value=1000, simple=False)
		self.projections_dimensions.set(key="y", native_unit="px", value=1000, simple=False)
		self.projections.add_subgroup(self.projections_dimensions)

		self.projections_pixelsize = Group("pixelsize")
		self.projections_pixelsize.set(key="x", native_unit="mm", value=0.1, simple=False)
		self.projections_pixelsize.set(key="y", native_unit="mm", value=0.1, simple=False)
		self.projections.add_subgroup(self.projections_pixelsize)

		self.projections_dark_field = Group("dark_field")
		self.projections_dark_field.set(key="number", value=0, simple=True)
		self.projections_dark_field.set(key="frame_average", value=None, simple=True)
		self.projections_dark_field.set(key="filename", value=None, native_unit="string", simple=True)
		self.projections_dark_field.set(key="projections_corrected", value=False, native_unit="bool", simple=True)
		self.projections.add_subgroup(self.projections_dark_field)

		self.projections_flat_field = Group("flat_field")
		self.projections_flat_field.set(key="number", value=0, simple=True)
		self.projections_flat_field.set(key="frame_average", value=None, simple=True)
		self.projections_flat_field.set(key="filename", value=None, native_unit="string", simple=True)
		self.projections_flat_field.set(key="projections_corrected", value=False, native_unit="bool", simple=True)
		self.projections.add_subgroup(self.projections_flat_field)

		self.output.add_subgroup(self.projections)

		# Tomogram
		self.tomogram = Group("tomogram")
		self.tomogram.set(key="filename", value="", native_unit="string", simple=True)
		self.tomogram.set(key="datatype", value="float32", native_unit="string", simple=True)
		self.tomogram.set(key="byteorder", value="little", native_unit="string", simple=True)

		self.tomogram_headersize = Group("headersize")
		self.tomogram_headersize.set(key="file", value=0, simple=True)
		self.tomogram_headersize.set(key="image", value=0, simple=True)
		self.tomogram.add_subgroup(self.tomogram_headersize)

		self.tomogram_dimensions = Group("dimensions")
		self.tomogram_dimensions.set(key="x", native_unit="px", value=1000, simple=False)
		self.tomogram_dimensions.set(key="y", native_unit="px", value=1000, simple=False)
		self.tomogram_dimensions.set(key="z", native_unit="px", value=1000, simple=False)
		self.tomogram.add_subgroup(self.tomogram_dimensions)

		self.tomogram_voxelsize = Group("voxelsize")
		self.tomogram_voxelsize.set(key="x", native_unit="mm", value=0.02, simple=False)
		self.tomogram_voxelsize.set(key="y", native_unit="mm", value=0.02, simple=False)
		self.tomogram_voxelsize.set(key="z", native_unit="mm", value=0.02, simple=False)
		self.tomogram.add_subgroup(self.tomogram_voxelsize)

		self.output.add_subgroup(self.tomogram)

		# Reconstruction
		self.reconstruction = Group("reconstruction")
		self.reconstruction.set(key="software", value="", native_unit="string", simple=True)
		self.reconstruction_settings = Group("settings")
		self.reconstruction.add_subgroup(self.reconstruction_settings)

		self.output.add_subgroup(self.reconstruction)

		self.acquisitionGeometry = Group("acquisitionGeometry")
		self.acquisitionGeometry.set(key="path_to_CTSimU_JSON", value="", native_unit="string", simple=True)
		self.output.add_subgroup(self.acquisitionGeometry)

		self.add_subgroup(self.output)