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

        self.new_subgroup("file")
        self.file.set(key="name", value="", native_unit="string", simple=True)
        self.file.set(key="description", value="", native_unit="string", simple=True)
        self.file.set(key="contact", value="", native_unit="string", simple=True)
        self.file.set(key="date_created", value="", native_unit="string", simple=True)
        self.file.set(key="date_changed", value="", native_unit="string", simple=True)
        self.file.set(key="file_type", value="CTSimU Metadata", native_unit="string", simple=True)

        self.file.new_subgroup("file_format_version")
        self.file.file_format_version.set(key="major", value=ctsimu_supported_metadata_version["major"], native_unit=None, simple=True)
        self.file.file_format_version.set(key="minor", value=ctsimu_supported_metadata_version["minor"], native_unit=None, simple=True)

        self.new_subgroup("output")
        self.output.set(key="system", value="", native_unit="string", simple=True)
        self.output.set(key="date_measured", value="", native_unit="string", simple=True)

        # Projections
        self.output.new_subgroup("projections")
        self.output.projections.set(key="filename", value="", native_unit="string", simple=True)
        self.output.projections.set(key="number", value=0, native_unit=None, simple=True)
        self.output.projections.set(key="frame_average", value=None, native_unit=None, simple=True)
        self.output.projections.set(key="max_intensity", value=None, native_unit=None, simple=True)
        self.output.projections.set(key="datatype", value="float32", native_unit="string", simple=True)
        self.output.projections.set(key="byteorder", value="little", native_unit="string", simple=True)

        self.output.projections.new_subgroup("headersize")
        self.output.projections.headersize.set(key="file", value=0, native_unit=None, simple=True)
        self.output.projections.headersize.set(key="image", value=0, native_unit=None, simple=True)

        self.output.projections.new_subgroup("dimensions")
        self.output.projections.dimensions.set(key="x", native_unit="px", value=None, simple=False)
        self.output.projections.dimensions.set(key="y", native_unit="px", value=None, simple=False)

        self.output.projections.new_subgroup("pixelsize")
        self.output.projections.pixelsize.set(key="x", native_unit="mm", value=None, simple=False)
        self.output.projections.pixelsize.set(key="y", native_unit="mm", value=None, simple=False)

        self.output.projections.new_subgroup("dark_field")
        self.output.projections.dark_field.set(key="number", native_unit=None, value=0, simple=True)
        self.output.projections.dark_field.set(key="frame_average", value=None, native_unit=None, simple=True)
        self.output.projections.dark_field.set(key="filename", value=None, native_unit="string", simple=True)
        self.output.projections.dark_field.set(key="projections_corrected", value=False, native_unit="bool", simple=True)

        self.output.projections.new_subgroup("flat_field")
        self.output.projections.flat_field.set(key="number", native_unit=None, value=0, simple=True)
        self.output.projections.flat_field.set(key="frame_average", value=None, native_unit=None, simple=True)
        self.output.projections.flat_field.set(key="filename", value=None, native_unit="string", simple=True)
        self.output.projections.flat_field.set(key="projections_corrected", value=False, native_unit="bool", simple=True)

        self.output.projections.new_subgroup("bad_pixel_map")
        self.output.projections.bad_pixel_map.set(key="filename", value=None, native_unit="string", simple=True)
        self.output.projections.bad_pixel_map.set(key="projections_corrected", value=False, native_unit="bool", simple=True)

        # Tomogram
        self.output.new_subgroup("tomogram")
        self.output.tomogram.set(key="filename", value="", native_unit="string", simple=True)
        self.output.tomogram.set(key="datatype", value="float32", native_unit="string", simple=True)
        self.output.tomogram.set(key="byteorder", value="little", native_unit="string", simple=True)

        self.output.tomogram.new_subgroup("headersize")
        self.output.tomogram.headersize.set(key="file", value=0, native_unit=None, simple=True)
        self.output.tomogram.headersize.set(key="image", value=0, native_unit=None, simple=True)

        self.output.tomogram.new_subgroup("dimensions")
        self.output.tomogram.dimensions.set(key="x", native_unit="px", value=None, simple=False)
        self.output.tomogram.dimensions.set(key="y", native_unit="px", value=None, simple=False)
        self.output.tomogram.dimensions.set(key="z", native_unit="px", value=None, simple=False)

        self.output.tomogram.new_subgroup("voxelsize")
        self.output.tomogram.voxelsize.set(key="x", native_unit="mm", value=None, simple=False)
        self.output.tomogram.voxelsize.set(key="y", native_unit="mm", value=None, simple=False)
        self.output.tomogram.voxelsize.set(key="z", native_unit="mm", value=None, simple=False)

        # Acquisition geometry
        self.new_subgroup("acquisition_geometry")
        self.acquisition_geometry.set(key="path_to_CTSimU_JSON", value=None, native_unit="string", simple=True)

        # Reconstruction
        self.new_subgroup("reconstruction")
        self.reconstruction.set(key="software", value=None, native_unit="string", simple=True)
        self.reconstruction.new_subgroup("settings")

        # Simulation
        self.new_subgroup("simulation")

        # aRTist-specific parameters. Not sure if we really need those in the toolbox,
        # as it is simulation-software agnostic.
        self.simulation.set(key="full_simulation", value=None, native_unit="bool", simple=True)
        self.simulation.set(key="compute_detector", value=None, native_unit="bool", simple=True)
        self.simulation.set(key="compute_xray_source", value=None, native_unit="bool", simple=True)
        self.simulation.set(key="load_samples", value=None, native_unit="bool", simple=True)
        self.simulation.set(key="set_multisampling", value=None, native_unit="bool", simple=True)
        self.simulation.set(key="set_scattering", value=None, native_unit="bool", simple=True)

        self.simulation.new_subgroup("multisampling")
        self.simulation.multisampling.set(key="source", value=None, native_unit="string", simple=True)
        self.simulation.multisampling.set(key="detector", value=None, native_unit="string", simple=True)

        self.simulation.new_subgroup("scattering")
        self.simulation.scattering.set(key="on", value=None, native_unit="bool", simple=True)
        self.simulation.scattering.set(key="image_interval", value=None, native_unit=None, simple=True)
        self.simulation.scattering.set(key="photons", value=None, native_unit=None, simple=True)

        self.simulation.new_subgroup("ctsimu_scenario") # Empty group, don't import any full scenario definition again.