# -*- coding: UTF-8 -*-
"""
Material component: formula and mass fraction.
Material: composition, density, and their drifts
"""
import numbers

from ..helpers import *
from .parameter import Parameter
from .drift import Drift

class MaterialComponent:
	"""Combination of a chemical formula and a mass fraction,
	used as a component for a full material (which may consist
	of a number of material components).

	Attributes
	----------
	formula : ctsimu.scenario.parameter.Parameter
		Chemical formula for the material component.
		A `Parameter` which represents a string (and
		may drift over time).

	mass_fraction : ctsimu.scenario.parameter.Parameter
		The component's mass fraction. A `Parameter` which
		represents a `float` number (and may drift over time).

	parent_material_id : str
		Identifier of the parent material (`material_id`
		as defined in the JSON scenario description).

	parent_material_name : str
		Trivial name of the parent material.
	"""

	def __init__(self, formula:str="Fe", mass_fraction:float=1, parent_material_id:str="", parent_material_name:str=""):
		"""All four attributes can be set when the material
		component is initialized.

		Parameters
		----------
		formula : str
			Chemical formula for the material component.
			A `Parameter` which represents a string (and
			may drift over time).

		mass_fraction : float
			The component's mass fraction as a `float` number.

		parent_material_id : str
			Identifier of the parent material (`material_id`
			as defined in the JSON scenario description).

		parent_material_name : str
			Trivial name of the parent material.
		"""

		self.parent_material_id = parent_material_id
		self.parent_material_name = parent_material_name
		self.formula = Parameter(native_unit="string", standard_value=formula)
		self.mass_fraction = Parameter(native_unit=None, standard_value=mass_fraction)

	def json_dict(self) -> dict:
		"""Create a CTSimU JSON dictionary for this material component.

		Returns
		-------
		json_dict : dict
		"""
		jd = dict()
		jd["formula"] = self.formula.json_dict()
		jd["mass_fraction"] = self.mass_fraction.json_dict()
		return jd

	def set_frame(self, frame:float, nFrames:int, reconstruction:bool=False) -> bool:
		"""Prepare the material component for the given `frame` number,
		considering potential drifts.

		Parameters
		----------
		frame : float
			Current frame number.

		nFrames : int
			Total number of frames in scan.

		reconstruction : bool
			`True` if only those drifts known to the reconstruction software
			must be taken into account, `False` if all drifts are taken
			into account.

		Returns
		-------
		changed : bool
			If the material component has changed from its previous state
			(due to drifts).
		"""

		value_changed = self.formula.set_frame(frame, nFrames, reconstruction) or self.mass_fraction.set_frame(frame, nFrames, reconstruction)
		return value_changed

	def set_from_json(self, json_object:dict):
		"""Set the material component from a CTSimU JSON object.

		Parameters
		----------
		json_object : dict
			CTSimU material component definition, as imported
			from a JSON file.

		Raises
		------
		ValueError
			If formula or mass fraction cannot be read from the `json_object`.
		"""

		if not self.formula.set_parameter_from_key(json_object, ["formula"]):
			raise ValueError(f"CTSimU MaterialComponent: set_from_json: Error reading a formula for material '{self.parent_material_id} ({self.parent_material_name})'.")

		if not self.mass_fraction.set_parameter_from_key(json_object, ["mass_fraction"]):
			raise ValueError(f"CTSimU MaterialComponent: set_from_json: Error reading a mass fraction for material '{self.parent_material_id} ({self.parent_material_name})'.")

		self.set_frame(0, 1, False)

	def set_from_json_legacy(self, json_object:dict):
		"""Set the material component from a legacy CTSimU JSON object
		(for file format version <= 1.0)

		Parameters
		----------
		json_object : dict
			CTSimU material component definition, as imported
			from a JSON file.

		Raises
		------
		ValueError
			If composition cannot be read from the `json_object`.
		"""

		# In file format version <= 1.0, the composition was
		# simply a string value, no mass fraction was defined.
		if not self.formula.set_parameter_from_key(json_object, ["composition"]):
			raise ValueError(f"CTSimU MaterialComponent: set_from_json_legacy: Error reading a composition for material '{self.parent_material_id} ({self.parent_material_name})'.")

		self.mass_fraction.reset()
		self.mass_fraction.set_standard_value(1.0)

		self.set_frame(0, 1, False)


class Material:
	"""Sample material.

	Attributes
	----------
	material_id : str
		Material ID, as defined in the CTSimU scenario description to
		reference the material definition.

	name : str
		Trivial name for the material.

	density : ctsimu.scenario.parameter.Parameter
		Mass density in g/cm³.

	composition : list
		List of `MaterialComponent` objects that define
		the composition (chemical formula and mass fracion) of
		the material.
	"""

	def __init__(self, material_id:str=None, name:str=None, density:float=0):
		"""The init constructor takes the following optional arguments.

		Parameters
		----------
		material_id : str
			Material ID, as defined in the CTSimU scenario description to
			reference the material definition.

		name : str
			Trivial name for the material.

		density : float
			Simple floating point value to define the mass density (in g/cm³).

		Returns
		-------
		changed : bool
			If the material has changed from its previous state
			(due to drifts).
		"""

		self.material_id = material_id
		self.name = name
		self.density = Parameter(native_unit="g/cm^3", standard_value=density)
		self.composition = list()

	def reset(self):
		"""Reset density to `0` (with no drifts),
		clear all material composition entries."""
		self.density.reset()
		self.density.set_standard_value(0.0)
		self.composition = list()

	def json_dict(self) -> dict:
		"""Create a CTSimU JSON dictionary for this material.

		Returns
		-------
		json_dict : dict
		"""
		jd = dict()
		jd["id"] = self.material_id
		jd["name"] = self.name
		jd["density"] = self.density.json_dict()
		jd["composition"] = list()
		for component in self.composition:
			jd["composition"].append(component.json_dict())

		return jd

	def set_frame(self, frame:float, nFrames:int, reconstruction:bool=False) -> bool:
		"""Prepare the material for the given `frame` number,
		considering potential drifts.

		Parameters
		----------
		frame : float
			Current frame number.

		nFrames : int
			Total number of frames in scan.
		"""

		density_changed = self.density.set_frame(frame, nFrames, reconstruction)

		composition_changed = False
		for comp in self.composition:
			composition_changed = composition_changed or comp.set_frame(frame, nFrames, reconstruction)

		changed = density_changed or composition_changed
		return changed

	def set_density(self, density):
		"""Set the density.

		Parameters
		----------
		density : float or ctsimu.scenario.parameter.Parameter
			The density (in g/cm³). Given as a simple numerical
			value (`float`), it will be converted to a Parameter
			with this as the standard value.
			It can also be given directly as a `Parameter` object.

		Raises
		------
		TypeError
			If the `density` is neither of type `float` nor an instance
			of `ctsimu.scenario.parameter.Parameter`.
		"""
		if isinstance(density, Parameter):
			# given as Parameter object
			self.density = density
		elif isinstance(density, numbers.Number):
			# given as float
			self.density.reset()
			self.density.set_standard_value(density)
		else:
			raise TypeError(f"CTSimU Material: set_density: invalid type of given density: '{type(density)}'. Must be 'float' or 'Parameter'.")

	def add_component(self, component:'MaterialComponent'):
		"""Add a material component to the material's components list.

		Parameters
		----------
		component : MaterialComponent
			Material component object.

		Raises
		------
		TypeError
			If the composition is not an instance
			of `ctsimu.scenario.material.MaterialComponent`.
		"""
		if isinstance(component, MaterialComponent):
			self.composition.append(component)
		else:
			raise TypeError(f"CTSimU Material: add_component: invalid type of given component: '{type(component)}'. Only a 'MaterialComponent' object will be accepted.")

	def set_from_json(self, json_object:dict):
		"""Set up the material from a CTSimU material definition.

		Parameters
		----------
		json_object : dict
			CTSimU material definition, as imported	from a JSON file.

		Raises
		------
		ValueError
			If the required properties are missing or not set in the
			`json_object` dictionary.
		"""

		self.reset()

		if json_exists_and_not_null(json_object, ["id"]):
			self.material_id = get_value(json_object, ["id"], fail_value=None)
		else:
			raise ValueError("CTSimU Material: set_from_json: missing id in material definition.")

		self.name = get_value(json_object, ["name"], fail_value=None)

		if not self.density.set_parameter_from_key(json_object, ["density"]):
			raise ValueError(f"CTSimU Material: set_from_json: error reading density of material '{self.material_id} ({self.name})'.")

		if json_exists_and_not_null(json_object, ["composition"]):
			# The composition should be an array since file format version 1.1:
			json_components = json_extract(json_object, ["composition"])
			if isinstance(json_components, list):
				for json_component in json_components:
					new_component = MaterialComponent(
						parent_material_id=self.material_id,
						parent_material_name=self.name
					)
					new_component.set_from_json(json_component)
					self.add_component(new_component)
			else:
				# Probably legacy definition?
				# Try composition given as single string.
				new_component = MaterialComponent(
					parent_material_id=self.material_id,
					parent_material_name=self.name
				)
				new_component.set_from_json_legacy(json_object)
				self.add_component(new_component)
		else:
			raise ValueError(f"CTSimU Material: set_from_json: error reading composition of material '{self.material_id} ({self.name})'.")

		self.set_frame(0, 1, False)