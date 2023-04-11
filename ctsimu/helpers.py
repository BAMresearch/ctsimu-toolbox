# -*- coding: UTF-8 -*-

import os
import json
import math
import numpy
from scipy import optimize, fft

def log(message:str):
    """Print an output message.

    Parameters
    ----------
    message : str
        Message to be printed.
    """
    print(message)

# ----------------
# JSON Handling
# ----------------
valid_native_units = [
        None, "mm", "rad", "deg", "s", "mA", "kV", "g/cm^3", "lp/mm", "bool", "string"
    ]

def is_valid_native_unit(native_unit:str) -> bool:
    if native_unit in valid_native_units:
        return True

    raise Exception(f"CTSimU: Not a valid native unit: '{native_unit}'. Valid options are: {valid_native_units}.")
    return False

def read_json_file(filename:str) -> dict:
    """Read a JSON file into a Python dictionary.

    Parameters
    ----------
    filename : str
        Filename of the JSON file.

    Returns
    -------
    jd : dict
        Dictionary representation of the JSON structure.
    """

    return json.load(filename)

def value_is_null(value) -> bool:
    """Check if a specific value is set to `null`."""
    if value is None:
        return True

    return False

def value_is_null_or_zero(value) -> bool:
    """Check if a specific value is set to `null` or 0."""
    if value is not None:
        if value != 0:
            return False

    return True

def object_value_is_null(json_obj:dict) -> bool:
    """Check if a JSON object has a `value` parameter and if this parameter is set to `null`."""
    if not isinstance(json_obj, dict):
        raise AttributeError("object_value_is_null() expects dict for the `json_obj`.")

    if "value" in json_obj:
        return value_is_null(json_obj["value"])

    return True

def object_value_is_null_or_zero(json_obj:dict) -> bool:
    """Check if a JSON object has a `value` parameter and if this parameter is set to `null` or 0."""
    if not isinstance(json_obj, dict):
        raise AttributeError("object_value_is_null_or_zero() expects dict for the `json_obj`.")

    if "value" in json_obj:
        return value_is_null_or_zero(json_obj["value"])

    return True

def get_value_or_none(dictionary:dict, *keys:str) -> float | str | bool | dict | list:
    """Get the dictionary value that is located at the given
    sequence of `keys`. If it cannot be found, return `None`.
    """
    if not isinstance(dictionary, dict):
        raise AttributeError("get_value_or_none() expects dict as first argument for the `dictionary`.")
    if len(keys) == 0:
        raise AttributeError("get_value_or_none() expects a list of keys, given as a sequence of strings.")

    current_element = dictionary

    if keys is None:
        return None
    else:
        for key in keys:
            if current_element != None:
                if key in current_element:
                    current_element = current_element[key]
                else:
                    return None
            else:
                return None

    return current_element

def get_value(dictionary:dict, keys:list=[], fail_value=None) -> float | str | bool | dict | list:
    """Get the specific value of the parameter that is
    located at the given sequence of `keys` in the JSON dictionary.

    If that element cannot be found or is set to `null`,
    return the `fail_value`.
    """
    if not isinstance(dictionary, dict):
        raise AttributeError("get_value() expects dict as first argument for the `dictionary`.")
    if len(keys) == 0:
        raise AttributeError("get_value() expects a list of keys as second argument for the `keys`.")

    result = get_value_or_none(dictionary, *keys)
    if result is None:
        return fail_value

    return result

def json_exists(dictionary:dict, keys:list=[]) -> bool:
    if not isinstance(dictionary, dict):
        raise AttributeError("json_exists() expects dict as first argument for the `dictionary`.")
    if len(keys) == 0:
        raise AttributeError("json_exists() expects a list of keys as second argument for the `keys`.")

    current_element = dictionary
    for key in keys:
        if key in current_element:
            current_element = current_element[key]
        else:
            return False

    return True

def json_isnull(dictionary:dict, keys:list=[]) -> bool:
    if not isinstance(dictionary, dict):
        raise AttributeError("json_isnull() expects dict as first argument for the `dictionary`.")
    if len(keys) == 0:
        raise AttributeError("json_isnull() expects a list of keys as second argument for the `keys`.")

    v = get_value(dictionary, keys)
    if v is None:
        return True

    return False

def json_exists_and_not_null(dictionary:dict, keys:list=[]) -> bool:
    if json_exists(dictionary, keys):
        if not json_isnull(dictionary, keys):
            return True

    return False

def json_extract(dictionary:dict, keys:list=[]) -> float | str | bool | dict | list:
    """Get the JSON sub-object that is located at a given sequence of `keys` in the JSON dictionary."""
    return get_value(dictionary, keys)

def json_extract_from_possible_keys(dictionary:dict, key_lists:list):
    """Searches the JSON object for each key sequence in the given list of key sequences. The first sequence that exists will return an extracted JSON object."""
    for keys in key_lists:
        if json_exists(dictionary, keys):
            return json_extract(dictionary, keys)

    return None

# ------------------------------------
# CTSimU conversions to native units
# ------------------------------------
# Each function supports the allowed units from the CTSimU file format specification.

def in_mm(value:float, unit:str) -> float:
    """Convert a length to mm."""
    if value is not None:
        if unit == "nm": return (value * 1.0e-6)
        if unit == "um": return (value * 1.0e-3)
        if unit == "mm": return value
        if unit == "cm": return (value * 10.0)
        if unit == "dm": return (value * 100.0)
        if unit == "m":  return (value * 1000.0)
    else:
        return None

    raise Exception(f"Not a valid unit of length: '{unit}'.")

def in_rad(value:float, unit:str="deg") -> float:
    """Convert an angle to radians."""
    if value is not None:
        if unit == "deg": return (value * math.pi / 180.0)
        if unit == "rad": return value
    else:
        return None

    raise Exception(f"Not a valid angular unit: '{unit}'.")

def in_deg(value:float, unit:str="rad") -> float:
    """Convert an angle to degrees."""
    if value is not None:
        if unit == "rad": return (value * 180.0 / math.pi)
        if unit == "deg": return value
    else:
        return None

    raise Exception(f"Not a valid angular unit: '{unit}'.")

def in_s(value:float, unit:str) -> float:
    """Convert a time to s."""
    if value is not None:
        if unit == "ms":  return (value * 1.0e-3)
        if unit == "s":   return value
        if unit == "min": return (value * 60.0)
        if unit == "h":   return (value * 3600.0)
    else:
        return None

    raise Exception(f"Not a valid unit of time: '{unit}'.")

def in_mA(value:float, unit:str) -> float:
    """Convert an electric current to mA."""
    if value is not None:
        if unit == "uA": return (value * 1.0e-3)
        if unit == "mA": return value
        if unit == "A":  return (value * 1000.0)
    else:
        return None

    raise Exception(f"Not a valid unit of electric current: '{unit}'.")

def in_kV(value:float, unit:str) -> float:
    """Convert a voltage to kV."""
    if value is not None:
        if unit == "V":  return (value * 1.0e-3)
        if unit == "kV": return value
        if unit == "MV": return (value * 1000.0)
    else:
        return None

    raise Exception(f"Not a valid unit of voltage: '{unit}'.")

def in_g_per_cm3(value:float, unit:str) -> float:
    """Convert a mass density to g/cm³."""
    if value is not None:
        if unit == "kg/m^3": return (value * 1.0e-3)
        if unit == "g/cm^3": return value
    else:
        return None

    raise Exception(f"Not a valid unit of mass density: '{unit}'.")

def in_lp_per_mm(value:float, unit:str) -> float:
    """Convert a resolution to line pairs per millimeter."""
    if value is not None:
        if unit == "lp/mm": return value
        if unit == "lp/cm": return (value * 0.1)
        if unit == "lp/dm": return (value * 0.01)
        if unit == "lp/m" : return (value * 0.001)
    else:
        return None

    raise Exception(f"Not a valid unit of length: '{unit}'.")

def from_bool(value) -> bool:
    """Convert `value` into a true boolean. (For example, `1` becomes `True`, and `0` becomes `False`.)
    """
    if value:
        return True

    return False

def convert_SNR_FWHM(SNR_or_FWHM:float, intensity:float) -> float:
    """Converts between SNR and Gaussian FWHM for a given intensity
    (i.e., more generally, the given distribution's mean value)."""

    return 2.0 * math.sqrt(2.0 * log(2.0)) * float(intensity) / float(SNR_or_FWHM)

def convert_to_native_unit(given_unit:str, native_unit:str, value:float) -> float | str | bool:
    """Check which native unit is requested, convert value accordingly."""
    if native_unit is None:
        return value
    else:
        if native_unit == "string": return value
        if native_unit == "mm":     return in_mm(value, given_unit)
        if native_unit == "s":      return in_s(value, given_unit)
        if native_unit == "deg":    return in_deg(value, given_unit)
        if native_unit == "rad":    return in_rad(value, given_unit)
        if native_unit == "mA":     return in_mA(value, given_unit)
        if native_unit == "kV":     return in_kV(value, given_unit)
        if native_unit == "g/cm^3": return in_g_per_cm3(value, given_unit)
        if native_unit == "lp/mm":  return in_lp_per_mm(value, given_unit)
        if native_unit == "bool":   return from_bool(value)

    raise Exception(f"Native unit '{native_unit}' is incompatible with the given unit '{given_unit}'.")
    return None

def json_convert_to_native_unit(native_unit:str, value_and_unit:dict, fallback_json_unit:str=None) -> float | str | bool:
    """Like the previous function `convert_to_native_unit`, but takes
    a JSON object `value_and_unit`, i.e. a dictionary that must contain a `value` and
    an associated `unit`. Checks which native unit is requested, converts
    JSON `value` accordingly. `fallback_json_unit` is used if the unit is not specified
    in the `value_and_unit` JSON object."""

    if native_unit is None:
        # No native unit given. Simply return the value.
        return get_value(value_and_unit, ["value"])
    elif native_unit == "bool":
        # This is not a value/unit dictionary, but a boolean.
        return from_bool(value_and_unit)
    elif native_unit == "string":
        if json_exists(value_and_unit, ["value"]):
            return get_value(value_and_unit, ["value"])
        else:
            if isinstance(value_and_unit, str):
                return value_and_unit

        raise Exception(f"Given value does not seem to be a string: {value_and_unit}")
    else:
        if json_exists(value_and_unit, ["value"]):
            value = get_value(value_and_unit, ["value"])
            unit = fallback_json_unit
            if json_exists(value_and_unit, ["unit"]):
                # The unit does not necessarily have to exist.
                # For example, in the case of strings it is clear
                # just from the native unit.
                unit = get_value(value_and_unit, ["unit"])

            return convert_to_native_unit(unit, native_unit, value)

    raise Exception(f"Failed to convert a value to {native_unit}: no valid value/unit pair is provided from the JSON object.")

def get_value_in_unit(native_unit:str, dictionary:dict, keys:list, fail_value=0) -> float | str | bool:
    """Takes a sequence of JSON keys from the given dictionary where
    a JSON object with a value/unit pair must be located.
    Returns the value of this JSON object in the requested `native_unit`.
    """
    if json_exists_and_not_null(dictionary, keys):
        value_unit_pair = json_extract(dictionary, keys)
        if (not object_value_is_null(value_unit_pair) or (native_unit=="string") or (native_unit=="bool")):
            value = json_convert_to_native_unit(native_unit, value_unit_pair)

            if value is not None:
                return value

    return fail_value

def add_filters_to_list(filter_list:list, json_object:dict, keys:list) -> list:
    """Add filters from a given key sequence in the json object to the given filter list."""
    if json_exists_and_not_null(json_object, keys):
        filters = json_extract(json_object, keys)
        if isinstance(filters, list):
            for f in filters:
                new_filter = Filter()
                new_filter.set_from_json(f)
                filter_list.append(new_filter)
        elif isinstance(filters, dict):
            # If no array is given, maybe just
            # one filter is defined as an object.
            new_filter = Filter()
            new_filter.set_from_json(f)
            filter_list.append(new_filter)

    return filter_list

# -----------------------
# Further little helpers
# -----------------------

def list_mean(l:list) -> float:
    """Mean value for a list of values.

    Parameters
    ----------
    l : list
        List of values, e.g. 'float'.

    Returns
    -------
    mean : float
        The mean of the value list.
    """
    return sum(l) / len(l)

def list_mean_and_stddev(l:list) -> tuple[float, float]:
    """Mean and standard deviation of a list of values.

    Parameters
    ----------
    l : list
        List of values, e.g. 'float'.

    Returns
    -------
    mean : float
        The mean of the value list.

    standard_deviation : float
        Root mean square deviation (RMSD) of the value list.
    """
    msqDev = 0
    mean = list_mean(l)
    for v in l:
        msqDev += math.pow(v - mean, 2)

    msqDev /= len(l)

    return mean, math.sqrt(msqDev)

def gaussian(x:float, mu:float, sigma:float, A:float) -> float:
    """Gaussian function.

    Parameters
    ----------
    x : float
        Input value (x axis).

    mu : float
        Distribution's mean value.

    sigma : float
        Distribution's standard deviation.

    A : float
        Distribution's amplitude, i.e., the maximum value.

    Returns
    -------
    y : float
        Output value (y axis).
        `y = A * exp((x-mu)**2 / (2*sigma**2))`
    """
    return A*numpy.exp(-(x-mu)*(x-mu)/(2.0*sigma*sigma))

def poly4(x, a:float, b:float, c:float, d:float, e:float) -> float:
    """ Fourth order polynomial, used for smoothing.

    Parameters
    ----------
    x : float
        Input value (x axis).

    a : float

    b : float

    c : float

    d : float

    e : float

    Returns
    -------
    y : float
        Output value (y axis).
        `y = ax^4 + bx^3 + cx^2 + dx + e`
    """
    return a*(x**4) + b*(x**3) + c*(x**2) + d*x + e

def divide_and_error(muA:float, muB:float, errA:float, errB:float) -> tuple[float, float]:
    """ Error propagation upon division; estimation of largest error. """
    value = muA / muB
    err = errA/abs(muB) + errB*abs(muA/(muB**2))
    return value, err

def divide_and_gaussian_error(muA:float, muB:float, sigmaA:float, sigmaB:float) -> tuple[float, float]:
    """ Gaussian error propagation upon division. """
    value = muA / muB
    uncertainty = math.sqrt((sigmaA**2)/(muB**2) + (sigmaB**2)*(muA**2)/(muB**4))
    return value, uncertainty

def ratios(values:list) -> list:
    """ Calculate ratio to preceding value, needed for step wedge evaluations. """
    results = []
    for v in range(1, len(values)):
        results.append(values[v] / values[v-1])

    return results


""" Unit conversions for values from CTSimU scenario descriptions (JSON files). """

def in_mm_json(jsonVal:dict) -> float:
    """ Convert JSON value/unit pair to mm. """
    if ("value" in jsonVal) and ("unit" in jsonVal):
        value = jsonVal["value"]
        unit  = jsonVal["unit"]

        if(unit == "mm"):
            return value
        elif(unit == "nm"):
            return (value * 1e-6)
        elif(unit == "um"):
            return (value * 1e-3)
        elif(unit == "cm"):
            return (value * 10)
        elif(unit == "dm"):
            return (value * 100)
        elif(unit == "m"):
            return (value * 1000)

        raise Exception(unit + " is not a valid unit of length.")
    else:
        raise KeyError("\"value\" or \"unit\" missing.")

def in_rad_json(jsonVal:dict) -> float:
    """ Convert JSON value/unit pair to radians. """
    if ("value" in jsonVal) and ("unit" in jsonVal):
        value = jsonVal["value"]
        unit  = jsonVal["unit"]

        if(unit == "rad"):
            return value
        elif(unit == "deg"):
            return ((value * math.pi) / 180.0)

        raise Exception(unit + " is not a valid angular unit.")
    else:
        raise KeyError("\"value\" or \"unit\" missing.")