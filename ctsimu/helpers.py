# -*- coding: UTF-8 -*-

import os
import json
import math
import csv
import numpy
from scipy import optimize, fft


ctsimu_supported_scenario_version = {
    "major": 1,
    "minor": 2
}

ctsimu_supported_metadata_version = {
    "major": 1,
    "minor": 2
}

ctsimu_valid_axes   = ["z", "y", "x", "w", "v", "u", "t", "s", "r"]
ctsimu_axis_strings = ["r", "s", "t", "u", "v", "w", "x", "y", "z"]
ctsimu_world_axis_designations  = ["x", "y", "z"]
ctsimu_local_axis_designations  = ["u", "v", "w"]
ctsimu_sample_axis_designations = ["r", "s", "t"]

openct_converter = {
    "datatype": {
        "uint8": "UInt8",
        "uint16": "UInt16",
        "uint32": "UInt32",
        "int8": "Int8",
        "int16": "Int16",
        "int32": "Int32",
        "float32": "Float32"
    },
    "endian": {
        "little": "Little",
        "big": "Big"
    }
}

cera_converter = {
    "datatype": {
        "uint16": "uint16",
        "float32": "float"
    }
}

def is_version_supported(supported_version:dict, version_to_test:dict) -> bool:
    """Test if the given version is supported by the toolbox.

    Parameters
    ----------
    supported_version : dict
        The version against which should be tested: a dictionary with the keys
        `"major"` and `"minor"`. Pass the global variable
        `ctsimu_supported_scenario_version` or `ctsimu_supported_metadata_version`.

    version_to_test : dict
        The version to test if it is supported. Must be a dictionary with
        the keys `"major"` and `"minor"`.

    Returns
    -------
    is_supported : bool
        If the given version number is supported by the toolbox scenario module.
    """

    if supported_version["major"] > version_to_test["major"]:
        return True

    if supported_version["major"] == version_to_test["major"]:
        if supported_version["minor"] >= version_to_test["minor"]:
            return True

    return False

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
        None, "mm", "rad", "deg", "s", "mA", "kV", "g/cm^3", "lp/mm", "deg/s", "bool", "string"
    ]
valid_dummy_units = ["px", "1/J", "C", "F", "K", "relative"] # units that are not converted
native_units_to_omit_in_json_file = [None, "bool", "string"]

def is_valid_native_unit(native_unit:str) -> bool:
    """Check if given string is a valid native unit.

    Parameters
    ----------
    native_unit : str
        The string to check if it is a valid native unit.

    Returns
    -------
    is_valid : bool
        `True` if the provided string is a valid native unit, `False` if not.
    """

    if (native_unit in valid_native_units) or (native_unit in valid_dummy_units):
        return True

    raise Exception(f"CTSimU: Not a valid native unit: '{native_unit}'. Valid options are: {valid_native_units}.")
    return False

def touch_directory(filename:str):
    """Create folder if it doesn't exist.

    Parameters
    ----------
    filename : str
        Complete path to a file. Can include filename.
    """

    folder = os.path.dirname(filename)

    if folder == "" or folder is None:
            folder = "."

    if not os.path.exists(folder):
        os.makedirs(folder)

    return folder

def join_dir_and_filename(directory:str, filename:str) -> str:
    """Joins directory and filename into a meaningful path.

    Parameters
    ----------
    directory: str
        Directory part of the path. Set to `None` for no
        directory.

    filename : str
        Filename part of the path.

    Returns
    -------
    full_path : str
        Fully joined path.
    """

    if directory is not None:
        if directory != "":
            return os.path.join(directory, filename)

    return filename

def abspath_of_referenced_file(filepath:str, referenced:str) -> str:
    """Get the absolute path of an external file `referenced` in `file`.

    Parameters
    ----------
    filepath : str
        The path of a file (e.g. scenario or metadata) in which
        an external file is referenced.

    referenced : str
        The path (relative to `filepath` or absolute)
        of the referenced file.

    Returns
    -------
    referenced_abs : str
        Absolute path of the referenced file.
    """
    if not os.path.isabs(referenced):
        if filepath is not None:
            file_abspath = os.path.abspath(filepath)
            file_absdir  = os.path.dirname(file_abspath)
            referenced_abspath = join_dir_and_filename(file_absdir, referenced)

            # Simplify relative dots away:
            referenced_abspath = os.path.abspath(referenced_abspath)
            return referenced_abspath
        else:
            # Assume that an abspath from the current working
            # directory is requested:
            return os.path.abspath(referenced)
    else:
        return referenced  # already an absolute path

def read_json_file(filename:str) -> dict:
    """Read a JSON file into a Python dictionary.

    Parameters
    ----------
    filename : str
        Filename of the JSON file.

    Returns
    -------
    dictionary : dict
        Dictionary representation of the JSON structure.
    """

    if os.path.isfile(filename):
        if os.path.exists(filename):
            with open(filename, 'r', encoding='utf-8') as f:
                json_dict = json.load(f)
                f.close()
                return json_dict
        else:
            raise Exception(f"File not found: '{filename}'")
    else:
        raise Exception(f"File not found: '{filename}'")

    raise Exception(f"Cannot read JSON file: '{filename}'")

def counter_format(n:int, zero_padding:bool=True) -> str:
    """Create a default counter format for sequentially numbered files.

    Parameters
    ----------
    n : int
        Number of files, starting from zero.

    zero_padding : bool
        Are the files zero-padded from the left?

    Returns
    -------
    counter_format : str
        Counter format string, e.g. '%05d' for `n`=20000 files.
        At least four digits will be assumed, even if `n` is smaller.
    """
    digits = 4

    # For anything bigger than 10000 projections (0000 ... 9999)
    # we need more filename digits:
    if n > 10000:
        digits = int(math.ceil(math.log10(float(n))))

    if zero_padding:
        pcformat = f"%0{int(digits)}d"
    else:
        pcformat = f"%d"

    return pcformat


def convert(converter_dict:dict, key:str) -> str:
    """Map a string to a different (converted) string.

    Parameters
    ----------
    converter_dict : dict
        Mapping dictionary.

    key : str
        Key to be converted.

    Returns
    -------
    mapped_value : str
        String that's mapped to the key.
    """
    if key is None:
        return None

    if key in converter_dict:
        return converter_dict[key]
    else:
        raise Exception(f"Cannot convert: key '{key}' not found in converter dictionary.")

def write_json_file(filename:str, dictionary:dict):
    """Write a JSON file from a given Python dictionary.

    Parameters
    ----------
    filename : str
        Filename of the JSON file.

    dictionary : dict
        Dictionary for the JSON file.
    """

    folder = touch_directory(filename)
    if os.path.exists(folder):
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(dictionary, f, ensure_ascii=False, indent="\t")
            f.close()
    else:
        raise Exception(f"Error writing JSON file. Directory does not exist: {folder}")

def read_csv_file(filename:str) -> dict:
    """Read a CSV file.

    Parameters
    ----------
    filename : str
        Filename of the CSV file to read.

    Returns
    -------
    data : list
        Each item in the returned list contains the values of one column from the CSV file.
    """

    values = []
    n_columns = 0

    with open(filename, newline='') as f:
        try:
            # Detect the CSV dialect: comma or tab-separated?
            dialect = csv.Sniffer().sniff(f.read(1024))
        except:
            dialect = None

        # Return to beginning
        f.seek(0)

        reader = csv.reader(f, dialect)
        for row in reader:
            if row[0].startswith('#'):
                # ignore commented lines
                continue

            if n_columns == 0:
                # Number of columns apparently not initialized yet.
                n_columns = len(row)
                for i in range(n_columns):
                    # Append an empty list for each column:
                    values.append([])

            for col, entry in enumerate(row):
                values[col].append(entry)

    return values

def value_is_null(value) -> bool:
    """Check if a specific JSON value is set to `null`.

    Parameters
    ----------
    value
        Value to check for nullness.

    Returns
    -------
    is_null : bool
        `True` if the value corresponds to a JSON `null`,
        `False` if it is something else.
    """
    if value is None:
        return True

    return False

def value_is_null_or_zero(value) -> bool:
    """Check if a specific value is set to `null` or `0`.

    Parameters
    ----------
    value
        Value to check for nullness or zeroness.

    Returns
    -------
    is_null_or_zero : bool
        `True` if the value corresponds to a JSON `null`
        or has the numerical value `0`, `False` if it is something else.
    """
    if value is not None:
        if value != 0:
            return False

    return True

def object_value_is_null(json_obj:dict) -> bool:
    """Check if a CTSimU JSON parameter object represents a
    `null` value, either because its `"value"` property is
    set to a JSON `null` or because it does not have
    a `"value"` property at all.

    Parameters
    ----------
    json_obj : dict
        Dictionary (as from a JSON structure) to check.

    Returns
    -------
    is_null : bool
        `True` if the object's `"value"` corresponds to a JSON `null`
        or if the object does not define a `"value"`.
        `False` if the object value is something else.

    Raises
    ------
    TypeError
        If `json_obj` is not a dictionary.
    """
    if json_obj is None:
        return True

    if not isinstance(json_obj, dict):
        return False

    if "value" in json_obj:
        return value_is_null(json_obj["value"])

    return True

def object_value_is_null_or_zero(json_obj:dict) -> bool:
    """Check if a CTSimU JSON parameter object represents
    `null` value or the numerical value `0`, either because
    its `"value"` property is set to a JSON `null` or `0`,
    or because it does not have a `"value"` property at all.

    Parameters
    ----------
    json_obj : dict
        Dictionary (as from a JSON structure) to check.

    Returns
    -------
    is_null_or_zero : bool
        `True` if the object's `"value"` corresponds to a JSON `null`
        or if the object does not define a `"value"`.
        `False` if the object value is something else.

    Raises
    ------
    TypeError
        If `json_obj` is not a dictionary.
    """
    if not isinstance(json_obj, dict):
        raise TypeError("object_value_is_null_or_zero() expects dict for the `json_obj`.")

    if "value" in json_obj:
        return value_is_null_or_zero(json_obj["value"])

    return True

def get_value_or_none(dictionary:dict, *keys:str) -> float | str | bool | dict | list:
    """Get the dictionary value that is located at the given
    sequence of `keys`. If it cannot be found, return `None`.

    Parameters
    ----------
    dictionary : dict
        Dictionary (as from a JSON structure).

    *keys : str
        Sequence of keys that identify a location in the dictionary tree.

    Returns
    -------
    value : float or str or bool or dict or list
        The entry in the `dictionary`, located by the sequence of `keys`,
        or `None` if the entry cannot be found.

    Raises
    ------
    TypeError
        If `dictionary` is not a Python `dict`.
    """
    if not isinstance(dictionary, dict):
        raise TypeError("get_value_or_none() expects dict as first argument for the `dictionary`.")

    if len(keys) == 0:
        return dictionary

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

def get_value(dictionary:dict, keys:list, fail_value=None) -> float | str | bool | dict | list:
    """Get the specific value of the parameter that is
    located at the given sequence of `keys` in the JSON dictionary.

    If that element cannot be found or is set to `null`,
    return the `fail_value`.

    Parameters
    ----------
    dictionary : dict
        Dictionary (as from a JSON structure).

    keys : list
        List of keys that identify a location in the dictionary tree.

    fail_value : float or str or bool or dict or list, optional
        Value to return if the given list of keys cannot be found
        in the dictionary.

    Returns
    -------
    value : float or str or bool or dict or list
        The entry in the `dictionary`, located by the sequence of `keys`
        or the `fail_value` if the entry cannot be found.

    Raises
    ------
    TypeError
        If `dictionary` is not a Python `dict`.
    """
    if not isinstance(dictionary, dict):
        raise TypeError("get_value() expects dict as first argument for the `dictionary`.")

    if len(keys) == 0:
        return dictionary

    result = get_value_or_none(dictionary, *keys)
    if result is None:
        return fail_value

    return result

def json_exists(dictionary:dict, keys:list) -> bool:
    """Check if the given key sequence can be found in the dictionary tree.

    Parameters
    ----------
    dictionary : dict
        Dictionary (as from a JSON structure).

    keys : list
        List of keys that identify a location in the dictionary tree.

    Returns
    -------
    exists : bool
        `True` if the given key sequence identifies an object or
        value in the dictionary, `False` if nothing can be found
        at the given sequence of keys.

    Raises
    ------
    TypeError
        If `dictionary` is not a Python `dict`.
    """

    if not isinstance(dictionary, dict):
        return False

    if len(keys) == 0:
        return True

    current_element = dictionary
    for key in keys:
        if key in current_element:
            current_element = current_element[key]
        else:
            return False

    return True

def json_isnull(dictionary:dict, keys:list) -> bool:
    """Check if the value at the given key sequence corresponds
    to a JSON `null`, or if the value cannot be found at all.

    Parameters
    ----------
    dictionary : dict
        Dictionary (as from a JSON structure).

    keys : list
        List of keys that identify a location in the dictionary tree.

    Returns
    -------
    is_null : bool
        `True` if the given key sequence identifies a `null` value
        or if the value cannot be found at all. `False` otherwise.

    Raises
    ------
    TypeError
        If `dictionary` is not a Python `dict`.
    """
    if not isinstance(dictionary, dict):
        raise TypeError("json_isnull() expects dict as first argument for the `dictionary`.")
    if len(keys) == 0:
        return True

    v = get_value(dictionary, keys)
    if v is None:
        return True

    return False

def json_exists_and_not_null(dictionary:dict, keys:list) -> bool:
    """Check if the value at the given key sequence
    exists and does not correspond to a JSON `null`.

    Parameters
    ----------
    dictionary : dict
        Dictionary (as from a JSON structure).

    keys : list
        List of keys that identify a location in the dictionary tree.

    Returns
    -------
    exists_and_not_null : bool
        `True` if the given key sequence exists and does
        not identify a `null` value, `False` otherwise.

    Raises
    ------
    TypeError
        If `dictionary` is not a Python `dict`.
    """
    if json_exists(dictionary, keys):
        if not json_isnull(dictionary, keys):
            return True

    return False

def json_extract(dictionary:dict, keys:list) -> float | str | bool | dict | list:
    """Get the JSON sub-object that is located
    at a given sequence of `keys` in the JSON dictionary.

    Parameters
    ----------
    dictionary : dict
        Dictionary (as from a JSON structure).

    keys : list
        List of keys that identify a location in the dictionary tree.

    Returns
    -------
    value : float or str or bool or dict or list
        The value, object or array located at the given key sequence,
        `None` if nothing can be found at the given key sequence.

    Raises
    ------
    TypeError
        If `dictionary` is not a Python `dict`.
    """
    return get_value(dictionary, keys)

def json_extract_from_possible_keys(dictionary:dict, key_lists:list):
    """Searches the JSON object for each key sequence in the given
    list of key sequences. The first sequence that exists will be used to
    extract and return the corresponding JSON object.

    Parameters
    ----------
    dictionary : dict
        Dictionary (as from a JSON structure).

    key_lists : list
        List of lists of keys that possibly specify a valid
        location in the dictionary tree.

    Returns
    -------
    value : float or str or bool or dict or list
        The value, object or array located at the first valid key sequence,
        `None` if nothing can be found at any of the given key sequences.

    Raises
    ------
    TypeError
        If `dictionary` is not a Python `dict`.
    """
    for keys in key_lists:
        if json_exists(dictionary, keys):
            return json_extract(dictionary, keys)

    return None

# ------------------------------------
# CTSimU conversions to native units
# ------------------------------------
# Each function supports the allowed units from the CTSimU file format specification.

def in_mm(value:float, unit:str) -> float:
    """Convert a length to mm.

    Parameters
    ----------
    value : float
        The value to convert.

    unit : str
        The value's current unit.
        Options are: `"nm"`, `"um"`, `"mm"`, `"cm"`, `"dm"`, `"m"`

    Returns
    -------
    value_in_mm : float
        The value converted to millimeters,
        or 'None' if the original value was `None`.

    Raises
    ------
    ValueError
        If the given `unit` is not a valid unit of length.
    """
    if value is not None:
        if unit == "nm": return (value * 1.0e-6)
        if unit == "um": return (value * 1.0e-3)
        if unit == "mm": return value
        if unit == "cm": return (value * 10.0)
        if unit == "dm": return (value * 100.0)
        if unit == "m":  return (value * 1000.0)
    else:
        return None

    raise ValueError(f"Not a valid unit of length: '{unit}'.")

def in_rad(value:float, unit:str="deg") -> float:
    """Convert an angle to radians.

    Parameters
    ----------
    value : float
        The value to convert.

    unit : str
        The value's current unit.
        Options are: `"rad"`, `"deg"`

    Returns
    -------
    value_in_rad : float
        The value converted to rad,
        or 'None' if the original value was `None`.

    Raises
    ------
    ValueError
        If the given `unit` is not a valid angular unit.
    """
    if value is not None:
        if unit == "deg": return (value * math.pi / 180.0)
        if unit == "rad": return value
    else:
        return None

    raise ValueError(f"Not a valid angular unit: '{unit}'.")

def in_deg(value:float, unit:str="rad") -> float:
    """Convert an angle to degrees.

    Parameters
    ----------
    value : float
        The value to convert.

    unit : str
        The value's current unit.
        Options are: `"rad"`, `"deg"`

    Returns
    -------
    value_in_deg : float
        The value converted to degrees,
        or 'None' if the original value was `None`.

    Raises
    ------
    ValueError
        If the given `unit` is not a valid angular unit.
    """
    if value is not None:
        if unit == "rad": return (value * 180.0 / math.pi)
        if unit == "deg": return value
    else:
        return None

    raise ValueError(f"Not a valid angular unit: '{unit}'.")

def in_s(value:float, unit:str) -> float:
    """Convert a time to seconds.

    Parameters
    ----------
    value : float
        The value to convert.

    unit : str
        The value's current unit.
        Options are: `"ms"`, `"s"`, `"min"`, `"h"`

    Returns
    -------
    value_in_s : float
        The value converted to seconds,
        or 'None' if the original value was `None`.

    Raises
    ------
    ValueError
        If the given `unit` is not a valid unit of time.
    """
    if value is not None:
        if unit == "ms":  return (value * 1.0e-3)
        if unit == "s":   return value
        if unit == "min": return (value * 60.0)
        if unit == "h":   return (value * 3600.0)
    else:
        return None

    raise ValueError(f"Not a valid unit of time: '{unit}'.")

def in_mA(value:float, unit:str) -> float:
    """Convert an electric current to mA.

    Parameters
    ----------
    value : float
        The value to convert.

    unit : str
        The value's current unit.
        Options are: `"uA"`, `"mA"`, `"A"`

    Returns
    -------
    value_in_mA : float
        The value converted to milliamps,
        or 'None' if the original value was `None`.

    Raises
    ------
    ValueError
        If the given `unit` is not a valid unit of electric current.
    """
    if value is not None:
        if unit == "uA": return (value * 1.0e-3)
        if unit == "mA": return value
        if unit == "A":  return (value * 1000.0)
    else:
        return None

    raise ValueError(f"Not a valid unit of electric current: '{unit}'.")

def in_kV(value:float, unit:str) -> float:
    """Convert a voltage to kV.

    Parameters
    ----------
    value : float
        The value to convert.

    unit : str
        The value's current unit.
        Options are: `"V"`, `"kV"`, `"MV"`

    Returns
    -------
    value_in_kV : float
        The value converted to kilovolts,
        or 'None' if the original value was `None`.

    Raises
    ------
    ValueError
        If the given `unit` is not a valid unit of electric voltage.
    """
    if value is not None:
        if unit == "V":  return (value * 1.0e-3)
        if unit == "kV": return value
        if unit == "MV": return (value * 1000.0)
    else:
        return None

    raise ValueError(f"Not a valid unit of electric voltage: '{unit}'.")

def in_deg_per_s(value:float, unit:str) -> float:
    """Convert an angular velocity to deg/s.

    Parameters
    ----------
    value : float
        The value to convert.

    unit : str
        The value's current unit.
        Options are: `"rad/s"`, `"rad/min"`, `"rad/h"`, `"deg/s"`,
        `"deg/min"`, `"deg/h"`

    Returns
    -------
    value_in_deg_per_s : float
        The value converted to degrees per second,
        or 'None' if the original value was `None`.

    Raises
    ------
    ValueError
        If the given `unit` is not a valid unit of angular velocity.
    """
    if value is not None:
        if unit == "rad/s":   return math.degrees(value)
        if unit == "rad/min": return math.degrees(value / 60.0)
        if unit == "rad/h":   return math.degrees(value / 3600.0)
        if unit == "deg/s":   return value
        if unit == "deg/min": return value / 60.0
        if unit == "deg/h":   return value / 3600.0
    else:
        return None

    raise ValueError(f"Not a valid unit of angular velocity: '{unit}'.")

def in_g_per_cm3(value:float, unit:str) -> float:
    """Convert a mass density to g/cm³.

    Parameters
    ----------
    value : float
        The value to convert.

    unit : str
        The value's current unit.
        Options are: `"kg/m^3"`, `"g/cm^3"`

    Returns
    -------
    value_in_g_per_cm3 : float
        The value converted to grams per cubic centimeter,
        or 'None' if the original value was `None`.

    Raises
    ------
    ValueError
        If the given `unit` is not a valid unit of mass density.
    """
    if value is not None:
        if unit == "kg/m^3": return (value * 1.0e-3)
        if unit == "g/cm^3": return value
    else:
        return None

    raise ValueError(f"Not a valid unit of mass density: '{unit}'.")

def in_lp_per_mm(value:float, unit:str) -> float:
    """Convert a resolution to line pairs per millimeter.

    Parameters
    ----------
    value : float
        The value to convert.

    unit : str
        The value's current unit.
        Options are: `"lp/mm"`, `"lp/cm"`, `"lp/dm"`, `"lp/m"`

    Returns
    -------
    value_in_lp_per_mm : float
        The value converted to line pairs per millimeter,
        or 'None' if the original value was `None`.

    Raises
    ------
    ValueError
        If the given `unit` is not a valid unit for resolution.
    """
    if value is not None:
        if unit == "lp/mm": return value
        if unit == "lp/cm": return (value * 0.1)
        if unit == "lp/dm": return (value * 0.01)
        if unit == "lp/m":  return (value * 0.001)
    else:
        return None

    raise ValueError(f"Not a valid unit for resolution: '{unit}'.")

def from_bool(value) -> bool:
    """Convert `value` into a true Python boolean.
    (For example, `1` becomes `True`, and `0` becomes `False`.)

    Parameters
    ----------
    value
        Any value that can be checked with a Python `if`.

    Returns
    -------
    b : bool
        `True` if the value evaluates to `True`,
        otherwise `False`.
    """
    if value:
        return True

    return False

def convert_SNR_FWHM(SNR_or_FWHM:float, mean:float) -> float:
    """Converts between SNR and Gaussian FWHM for a
    given `mean` value of the Gaussian distribution (e.g., intensity).

    Parameters
    ----------
    SNR_or_FWHM : float
        SNR value (signal to noise ratio) or
        FWHM value (full width at half maximum).

    mean : float
        Mean value of the Gaussian distribution.

    Returns
    -------
    FWHM_or_SNR : float
        FWHM if an SNR was given, or the SNR if an FWHM was given.
    """

    return 2.0 * math.sqrt(2.0 * log(2.0)) * float(mean) / float(SNR_or_FWHM)

def convert_to_native_unit(given_unit:str, native_unit:str, value:float) -> float | str | bool:
    """Check which native unit is requested, convert value accordingly.

    If the `native_unit` is set to None, the `value` will simply
    be returned as it is (i.e., unaltered).

    Parameters
    ----------
    given_unit : str
        The unit of the given value.

    native_unit : str
        The native unit into which the value must be converted.
        Possible values: `None`, `"mm"`, `"rad"`, `"deg"`, `"s"`, `"mA"`,
        `"kV"`, `"g/cm^3"`, `"lp/mm"`, `"bool"`, `"string"`.

    value : float or str or bool
        Value to convert.

    Returns
    -------
    value_in_native_unit : float or str or bool
        Converted value.

    Raises
    ------
    ValueError
        If the `given_unit` is not compatible with the `native_unit`.
        For example, the function cannot convert a length into a time.
    """
    if native_unit is None:
        return value
    elif native_unit in valid_dummy_units:
        return value
    else:
        if native_unit == "string": return value
        if native_unit == "mm":     return in_mm(value, given_unit)
        if native_unit == "s":      return in_s(value, given_unit)
        if native_unit == "deg":    return in_deg(value, given_unit)
        if native_unit == "rad":    return in_rad(value, given_unit)
        if native_unit == "mA":     return in_mA(value, given_unit)
        if native_unit == "kV":     return in_kV(value, given_unit)
        if native_unit == "deg/s":  return in_deg_per_s(value, given_unit)
        if native_unit == "g/cm^3": return in_g_per_cm3(value, given_unit)
        if native_unit == "lp/mm":  return in_lp_per_mm(value, given_unit)
        if native_unit == "bool":   return from_bool(value)

    raise ValueError(f"Native unit '{native_unit}' is incompatible with the given unit '{given_unit}'.")
    return None

def json_convert_to_native_unit(native_unit:str, value_and_unit:dict, fallback_json_unit:str=None) -> float | str | bool:
    """Convert a value/unit pair from a dictionary into
    the requested `native_unit`.

    Works like the function `convert_to_native_unit`, but takes
    a JSON object `value_and_unit`, i.e. a dictionary that must contain
    a `value` and an associated `unit`. Checks which native unit is requested,
    and converts the JSON `value` accordingly.

    `fallback_json_unit` is used if the unit is not specified
    in the `value_and_unit` JSON object.

    Parameters
    ----------
    native_unit : str
        The native unit into which the value will be converted.
        Possible values: `None`, `"mm"`, `"rad"`, `"deg"`, `"s"`, `"mA"`,
        `"kV"`, `"g/cm^3"`, `"lp/mm"`, `"bool"`, `"string"`.

    value_and_unit : dict
        Dictionary that specifies a `"value"` and a `"unit"`,
        as imported from a CTSimU JSON parameter object.

    fallback_json_unit : str, optional
        The fallback unit to be used if the `value_and_unit` dictionary
        does not specify a unit.

    Returns
    -------
    value_in_native_unit : float or str or bool
        Converted value.

    Raises
    ------
    ValueError
        If the unit given in the dictionary is not compatible
        with the requested `native_unit`. For example, the function
        cannot convert a length into a time.
        Also, if no valid value/unit pair is provided in the given
        `value_and_unit` dictionary.
    """

    if native_unit is None:
        # No native unit given. Simply return the value.
        return get_value(value_and_unit, ["value"])
    elif native_unit == "bool":
        # This is not a value/unit dictionary, but a boolean.
        return from_bool(value_and_unit)
    elif native_unit == "string":
        if json_exists(value_and_unit, ["value"]):
            return str(get_value(value_and_unit, ["value"]))
        else:
            if isinstance(value_and_unit, str):
                return value_and_unit

        raise TypeError(f"Given value does not seem to be a string: {value_and_unit}")
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

    raise ValueError(f"Failed to convert a value to '{native_unit}': no valid value/unit pair is provided from the JSON object.")

def get_value_in_native_unit(native_unit:str, dictionary:dict, keys:list, fail_value=None) -> float | str | bool:
    """Get a parameter value in the `native_unit`,
    located at a sequence of `keys` in a given `dictionary`.

    Parameters
    ----------
    native_unit : str
        The native unit into which the value will be converted.
        Possible values: `None`, `"mm"`, `"rad"`, `"deg"`, `"s"`, `"mA"`,
        `"kV"`, `"g/cm^3"`, `"lp/mm"`, `"bool"`, `"string"`.

    dictionary : dict
        Dictionary representation of the JSON structure.

    keys : list
        List of strings that identify the location of the paramter in
        the dictionary.

    fail_value : float or str or bool
        The value that will be returned if nothing can be found
        at the given sequence of `keys`.

    Returns
    -------
    value_in_native_unit : float or str or bool
        Converted value or `fail_value` if nothing could be found.

    Raises
    ------
    ValueError
        If the unit given in the dictionary is not compatible
        with the requested `native_unit`. For example, the function
        cannot convert a length into a time.
        Also, if no valid value/unit pair is provided in the given
        `value_and_unit` dictionary.
    """
    if json_exists_and_not_null(dictionary, keys):
        value_unit_pair = json_extract(dictionary, keys)
        if (not object_value_is_null(value_unit_pair) or (native_unit=="string") or (native_unit=="bool")):
            value = json_convert_to_native_unit(native_unit, value_unit_pair)

            if value is not None:
                return value

    return fail_value

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