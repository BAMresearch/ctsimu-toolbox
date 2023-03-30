# -*- coding: UTF-8 -*-
"""
Helper functions and variables, especially for handling CTSimU JSON objects as well as unit conversions and other types of conversions.
"""

import json

def read_json_file(filename):
    return json.load(filename)

def value_is_null_or_zero(value):
    """ Checks if a specific value is set to `null` or 0. """
    if value is None:
        return True

    if value == 0:
        return True

    return False

def object_value_is_null(json_obj):
    """ Checks if a JSON object has a `"value"` parameter  and if this parameter is set to `null`. """
    if "value" in json_obj:
        if json_obj["value"] is None:
            return True

        return False

    return True

def object_value_is_null_or_zero(json_obj):
    """ Checks if a JSON object has a `"value"` parameter and if this parameter is set to `null` or 0. """
    if object_value_is_null(json_obj):
        return True

    return value_is_null_or_zero(json_obj["value"])

def get_value(dictionary, keys, fail_value=None):
    currentElement = dictionary

    if keys is None:
        return fail_value
    else:
        for key in keys:
            if currentElement is not None:
                if key in currentElement:
                    currentElement = currentElement[key]
                else:
                    return fail_value
            else:
                return fail_value

    return currentElement

def json_exists(dictionary, keys):
    """Checks if the key sequency exists."""
    currentElement = dictionary

    if keys is None:
        return False
    else:
        for key in keys:
            if currentElement is not None:
                if key in currentElement:
                    currentElement = currentElement[key]
                else:
                    return False
            else:
                return False

    return True

def json_exists_and_not_null(dictionary, keys):
    """Checks if the key sequence exists and its final value is not `null`."""
    if get_value(dictionary, keys) is not None:
        return True

    return False

def json_extract(dictionary, keys, fail_value=None):
    """Get the JSON sub-object that is located by a given sequence of `keys` in the JSON dictionary."""

    # For the Python toolbox, this is the same as get_value.
    # Function remains for comparibility with the Tcl code of the aRTist module.
    return get_value(dictionary, keys, fail_value)

def json_extract_from_possible_keys(dictionary, key_sequences, fail_value=None):
    """Searches the JSON object for each key sequence in the given list of key_sequences. The first sequence that exists will return an extracted JSON object."""
    for keys in key_sequences:
        if json_exists(dictionary, keys):
            return json_extract(dictionary, keys)

    return fail_value

"""
Unit Conversion
-----------------------------
Unit conversion functions take a JSON object that must
contain a `"value"` and a `"unit"`. Each function supports
the allowed units from the CTSimU file format specification.
"""

