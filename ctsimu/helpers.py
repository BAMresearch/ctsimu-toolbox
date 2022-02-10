# -*- coding: UTF-8 -*-

import math
import numpy
from scipy import optimize, fft

def log(message):
    """Print an output message."""
    print(message)

def getFieldOrNone(dictionary, *fields):
    currentElement = dictionary

    if fields is None:
        return None
    else:
        for field in fields:
            if currentElement != None:
                if field in currentElement:
                    currentElement = currentElement[field]
                else:
                    return None
            else:
                return None

    return currentElement

def rad2deg(rad):
    return rad * 180.0 / math.pi

def deg2rad(deg):
    return deg * math.pi / 180.0

def listMean(l):
    return sum(l) / len(l)

def listMeanAndStdDev(l):
    msqDev = 0
    mean = listMean(l)
    for v in l:
        msqDev += math.pow(v - mean, 2)

    msqDev /= len(l)

    return mean, math.sqrt(msqDev)

# Convert to unit, assuming angle is in deg:
def deg2x(angle, unit='deg'): 
    if unit == 'deg':
        return angle
    elif unit == 'rad':
        if angle != None:
            return ballTools.deg2rad(angle)
        else:
            return None
    else:
        raise Exception("Angle unit must be 'deg' or 'rad'.")

# Convert to unit, assuming angle is in rad:
def rad2x(angle, unit='rad'):
    if unit == 'rad':
        return angle
    elif unit == 'deg':
        if angle != None:
            return ballTools.rad2deg(angle)
        else:
            return None
    else:
        raise Exception("Angle unit must be 'deg' or 'rad'.")

def gaussian(x, mu, sigma, A):
    return A*numpy.exp(-(x-mu)*(x-mu)/(2.0*sigma*sigma))

def poly4(x, a, b, c, d, e):
    """ Fourth order polynomial, used for smoothing during MTF calculation
        and for smoothing of the Monte-Carlo grey value profiles. """
    return a*(x**4) + b*(x**3) + c*(x**2) + d*x + e

def divideAndError(muA, muB, errA, errB):
    """ Error propagation upon division; estimation of largest error. """
    value = muA / muB
    err = errA/abs(muB) + errB*abs(muA/(muB**2))
    return value, err

def divideAndGaussianError(muA, muB, sigmaA, sigmaB):
    """ Gaussian error propagation upon division. """
    value = muA / muB
    uncertainty = math.sqrt((sigmaA**2)/(muB**2) + (sigmaB**2)*(muA**2)/(muB**4))
    return value, uncertainty

def ratios(values):
    """ Calculate ratio to preceding value, needed for step wedge evaluations. """
    results = []
    for v in range(1, len(values)):
        results.append(values[v] / values[v-1])

    return results


""" Unit conversions for values from CTSimU scenario descriptions (JSON files). """

def in_mm(jsonVal):
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

def in_rad(jsonVal):
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