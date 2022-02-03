# -*- coding: UTF-8 -*-
"""Calculate the MTF from a given ESF or LSF."""

import math
import numpy
from scipy import optimize, fft

def MTF(positions, ESF=None, LSF=None, ESFsmoothingWidth=8):
    """ Calculate the MTF for given edge spread function (ESF) or line spread function (LSF) and positions. """

    if not(ESF is None):
        # Smooth the ESF using a Gaussian-weighted 4th order polynomial fit:
        nSamples = len(ESF)
        smoothedESF = numpy.zeros(nSamples, dtype=numpy.float64)

        w = ESFsmoothingWidth  # number of samples for an interpolation fit, in each direction

        for pos in range(nSamples):
            # Fit to weighted sub-segment:
            startIdx   = int(max(0, pos-w))
            stopIdx    = int(min(pos+w+1, nSamples))
            centralIdx = int(min(pos, startIdx+w))

            fitPositions = positions[startIdx:stopIdx]
            fitESF       = ESF[startIdx:stopIdx]

            minVal = numpy.amin(fitESF)
            maxVal = numpy.amax(fitESF)
            delta  = maxVal-minVal

            smoothedValue = 0

            #if delta > 0.1:  # skip if all elements are (roughly) equal
            sigma = numpy.zeros(stopIdx-startIdx, dtype=numpy.float64)
            for i in range(stopIdx-startIdx):
                g = float(i - (centralIdx - startIdx))
                sigma[i] = math.exp(math.pow(2.0*g/w, 2.0))

            try:
                popt, pcov = optimize.curve_fit(f=poly4, xdata=fitPositions, ydata=fitESF, p0=None, sigma=sigma)

                smoothedValue = poly4(x=positions[pos], a=popt[0], b=popt[1], c=popt[2], d=popt[3], e=popt[4])
            except Exception:
                #print("Exception at pos. {}".format(pos-nSamples/2.0))
                smoothedValue = numpy.mean(fitESF)
                
            smoothedESF[pos] = smoothedValue

            print("\rSmoothing ESF... {:.1f}%".format(100.0*pos/nSamples), end='')

        print("\rSmoothing ESF... 100%  ")

        print("Calculating LSF and MTF...")
        # Differentiate smoothed line to get line spread function (LSF)
        diffStepSize = 2.0*(positions[1]-positions[0])

        LSF = numpy.zeros(nSamples, dtype=numpy.float64)
        for pos in range(1, nSamples-1):
            LSF[pos] = (smoothedESF[pos+1] - smoothedESF[pos-1]) / diffStepSize

    elif not(LSF is None):
        nSamples = len(LSF)
        smoothedESF = None
    else:
        raise Exception("MTF: either an ESF or an LSF must be passed.")


    # Subtract background in long tails (10 px at the very left and right of the LSF):
    background = (numpy.mean(LSF[0:10]) + numpy.mean(LSF[-10:])) / 2.0
    if abs(background) > 1e-3:
        LSF -= background
        print("Subtracted LSF background level: {}".format(background))

    # Normalize LSF maximum to 1:
    LSFmax = numpy.amax(LSF)
    if LSFmax != 0:
        LSF /= LSFmax

    # Apply a von Hann window before Fourier transform:
    hann_window_width = nSamples
    hann_window = numpy.hanning(hann_window_width)

    LSF_windowed = LSF * hann_window

    # Normalize smoothed LSF maximum to 1:
    LSFsmoothed_max = numpy.amax(LSF_windowed)
    if LSFsmoothed_max != 0:
        LSF_windowed /= LSFsmoothed_max

    # Pad zeros to nearest full power of 2 (or next-nearest)
    power_of_two = math.floor(math.log(nSamples, 2)) + 1
    if (2**power_of_two - nSamples) < 20:
        # Pad at least 20 zeros or increase power of two:
        power_of_two += 1

    nMTFsamples = 2**power_of_two
    LSF_for_MTF = numpy.zeros(nMTFsamples, dtype=LSF_windowed.dtype)
    LSF_for_MTF[0:len(LSF_windowed)] = LSF_windowed

    print("{} samples before zero-padding, {} samples after.".format(len(LSF_windowed), len(LSF_for_MTF)))

    mtf = numpy.absolute(fft.fft(LSF_for_MTF))

    # Normalize value at zero frequency to 1:
    if mtf[0] != 0:
        mtf = mtf / mtf[0]

    stepsize = 1
    if len(positions) > 1:
        stepsize = positions[1] - positions[0]

    # Calculate the Nyquist frequency and frequency position axis:
    fnyq = 1.0 / (2.0*stepsize)
    mtfpos = numpy.linspace(start=0, stop=2.0*fnyq, num=nMTFsamples, endpoint=False, dtype=numpy.float64)

    return smoothedESF, LSF, LSF_windowed, mtf, mtfpos


def MTFfreq(MTFpos, MTF, modulation=0.2):
    """ Return the frequency for the requested modulation height. """

    mtf_freq = 0
    for f in range(len(MTFpos)):
        if MTF[f] < modulation:
            if f > 0:
                # Linear interpolation:
                x0 = MTFpos[f-1]
                x1 = MTFpos[f]
                y0 = MTF[f-1]
                y1 = MTF[f]

                m = (y1-y0)/(x1-x0)
                n = y0 - m*x0
                mtf_freq = (modulation-n)/m
            else:
                mtf_freq = MTFpos[f]
            
            break

    return mtf_freq
