# -*- coding: UTF-8 -*-
# From Bendix Hartlaub's iSRb tool

from warnings import warn

import numpy as np
#from skimage.measure import profile_line   # toolbox implements its own line profle measurement routine
from scipy.optimize import curve_fit, minimize
from scipy.signal import find_peaks, peak_widths
from time import time

import copy

class OrderError(Exception):
    '''Raise on inappropiate execution order. Correct order is:
    __init__(), profile(), calc_dips(), interpolate()'''
    def __init__(self, message):
        super().__init__(message)
        return

class ResultError(Exception):
    '''Raise on bad calculation results, if unable to continue'''
    def __init__(self, message):
        super().__init__(message)
        return

class Interpolation(object):
    '''Manage Interpolation
    The interpolation process is separated into 4 parts:
    __init__()    - open image
    profile()     - measure profile in image
    calc_dips()   - calculate dips
    interpolate() - calculate isrb from dips
    
    Each part (except _init__) can be re-executed to estimate suitable 
    variables. The Interpolation object saves parameters additional to 
    the parameters needed for further calculations. They simplify
    finding fitting parameters for the calculation.
    '''
    
    def __init__(self, im, pixelsize, SOD=1000, SDD=1000, wire_length = 15,
        wire_spacing = [0.8, 0.63, 0.5, 0.4, 0.32, 0.25, 0.2, 0.16, 0.13, 0.1, 0.08, 0.063, 0.05]):
        '''Create the Interpolation instance with the image representing 
        array and other necessary parameters.
        
        Parameters
        ----------
        im : array-like, shape (n, m)
            The array representing the gray/intensity scale image 
            (background is assumed to have high intensity, bright 
            gray scales).
        pixelsize : scalar
            The pixelsize of the image in milimeters. Only square 
            pixels are possible.
        SOD : float, optional
            Source Object Distance of the scan setup in milimeters.
        SDD : float, optional
            Source Detector Distance of the scan setup in milimeters.
        wire_length : scalar
            Length of the wires in mm
        wire_spacing : array-like, shape (n, )
            Spacing that the wire pairs hold, acending (matches the 
            wire diameter)
            '''
        
        # calculating the scaling of the duplex wires
        if SDD == 0:
            raise ValueError('SDD cannot be 0')
        if SOD == 0:
            raise ValueError('SOD cannot be 0')
        self.scale = SDD/SOD
        if self.scale < 1:
            raise ValueError('SOD cannot be larger than SDD')
        
        # initializing required variables
        self.im = np.asarray(im, dtype = np.float64)
        self.pxsize = float(pixelsize)
        
        # the space inbetween that the wire pairs hold and the wire length
        self.wire_spacing = np.asarray(wire_spacing)
        self.wire_len = wire_length
        
        # variables that will later be set
        self.measure = None # ndarray
        self.dips = None    # ndarray
        self.dip20 = None   # scalar
        self.criticalIndex = None   # Index of last wire pair with dip > 20%.
        
        # variable not necessary to pass between functions but useful to investigate
        self.coords = None  # tuple of profile line properties
        self.peaks = None   # tuples of ndarray, shape (2, x)
        self.max_peaks = None
        self.despike = None
        self.bg = None
        self.dipsoi = None
        self.inter = None

        # Interpolation fit results:
        self.a = 0
        self.b = 0
        self.c = 0

        return
    
    
    
    def quadratic(x, a, b, c):
        '''quadratic function'''
        
        return a*x**2 + b*x + c
    
    def inverted_quadratic(y, a, b, c):
        '''solve quadratic function ax^2+bx+c=y
        
        Returns
        -------
        tuple
            Both possible solutions for x'''
        
        return (-b/(2*a) + np.sqrt((b/a)**2 /4 - (c-y)/a),
                -b/(2*a) - np.sqrt((b/a)**2 /4 - (c-y)/a))
    
    
    
    def _bivariance(self, phi, rho, start, def_kwargs):
        '''Calculate the variance of the vertically calculated variance 
        of the region of interest
        
        Parameters
        ----------
        phi : scalar
            Angle between the profile line and the horizontal plane
        rho : scalar
            Length of the profile line (pixel coordinates)
        start : array-like, shape (2, )
            Pixel coordinate of the profile line start
        def_kwargs : dict
            Options specifying the profile line
        
        Returns
        -------
        bivar : scalar
            variance of vertical variance
        '''
        
        # scipy minimization passes arrays
        phi = phi[0]
        
        def_kwargs['reduce_func'] = np.var
        stop = np.array(start) + rho * np.array((np.cos(phi), np.sin(phi)))
        # matrix indeces differ from cartesian coordinates
        var = profile_line(self.im, start[::-1], stop[::-1], **def_kwargs)
        bivar = np.var(var)
        return bivar
    
    
    
    def profile(self, start, stop, polar_coord = False, optimize = False, rel_width = 0.6):
        '''Measure intensity along a profile line (area) in the loaded image.
        
        Parameters
        ----------
        start : array-like, shape (2, )
            Pixel coordinate of the profile line start.
        stop : array-like, shape (2, )
            Coordinate of the profile line stop. Look at polar_coord for more info.
        polar_coord : bool, optional
            True: 'stop' is treated as polar coordinate (rho, phi), relative to 'start'.
                  Angle in degree, counter-clockwise.
            False: 'stop' is treated as cartesian coordinate (x, y)
        optimize : bool, optional
            If True, optimize the angle (phi) of 'stop' (regardless of coordinate type).
            Optimize by minimizig the variance of vertical variance of the region of 
            interest with the Powell's Method.
        rel_width : scalar, optional
            Relative portion of the wire length used for measuring. values ranging
            from 0.3 to 0.6 are recommended.
        
        Returns
        -------
        None
            Write the intensity profile along the scan line to
            self.measure : ndarray, shape (n, )
            on execution
        
        General Settings
        ----------------
         - linewidth = width * wire_length
         - bi-quadratic filtering for off-pixel coordinates
         - nearest filter for coordinates outside of the image
         - arithmetic mean for aggregation of pixels perpendicular to the line
        
        Notes
        -----
        The possibility of non-square pixels:
        The transformation between pixel-coordinates and cartesian 
        distance-coordinates is not a conform map, if the pixels are not squares. 
        Therefore the measured pixels perpendicular to the profile line in the pixel 
        coordinates would not be perpendicular in distance coordinates.
        
        The measurement is based on skimage.measure.profile_line(). For
        further informations, check out
        https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.profile_line'''
        
        try:
            rel_width = float(rel_width)
        except (TypeError, ValueError):
            raise TypeError("'width' must be float type.")
        if rel_width < 0.3 or rel_width > 0.6:
            warn("Expected 0.3 <= 'width' <= 0.6 but width = {0:.2f}".format(rel_width))
        
        def_kwargs = {'linewidth' : int(np.rint(self.scale * self.wire_len/self.pxsize * rel_width)),
                      'order' : 2,
                      'mode' : 'nearest'}
        
        if optimize:
            print('optimizing')
            t = time()
            if polar_coord:
                rho, phi = stop[0], -np.deg2rad(stop[1])
            else:
                rho = np.linalg.norm(stop-start)
                # geometric scalar product
                phi = np.arccos((stop[0]-start[0])/rho)
            
            # minimization
            sol = minimize(self.bivariance, phi, method='Powell', args=(rho, start, def_kwargs.copy()))
            phi = sol['x']
            print('optimizing done in {:.6f}ms \n'.format((time()-t)*1000))
            print('phi = {:.6f}, bivar = {:.6f}, rho = {:.6f} \n'.format(-np.rad2deg(phi), sol['fun'], rho))
            
            if polar_coord:
                stop = start + rho*np.array((np.cos(phi), np.sin(phi)))
                stop = np.asarray(np.rint(stop), dtype=int)
        
        else:
            if polar_coord:
                rho, phi = stop[0], -np.deg2rad(stop[1])
                stop = start + rho*np.array((np.cos(phi), np.sin(phi)))
                stop = np.asarray(np.rint(stop), dtype=int)
            else:
                pass
        self.coords = (start, stop, def_kwargs['linewidth'])
        
        # saving variables
        # matrix indeces differ from cartesian coordinates
        measure = profile_line(self.im, start[::-1], stop[::-1], **def_kwargs)
        self.measure = measure
        self.ind = np.arange(0, self.measure.size)
        return
        
        
        
    def calc_dips(self, bg_func = quadratic, prominence=10, height=None, width=None, rel_height=0.9):
        '''Calculate the dips in the profile.
        Find downwards orientated peaks (lower peaks, minimum peaks).
        Calculate the widths from the prominence data.
        Order peaks into pairs.
        Estimate background values.
        Calculate dips from pairs.
        
        Parameters
        ----------
        bg_func : callable, optional
            f(x, *args) -> y
            Scalar function that approximates the background values of the 
            measurement. 'args' will be calculated by regression.
        prominence : float, optional
            Peak prominence of the lower peaks to find. The absolute difference
            between the peak and its contour line.
            Reference in 'scipy.signal.peak_prominence()'.
        height : float, optional
            Peak height of the lower peaks to find (peaks lying above this value 
            will be excluded).
            Reference in 'scipy.signal.find_peaks()'.
        width : float, optional
            Peak width of the lower peaks to find.
            Reference in 'scipy.signal.peak_width()'.
        rel_height : float, optional
            Relative height to measure the width of the peaks at.
            Reference in 'scipy.signal.peak_widths()'.
        
        Returns
        -------
        None
            Write found dips to self.dips : ndarray, shape (m, )
            Write found minima to self.min : array-like, shape (2, n)
            Write found background values to self.bg : array-like, shape (2, o)'''
        
        
        if self.measure is None:
            raise OrderError("no profile was measured")
                
        if height is not None: 
            # self.measure is negated and so is the height
            height = -height
        
        
        # finding minimum peaks
        
        peaks, prop = find_peaks(-self.measure, prominence=prominence, height=height, width=width)
        prominence_data = (prop['prominences'], prop['left_bases'], prop['right_bases'])
        
        if peaks.size == 0:
            # raised when no peaks are found
            raise ResultError('No lower peaks were found')
        
        ### fitting background
        # cutting out peaks
        widths, _, _, _= peak_widths(-self.measure, peaks, rel_height=rel_height, prominence_data=prominence_data)
        
        mask_bg = np.ones(self.ind.shape, dtype=bool)
        intervals = ()
        for i in range(len(widths)):
            # limits should not exceed the array limits 0...len
            lim1 = int(np.rint(peaks[i]-widths[i]))
            if lim1 < 0:
                lim1 = 0
            elif lim1 >= self.measure.size:
                lim1 = self.measure.size-1
            
            lim2 = int(np.rint(peaks[i]+widths[i]))
            if lim2 < 0:
                lim2 = 0
            elif lim2 >= self.measure.size:
                lim2 = self.measure.size-1
            
            mask_bg[lim1 : lim2+1] = False
        
        # fitting
        popt_bg, _ = curve_fit(bg_func, self.ind[mask_bg], self.measure[mask_bg])
        bg_val = bg_func(self.ind, *popt_bg)
        
        
        # ordering minima into pairs.
        dist = peaks[1:] - peaks[:-1]  # pairwise distances: #1-#0, #2-#1, #3-#2, etc.
        dist_max = 1.05*dist[0]
        
        dips = []
        max_peaks = []
        for i in range(len(dist)):
            
            # end of array
            #if i == len(dist)-1:
            #    dips.append(0)
                
            # a pair
            if dist[i] <= dist_max: #weakest point of algorithm!!!
                lim1, lim2 = peaks[i], peaks[i+1]
                
                # the minima
                A = abs(bg_val[lim1] - self.measure[lim1])
                B = abs(bg_val[lim2] - self.measure[lim2])
                
                # the intermediate maximum
                C_pos = np.argmax(self.measure[lim1: lim2+1]) + lim1
                max_peaks.append(C_pos)
                C = abs(bg_val[C_pos] - self.measure[C_pos])
                
                dips.append(100 * (A + B - 2*C) / (A + B))
                
            # segmentations of Non-Pairs are not needed

        while len(dips) < len(self.wire_spacing):
            dips.append(0)
        
        # saving variables
        self.dips = np.array(dips)
        
        self.peaks = (peaks, self.measure[peaks])
        self.max_peaks = (max_peaks, self.measure[max_peaks])
        self.despike = (self.ind[mask_bg], self.measure[mask_bg])
        self.bg = (self.ind, bg_val)
        return
    
    
    def interpolate(self):
        '''Choose important dips, perform quadradic interpolation and
        choose the correct solution for the 20% dip.
        
        Returns
        -------
        None
            Write result to self.dip20 : scalar
            Write found dips of interest to self.dipsoi : array-like, shape (2, n)
            Write interpolation fit to self.inter : array-like, shape (2, m)'''
        
        if self.dips is None:
            raise OrderError("'segmentation()' has to be executed at least once")
        
        if self.dips.size == 0:
            raise ResultError('no dips were found')
        
        if all(self.dips < 20):
            raise ResultError("no dip > 20 was found, interpolation not possible")
        elif all(self.dips > 20):
            warn('no dip < 20 was found, a dip of 0 is used instead')

        # clean up:
        use_dips = copy.deepcopy(self.dips)
        use_wire_spacing = copy.deepcopy(self.wire_spacing)
        print("  {n} Dips found.".format(n=len(self.dips)))
        success = False
        while not success:
            success = True
            i = 1
            while i < len(use_dips):
                print("  Dip {i}: {d}".format(i=i, d=use_dips[i]))
                
                # The following, more shallow dip must not be deeper by more than 5% of previous dip.
                # Otherwise, there seems to be an order problem, i.e. the dips should become more shallow
                # with each iteration. 
                if (use_dips[i] - use_dips[i-1]) > 5:  
                    # This dip seems to be invalid. Drop it.
                    print("  Dropping {idx}.".format(idx=(i-1)))
                    use_dips = np.delete(use_dips, i-1)
                    use_wire_spacing = np.delete(use_wire_spacing, i-1)
                    i = i - 1
                    success = False

                i = i + 1


        # finding the 20% crossing element and choosing the important neighboring values
        index = np.arange(0, (len(use_dips)))
        for i in index:
            lower = i-2   # including
            upper = i+2   # including

            while lower < 0:
                lower += 1
            while upper >= len(use_dips):
                upper -= 1

            # Do not include nearest or next-nearest neighbor dips
            # with less than 1.5% modulation depth,
            # as this will give unpleasant fits:
            if (upper-i) >= 2:
                if use_dips[upper] < 1.5:  # %
                    upper -= 1

            if (upper-i) >= 1:
                if use_dips[upper] < 1.5:  # %
                    upper -= 1

            pos = index[lower:(upper+1)]

            if i < len(use_dips)-1:
                if use_dips[i+1] < 20:
                    self.criticalIndex = i
                    break
            else:
                self.criticalIndex = i
                break
        
        dists = use_wire_spacing[pos] #* self.scale
        dips = use_dips[pos]
        popt, _ = curve_fit(Interpolation.quadratic, dists, dips)
        self.a = popt[0]
        self.b = popt[1]
        self.c = popt[2]
        
        dips20 = Interpolation.inverted_quadratic(20, *popt)
        #for i in dips20:
        #    if i>= min(dists) and i<=max(dists):
        #        self.dip20 = i

        # Depending on the curvature of the interpolation function,
        # we choose either the left or the right zero-crossing of the quadratic function
        # as the iSRb:
        if self.a >= 0:
            self.dip20 = max(dips20)
        else:
            self.dip20 = min(dips20)
        
        if isinstance(self.dip20, type(None)):
            warn("could not estimate a single 20%-dip. The found values are {0:.4f} and {1:.4f}".format(*dips20))
        
        self.dipsoi = (dists, dips)
        lowerBound = min(dists)
        lowerBound = min(lowerBound, self.dip20)
        upperBound = max(dists)
        upperBound = max(upperBound, self.dip20)
        x = np.linspace(lowerBound, upperBound, 100)
        self.inter = (x, Interpolation.quadratic(x, *popt))
        return