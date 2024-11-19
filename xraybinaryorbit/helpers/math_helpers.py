import numpy as np
import pandas as pd
from scipy.integrate import quad
from pyswarm import pso
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import warnings
import tkinter as tk
from tkinter import messagebox
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks, peak_widths, peak_prominences, find_peaks_cwt
from scipy.integrate import quad
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import inspect
import math

def _gaussian(x, mean, sigma):
    """
    Compute the Gaussian (normal) distribution for the given input.

    This function evaluates the Gaussian (normal) distribution for a given
    point `x`, with a specified mean and standard deviation (sigma).

    Parameters
    ----------
    x : float or array-like
        The point(s) at which to evaluate the Gaussian function.
    mean : float
        The mean (center) of the Gaussian distribution.
    sigma : float
        The standard deviation (width) of the Gaussian distribution.

    Returns
    -------
    float or array-like
        The value of the Gaussian function at `x`, based on the specified mean and sigma.
    
    Notes
    -----
    The Gaussian function is given by:
    .. math::
        f(x) = \\frac{1}{\\sigma \\sqrt{2 \\pi}} \\exp \\left(-\\frac{(x - \\mu)^2}{2 \\sigma^2}\\right)
    where `mu` is the mean and `sigma` is the standard deviation.

    """
    
    return np.exp(-(x - mean)**2 / (2 * sigma**2)) / (sigma * np.sqrt(2 * np.pi))

#...................................................  Prepare timebins for fitting functions
def _time_pairs(edges):
    """
    Create consecutive pairs of values from the input edges array.

    This function takes an array of edges (e.g., time bin edges) and creates
    consecutive pairs of adjacent elements. The result is a list of tuples, where
    each tuple contains two consecutive values from the edges array.

    Parameters
    ----------
    edges : array-like
        A list or array of values representing the edges. These could be time bin
        edges or any other set of ordered values.

    Returns
    -------
    pairs_ : list of tuples
        A list of consecutive pairs of the form (edge[i], edge[i+1]) from the input `edges`.
    """
    edge_pairs = [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]
    pairs_ = []

    for pair in edge_pairs:
        pairs_.append(pair)

    return pairs_

#...................................................PCHIP (Piecewise Cubic Hermite Interpolating Polynomial) is a method used for interpolation, particularly for interpolating data points with unevenly spaced independent variable (x) values.
def _interpolate_pchip(tmin, tmax, tinterval, t, y, sy):
    """
    Interpolate data using Piecewise Cubic Hermite Interpolating Polynomial (PCHIP).

    This function uses the PCHIP method to interpolate both the data and its associated
    uncertainties (errors) over a new range of time values. PCHIP is ideal for maintaining
    the shape of data and avoiding oscillations that other interpolation methods may introduce.

    Parameters
    ----------
    tmin : float
        The minimum value of the new time grid for interpolation.
    tmax : float
        The maximum value of the new time grid for interpolation.
    tinterval : float
        The interval between consecutive time points in the new time grid.
    t : array-like
        The original time data points at which the data `y` and its uncertainties `sy` are defined.
    y : array-like
        The data values corresponding to the time points `t`.
    sy : array-like
        The uncertainties (errors) associated with the data points `y`.

    Returns
    -------
    t_new : array-like
        The new time points from `tmin` to `tmax` with step size `tinterval`.
    new_data : array-like
        The interpolated data values at the new time points `t_new`.
    snew_data : array-like
        The interpolated uncertainty values at the new time points `t_new`.

    """
    
    t_new = np.arange(tmin, tmax, tinterval)
    
    # Interpolate the data using PCHIP
    cs = PchipInterpolator(t, y)
    scs = PchipInterpolator(t, sy)
    
    new_data = cs(t_new)
    snew_data = scs(t_new)
    
    return t_new, new_data, snew_data



#................................................... Weighted to the error chi squared
def _chi_squared_weighted(y_data, y_err, y_pred):
    """
    Compute the weighted chi-squared statistic.

    This function calculates the chi-squared statistic, which is a measure of
    the goodness of fit between observed data (`y_data`) and predicted data
    (`y_pred`), weighted by the uncertainties (`y_err`). The lower the chi-squared
    value, the better the fit of the model to the data.

    Parameters
    ----------
    y_data : array-like
        The observed data (dependent variable).
    y_err : array-like
        The uncertainties (errors) associated with the observed data.
    y_pred : array-like
        The predicted data (from a model or fit).

    Returns
    -------
    chi_squared : float
        The weighted chi-squared statistic, which quantifies how well the predicted
        data matches the observed data, taking into account the uncertainties.
    """
    y_data = np.asarray(y_data)
    y_err = np.asarray(y_err)
    y_pred = np.asarray(y_pred)

    # Calculate the weighted chi-squared statistic
    chi_squared = np.sum(((y_data - y_pred) / y_err) ** 2)
    
    return chi_squared

    
def _chi_squared(y_data, y_pred):
    """
    Compute the unweighted chi-squared statistic.

    This function calculates the chi-squared statistic, which is a measure of
    the goodness of fit between observed data (`y_data`) and predicted data
    (`y_pred`) without taking into account any uncertainties. This is a simpler
    version of the chi-squared test where all uncertainties are assumed to be equal.

    Parameters
    ----------
    y_data : array-like
        The observed data (dependent variable).
    y_pred : array-like
        The predicted data (from a model or fit).

    Returns
    -------
    chi_squared : float
        The unweighted chi-squared statistic, which quantifies how well the predicted
        data matches the observed data.
    """

    y_data = np.asarray(y_data)
    y_pred = np.asarray(y_pred)

    # Calculate the unweighted chi-squared statistic
    chi_squared = np.sum((y_data - y_pred) ** 2)
    
    return chi_squared


#................................................... scale data for plotting purposes
def scale(x, y):
    """
    Scale the `x` data to match the range of the `y` data. The purpose is facilitate creation of plots.

    This function scales the values in the `x` array to match the range of the `y` array.
    It linearly transforms the `x` values such that they span the same range as `y`.

    Parameters
    ----------
    x : array-like
        The input data to be scaled.
    y : array-like
        The data whose range is used for scaling `x`.

    Returns
    -------
    x_new : array-like
        The scaled version of `x`, with values transformed to the range of `y`.
    """

    x_new = ((max(y) - min(y)) / (max(x) - min(x))) * (x - max(x)) + max(y)
    return x_new



def _orbital_phase_to_time(ph, iphase, semimajor, orbitalperiod, eccentricity, periapsis,
                           Rstar, Mstar1, Mstar2, precision=0.01):
    """
    Converts an orbital phase array to a time array for a compact object orbiting a companion star.

    The compact object moves faster at periastron than at apoastro due to the conservation of angular momentum
    and Kepler's second law of planetary motion. This law states that the object sweeps out equal areas in
    equal times, requiring higher velocity when the object is closer to the star (periapsis).

    Parameters
    ----------
    ph : array-like
        Orbital phase array, representing the phase of the object in its orbit.
    iphase : float
        The initial orbital phase.
    semimajor : float
        Semi-major axis of the orbit.
    orbitalperiod : float
        Orbital period of the compact object around the companion star.
    eccentricity : float
        Orbital eccentricity, describing how elongated the orbit is.
    periapsis : float
        Argument of periapsis, representing the orientation of the elliptical orbit.
    Rstar : float
        Radius of the companion star.
    Mstar1 : float
        Mass of the compact object.
    Mstar2 : float
        Mass of the companion star.
    precision : float, optional
        Resolution for the phase array. Default is 0.01.

    Returns
    -------
    ph : array-like
        The input orbital phase array.
    time : array-like
        Time array corresponding to the orbital phase.
    W : array-like
        Angular velocity array corresponding to the orbital phase, accounting for changes due to eccentricity.

    """

    #.............................Load parameters
    if len(ph)>1:
        th_ = np.arange(-precision+iphase, max(ph)+iphase, precision)
        th = np.arange(iphase, max(ph)+iphase, precision)

    else:
        th_ = np.arange(-precision+iphase, 1+ph+iphase, precision)
        th = np.arange(iphase, 1+ph+iphase, precision)
        
    number_of_orbits=max(th_-min(th_))
    
    abar = semimajor * max(Mstar1, Mstar2) / (Mstar1 + Mstar2)
    orbital_period_s = orbitalperiod*24*60*60

    def integrand(theta):
        return 0.5 * (abar * (1 - eccentricity ** 2) / (1 + eccentricity * np.cos((theta - periapsis / 360) * 2 * np.pi))) ** 2

    time_ = []
    tprev = 0
    w_ = []

    for i in range(len(th_) - 1):
        
        area_, _ = quad(integrand, th_[i], th_[i + 1])
        tprev += area_
        time_.append(tprev)
        w_.append(2*np.pi* abs(th_[i+1]- th_[i])/area_)

    constant = number_of_orbits * orbital_period_s/max(time_)
    time = np.array(time_) * constant
    W = np.array(w_) / constant

    time_interpolator = interp1d(th, time, kind='cubic', fill_value="extrapolate")
    time = time_interpolator(ph)
    
    w_interpolator = interp1d(th, W, kind='cubic', fill_value="extrapolate")
    W = w_interpolator(ph)
        
    return ph, time, W
# TIME TO PHASE ###########################################################################

def _orbital_time_to_phase(t, iphase, semimajor, orbitalperiod, eccentricity, periapsis,
                           Rstar, Mstar1, Mstar2, precision=0.01):
    """
    Converts an orbital time array to a phase array for a compact object orbiting a companion star.

    The compact object moves faster at periastron than at apoastro due to the conservation of angular momentum
    and Kepler's second law of planetary motion. According to this law, objects sweep out equal areas in equal times,
    resulting in faster motion when closer to the star (at periapsis).

    Parameters
    ----------
    t : array-like
        Orbital time array.
    iphase : float
        Initial orbital phase.
    semimajor : float
        Semi-major axis of the orbit.
    orbitalperiod : float
        Orbital period of the compact object.
    eccentricity : float
        Orbital eccentricity, representing the elongation of the orbit.
    periapsis : float
        Argument of periapsis, describing the orientation of the elliptical orbit.
    Rstar : float
        Radius of the companion star.
    Mstar1 : float
        Mass of the compact object.
    Mstar2 : float
        Mass of the companion star.
    precision : float, optional
        Resolution for the phase array. Default is 0.01.

    Returns
    -------
    ph : array-like
        Orbital phase array corresponding to the input time array.
    W : array-like
        Angular velocity array corresponding to the orbital phase, adjusted for the varying speed due to orbital eccentricity.

    """
    #.............................Load parameters
    orbital_period_s = 24*60*60*orbitalperiod
        
    if len(t)>1:
        number_of_orbits = (max(t)-min(t))/orbital_period_s+10
    else:
        number_of_orbits = 10
    
    th_ = np.arange(-precision+iphase, number_of_orbits+iphase, precision)
    
    th = np.arange(iphase, number_of_orbits+iphase, precision)
    
    abar = semimajor * max(Mstar1, Mstar2) / (Mstar1 + Mstar2)

    #.............................Decipher relation between th and time
    abar = semimajor * max(Mstar1, Mstar2) / (Mstar1 + Mstar2)
    
    def integrand(theta):
        return 0.5 * (abar * (1 - eccentricity ** 2) / (1 + eccentricity * np.cos((theta - periapsis / 360) * 2 * np.pi))) ** 2

    time_ = []
    tprev = 0
    w_ = []

    for i in range(len(th_) - 1):
        
        area_, _ = quad(integrand, th_[i], th_[i + 1])
        tprev += area_
        time_.append(tprev)
        w_.append(2*np.pi* abs(th_[i+1]- th_[i])/area_)

    constant = number_of_orbits * orbital_period_s/max(time_)
    time = np.array(time_) * constant
    W_to_interpolate  = np.array(w_) / constant
    
    if (len(time) != (len(th))):
        time=np.array(time[0:len(th)])
        
    if (len(W_to_interpolate) != (len(th))):
        W_to_interpolate=np.array(W_to_interpolate[0:len(th)])
                    
    #.............................Now that we know the relation between time, W and phase respectively, interpolate to obtain phase from our input time
    times_to_interpolate = t-min(t)
    
    phase_interpolator = interp1d(time, th, kind='cubic', fill_value="extrapolate")
    phase = phase_interpolator(times_to_interpolate)
    
    w_interpolator = interp1d(time, W_to_interpolate, kind='cubic', fill_value="extrapolate")
    W = w_interpolator(times_to_interpolate)
    
    return phase, t, W
