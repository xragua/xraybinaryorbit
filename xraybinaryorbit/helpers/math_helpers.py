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

# PHASE TO TIME ###########################################################################

def _orbital_phase_to_time(
    ph,
    iphase,
    semimajor,
    orbitalperiod,
    eccentricity,
    periapsis,
    Rstar,
    Mstar1,
    Mstar2,
    precision=0.01,
):
    """
    Converts an orbital phase array to a time array.

    The phase-time relation is calculated using Kepler's second law.
    When the orbital period is supplied directly, the spatial scale of
    the orbit and the component masses do not enter the conversion.

    Parameters
    ----------
    ph : array-like
        Orbital phase array.
    iphase : float
        Initial orbital phase.
    semimajor : float
        Retained for backward compatibility. It is not used in this
        calculation when orbitalperiod is supplied.
    orbitalperiod : float
        Orbital period in days.
    eccentricity : float
        Orbital eccentricity.
    periapsis : float
        Argument of periapsis in degrees.
    Rstar : float
        Retained for backward compatibility. Not used.
    Mstar1 : float
        Retained for backward compatibility. Not used.
    Mstar2 : float
        Retained for backward compatibility. Not used.
    precision : float, optional
        Resolution of the internal phase grid. Default is 0.01.

    Returns
    -------
    ph : array-like
        Input orbital phase array.
    time : array-like
        Time corresponding to each orbital phase, in seconds.
    W : array-like
        Orbital angular velocity corresponding to each phase, in rad s^-1.
    """

    # Internal phase grid
    if len(ph) > 1:
        th_ = np.arange(
            -precision + iphase,
            max(ph) + iphase,
            precision,
        )

        th = np.arange(
            iphase,
            max(ph) + iphase,
            precision,
        )

    else:
        th_ = np.arange(
            -precision + iphase,
            1 + ph + iphase,
            precision,
        )

        th = np.arange(
            iphase,
            1 + ph + iphase,
            precision,
        )

    number_of_orbits = max(th_ - min(th_))
    orbital_period_s = orbitalperiod * 24 * 60 * 60

    # Relative orbital radius. Its constant spatial scale is unnecessary
    # because the swept areas are subsequently normalized to orbitalperiod.
    def integrand(theta):
        relative_radius = (
            (1 - eccentricity**2)
            / (
                1
                + eccentricity
                * np.cos(
                    (theta - periapsis / 360)
                    * 2
                    * np.pi
                )
            )
        )

        return 0.5 * relative_radius**2

    time_ = []
    tprev = 0
    w_ = []

    for i in range(len(th_) - 1):
        area_, _ = quad(
            integrand,
            th_[i],
            th_[i + 1],
        )

        tprev += area_
        time_.append(tprev)

        w_.append(
            2
            * np.pi
            * abs(th_[i + 1] - th_[i])
            / area_
        )

    # Normalize the swept-area coordinate to the supplied orbital period
    constant = (
        number_of_orbits
        * orbital_period_s
        / max(time_)
    )

    time = np.asarray(time_) * constant
    W = np.asarray(w_) / constant

    time_interpolator = interp1d(
        th,
        time,
        kind="cubic",
        fill_value="extrapolate",
    )

    time = time_interpolator(ph)

    w_interpolator = interp1d(
        th,
        W,
        kind="cubic",
        fill_value="extrapolate",
    )

    W = w_interpolator(ph)

    return ph, time, W


# TIME TO PHASE ###########################################################################
def _true_to_eccentric_anomaly_unwrapped(true_anomaly, eccentricity):
    """
    Convert unwrapped true anomaly into unwrapped eccentric anomaly.

    Angles are expressed in radians and may cover multiple orbital cycles.
    """

    true_anomaly = np.asarray(true_anomaly, dtype=float)

    orbit_number = np.floor(
        true_anomaly / (2.0 * np.pi)
    )

    true_anomaly_mod = (
        true_anomaly
        - orbit_number * 2.0 * np.pi
    )

    eccentric_anomaly_mod = np.arctan2(
        np.sqrt(1.0 - eccentricity**2)
        * np.sin(true_anomaly_mod),
        eccentricity
        + np.cos(true_anomaly_mod),
    )

    eccentric_anomaly_mod = np.mod(
        eccentric_anomaly_mod,
        2.0 * np.pi,
    )

    return (
        eccentric_anomaly_mod
        + orbit_number * 2.0 * np.pi
    )


def _eccentric_to_true_anomaly_unwrapped(
    eccentric_anomaly,
    eccentricity,
):
    """
    Convert unwrapped eccentric anomaly into unwrapped true anomaly.

    Angles are expressed in radians and may cover multiple orbital cycles.
    """

    eccentric_anomaly = np.asarray(
        eccentric_anomaly,
        dtype=float,
    )

    orbit_number = np.floor(
        eccentric_anomaly / (2.0 * np.pi)
    )

    eccentric_anomaly_mod = (
        eccentric_anomaly
        - orbit_number * 2.0 * np.pi
    )

    true_anomaly_mod = np.arctan2(
        np.sqrt(1.0 - eccentricity**2)
        * np.sin(eccentric_anomaly_mod),
        np.cos(eccentric_anomaly_mod)
        - eccentricity,
    )

    true_anomaly_mod = np.mod(
        true_anomaly_mod,
        2.0 * np.pi,
    )

    return (
        true_anomaly_mod
        + orbit_number * 2.0 * np.pi
    )


def _solve_kepler_equation_unwrapped(
    mean_anomaly,
    eccentricity,
    tolerance=1e-13,
    max_iterations=50,
):
    """
    Solve Kepler's equation:

        M = E - e sin(E)

    for unwrapped mean anomalies covering any number of orbital cycles.
    """

    mean_anomaly = np.asarray(
        mean_anomaly,
        dtype=float,
    )

    orbit_number = np.floor(
        mean_anomaly / (2.0 * np.pi)
    )

    mean_anomaly_mod = (
        mean_anomaly
        - orbit_number * 2.0 * np.pi
    )

    if eccentricity < 0.8:
        eccentric_anomaly_mod = (
            mean_anomaly_mod.copy()
        )
    else:
        eccentric_anomaly_mod = np.full_like(
            mean_anomaly_mod,
            np.pi,
            dtype=float,
        )

    converged = np.zeros(
        mean_anomaly_mod.shape,
        dtype=bool,
    )

    for _ in range(max_iterations):

        residual = (
            eccentric_anomaly_mod
            - eccentricity
            * np.sin(eccentric_anomaly_mod)
            - mean_anomaly_mod
        )

        derivative = (
            1.0
            - eccentricity
            * np.cos(eccentric_anomaly_mod)
        )

        correction = residual / derivative

        eccentric_anomaly_mod -= correction

        converged = (
            np.abs(correction) < tolerance
        )

        if np.all(converged):
            break

    # Guaranteed fallback for any element for which Newton's
    # method did not reach the requested tolerance.
    if not np.all(converged):

        failed_indices = np.flatnonzero(
            ~converged.ravel()
        )

        flat_mean = mean_anomaly_mod.ravel()
        flat_eccentric = (
            eccentric_anomaly_mod.ravel()
        )

        for index in failed_indices:

            target_mean_anomaly = flat_mean[index]

            def kepler_residual(E):
                return (
                    E
                    - eccentricity * np.sin(E)
                    - target_mean_anomaly
                )

            flat_eccentric[index] = brentq(
                kepler_residual,
                0.0,
                2.0 * np.pi,
                xtol=tolerance,
            )

        eccentric_anomaly_mod = (
            flat_eccentric.reshape(
                mean_anomaly_mod.shape
            )
        )

    return (
        eccentric_anomaly_mod
        + orbit_number * 2.0 * np.pi
    )







def _orbital_phase_to_time(
    ph,
    iphase,
    semimajor,
    orbitalperiod,
    eccentricity,
    periapsis,
    Rstar,
    Mstar1,
    Mstar2,
    precision=0.01,
):
    """
    Convert orbital phase into elapsed time using Kepler's equation.

    The intentional phase convention is preserved:

        time = 0  <->  phase = iphase - precision

    Parameters retained for backward compatibility
    ------------------------------------------------
    semimajor, Rstar, Mstar1 and Mstar2 are not required when
    orbitalperiod is supplied directly.

    Returns
    -------
    ph : numpy.ndarray
        Input orbital phases.

    time : numpy.ndarray
        Elapsed time in seconds relative to iphase - precision.

    W : numpy.ndarray
        Instantaneous angular velocity dnu/dt in rad s^-1.
    """

    ph = np.atleast_1d(
        np.asarray(ph, dtype=float)
    )

    if ph.size == 0:
        raise ValueError(
            "ph must contain at least one value."
        )

    if not np.all(np.isfinite(ph)):
        raise ValueError(
            "ph contains non-finite values."
        )

    if orbitalperiod <= 0:
        raise ValueError(
            "orbitalperiod must be positive."
        )

    if not 0 <= eccentricity < 1:
        raise ValueError(
            "eccentricity must satisfy "
            "0 <= eccentricity < 1."
        )

    if precision <= 0:
        raise ValueError(
            "precision must be positive."
        )

    orbital_period_s = (
        orbitalperiod * 86400.0
    )

    mean_motion = (
        2.0 * np.pi / orbital_period_s
    )

    periapsis_phase = periapsis / 360.0

    # Intentional convention inherited from the previous code.
    reference_phase = iphase - precision

    reference_true_anomaly = (
        2.0
        * np.pi
        * (
            reference_phase
            - periapsis_phase
        )
    )

    reference_eccentric_anomaly = (
        _true_to_eccentric_anomaly_unwrapped(
            reference_true_anomaly,
            eccentricity,
        )
    )

    reference_mean_anomaly = (
        reference_eccentric_anomaly
        - eccentricity
        * np.sin(
            reference_eccentric_anomaly
        )
    )

    true_anomaly = (
        2.0
        * np.pi
        * (
            ph
            - periapsis_phase
        )
    )

    eccentric_anomaly = (
        _true_to_eccentric_anomaly_unwrapped(
            true_anomaly,
            eccentricity,
        )
    )

    mean_anomaly = (
        eccentric_anomaly
        - eccentricity
        * np.sin(eccentric_anomaly)
    )

    time = (
        mean_anomaly
        - reference_mean_anomaly
    ) / mean_motion

    # Exact instantaneous angular velocity:
    #
    # dnu/dt =
    # n (1 + e cos(nu))² / (1 - e²)^(3/2)
    W = (
        mean_motion
        * (
            1.0
            + eccentricity
            * np.cos(true_anomaly)
        )**2
        / (
            1.0 - eccentricity**2
        )**1.5
    )

    return ph, time, W


def _orbital_time_to_phase(
    t,
    iphase,
    semimajor,
    orbitalperiod,
    eccentricity,
    periapsis,
    Rstar,
    Mstar1,
    Mstar2,
    precision=0.01,
):
    """
    Convert elapsed orbital time into orbital phase using Kepler's
    equation.

    The intentional phase convention is preserved:

        min(t)  <->  phase = iphase - precision

    Parameters retained for backward compatibility
    ------------------------------------------------
    semimajor, Rstar, Mstar1 and Mstar2 are not required when
    orbitalperiod is supplied directly.

    Returns
    -------
    phase : numpy.ndarray
        Unwrapped orbital phase.

    t : numpy.ndarray
        Input time array.

    W : numpy.ndarray
        Instantaneous angular velocity dnu/dt in rad s^-1.
    """

    t = np.atleast_1d(
        np.asarray(t, dtype=float)
    )

    if t.size == 0:
        raise ValueError(
            "t must contain at least one value."
        )

    if not np.all(np.isfinite(t)):
        raise ValueError(
            "t contains non-finite values."
        )

    if orbitalperiod <= 0:
        raise ValueError(
            "orbitalperiod must be positive."
        )

    if not 0 <= eccentricity < 1:
        raise ValueError(
            "eccentricity must satisfy "
            "0 <= eccentricity < 1."
        )

    if precision <= 0:
        raise ValueError(
            "precision must be positive."
        )

    orbital_period_s = (
        orbitalperiod * 86400.0
    )

    mean_motion = (
        2.0 * np.pi / orbital_period_s
    )

    elapsed_time = t - np.min(t)

    periapsis_phase = periapsis / 360.0

    # Intentional convention inherited from the previous code.
    reference_phase = iphase - precision

    reference_true_anomaly = (
        2.0
        * np.pi
        * (
            reference_phase
            - periapsis_phase
        )
    )

    reference_eccentric_anomaly = (
        _true_to_eccentric_anomaly_unwrapped(
            reference_true_anomaly,
            eccentricity,
        )
    )

    reference_mean_anomaly = (
        reference_eccentric_anomaly
        - eccentricity
        * np.sin(
            reference_eccentric_anomaly
        )
    )

    mean_anomaly = (
        reference_mean_anomaly
        + mean_motion * elapsed_time
    )

    eccentric_anomaly = (
        _solve_kepler_equation_unwrapped(
            mean_anomaly,
            eccentricity,
        )
    )

    true_anomaly = (
        _eccentric_to_true_anomaly_unwrapped(
            eccentric_anomaly,
            eccentricity,
        )
    )

    phase = (
        periapsis_phase
        + true_anomaly / (2.0 * np.pi)
    )

    W = (
        mean_motion
        * (
            1.0
            + eccentricity
            * np.cos(true_anomaly)
        )**2
        / (
            1.0 - eccentricity**2
        )**1.5
    )

    return phase, t, W