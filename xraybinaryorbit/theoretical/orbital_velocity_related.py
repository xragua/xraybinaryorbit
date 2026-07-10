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

from ..helpers.data_helpers import _manage_parameters,_define_x_y_sy,_copy_fields, _load_values_to_interface, _manage_parameters,_load_bounds_to_interface, _manage_bounds

from ..helpers.math_helpers import _gaussian,_time_pairs,_interpolate_pchip,_chi_squared_weighted,_chi_squared,_orbital_phase_to_time,_orbital_time_to_phase

c = 299792458

msun = (1.98847*10**30)*1000 #gr
rsun_m = 696340*1000 #
rsun_cm = 696340*1000*100 #cm

kev_ams = 1.23984193

na = 6.02214076*10**23/1.00797
mu = 0.5
mp = 1.67E-24

# PHASE TO TIME ###########################################################################
def orbital_phase_to_time(
    ph,
    precision=0.01,
    load_directly=False,
    parameter_list=None,
):
    """Convert orbital phase to time using the shared helper implementation.

    The spatial scale and component masses are unnecessary for this conversion
    when the orbital period is supplied directly, so the public interface keeps
    only the four parameters required by the calculation.

    Parameters
    ----------
    ph : array-like
        Orbital phase array.
    precision : float, optional
        Resolution of the internal phase grid. Default is 0.01.
    load_directly : bool, optional
        Passed to ``_manage_parameters``.
    parameter_list : sequence, optional
        Values for ``iphase``, ``orbitalperiod``, ``eccentricity`` and
        ``periapsis``.

    Returns
    -------
    ph : numpy.ndarray
        Input orbital phases.
    time : numpy.ndarray
        Time corresponding to each phase, in seconds.
    W : numpy.ndarray
        Orbital angular velocity, in rad s^-1.
    """

    parameter_names = [
        "iphase",
        "orbitalperiod",
        "eccentricity",
        "periapsis",
    ]

    iphase, orbitalperiod, eccentricity, periapsis = _manage_parameters(
        parameter_names,
        "phase_time",
        load_directly=load_directly,
        parameter_list=parameter_list,
    )

    ph = np.atleast_1d(np.asarray(ph, dtype=float))

    return _orbital_phase_to_time(
        ph,
        iphase,
        None,  # semimajor: unused when orbitalperiod is supplied
        orbitalperiod,
        eccentricity,
        periapsis,
        None,  # Rstar: unused
        None,  # Mstar1: unused
        None,  # Mstar2: unused
        precision=precision,
    )


# TIME TO PHASE ###########################################################################
def orbital_time_to_phase(
    t,
    precision=0.01,
    load_directly=False,
    parameter_list=None,
):
    """Convert orbital time to phase using the shared helper implementation.

    The spatial scale and component masses are unnecessary for this conversion
    when the orbital period is supplied directly, so the public interface keeps
    only the four parameters required by the calculation.

    Parameters
    ----------
    t : array-like
        Orbital time array in seconds.
    precision : float, optional
        Resolution of the internal phase grid. Default is 0.01.
    load_directly : bool, optional
        Passed to ``_manage_parameters``.
    parameter_list : sequence, optional
        Values for ``iphase``, ``orbitalperiod``, ``eccentricity`` and
        ``periapsis``.

    Returns
    -------
    phase : numpy.ndarray
        Orbital phase corresponding to each input time.
    t : numpy.ndarray
        Input time array, in seconds.
    W : numpy.ndarray
        Orbital angular velocity, in rad s^-1.
    """

    parameter_names = [
        "iphase",
        "orbitalperiod",
        "eccentricity",
        "periapsis",
    ]

    iphase, orbitalperiod, eccentricity, periapsis = _manage_parameters(
        parameter_names,
        "time_phase",
        load_directly=load_directly,
        parameter_list=parameter_list,
    )

    t = np.atleast_1d(np.asarray(t, dtype=float))

    return _orbital_time_to_phase(
        t,
        iphase,
        None,  # semimajor: unused when orbitalperiod is supplied
        orbitalperiod,
        eccentricity,
        periapsis,
        None,  # Rstar: unused
        None,  # Mstar1: unused
        None,  # Mstar2: unused
        precision=precision,
    )