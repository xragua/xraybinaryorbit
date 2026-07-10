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

# CONIC ORBIT #############################################################################
def doppler_orbit_theoretical(
    t,
    units="keV",
    show_plot=False,
    precision_for_phase=0.01,
    load_directly=False,
    parameter_list=None,
):
    """Compute the Doppler modulation produced by the compact object's orbit.

    The package convention is assumed to be:

    - ``semimajor``: relative semimajor axis, ``a_rel = a1 + a2``, in units
      of the companion-star radius ``Rstar``.
    - ``Mstar1``: mass of the compact object that carries the observed
      feature.
    - ``Mstar2``: mass of the companion star.
    - ``periapsis`` and ``inclination``: angles in degrees.

    The compact object's barycentric semimajor axis is therefore

        a1 = a_rel * Mstar2 / (Mstar1 + Mstar2).

    For an eccentric orbit, the line-of-sight orbital velocity contains both
    the radial polar component, ``dR/dt``, and the tangential component,
    ``R*dtheta/dt``. The angular velocity returned by
    ``_orbital_time_to_phase`` is used for both terms.

    Parameters
    ----------
    t : array-like
        Time array in seconds.
    units : {"keV", "s", "angstrom"}, optional
        Units of the returned shifted feature.
    show_plot : bool, optional
        If True, save and display the orbital/Doppler diagnostic plot.
    precision_for_phase : float, optional
        Phase resolution passed to ``_orbital_time_to_phase``.
    load_directly : bool, optional
        Passed to ``_manage_parameters``.
    parameter_list : sequence, optional
        Parameter values passed to ``_manage_parameters``.

    Returns
    -------
    t : numpy.ndarray
        Input times in seconds.
    x : numpy.ndarray
        Orbital phase corresponding to ``t``.
    equation : numpy.ndarray
        Doppler-shifted feature in the requested units.
    """

    parameter_names = [
        "iphase",
        "semimajor",
        "orbitalperiod",
        "eccentricity",
        "periapsis",
        "inclination",
        "Rstar",
        "Mstar1",
        "Mstar2",
        "wind_vel",
        "feature",
    ]

    fixed_values = _manage_parameters(
        parameter_names,
        "orbit",
        load_directly=load_directly,
        parameter_list=parameter_list,
    )

    (
        iphase,
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        inclination,
        Rstar,
        Mstar1,
        Mstar2,
        wind_vel,
        feature,
    ) = fixed_values

    t = np.asarray(t, dtype=float)

    if t.ndim != 1 or t.size == 0:
        raise ValueError("t must be a non-empty one-dimensional array.")
    if not np.all(np.isfinite(t)):
        raise ValueError("t contains non-finite values.")
    if units not in {"keV", "s", "angstrom"}:
        raise KeyError(
            f"Invalid unit: {units}. Valid units are 'keV', 's', or "
            "'angstrom'."
        )
    if not 0 <= eccentricity < 1:
        raise ValueError("eccentricity must satisfy 0 <= eccentricity < 1.")
    if Mstar1 <= 0 or Mstar2 <= 0:
        raise ValueError("Mstar1 and Mstar2 must be positive.")
    if Rstar <= 0 or semimajor <= 0 or orbitalperiod <= 0:
        raise ValueError("Rstar, semimajor, and orbitalperiod must be positive.")
    if precision_for_phase <= 0:
        raise ValueError("precision_for_phase must be positive.")

    # Preserve the original feature convention used by this module.
    rest_feature = {
        "keV": kev_ams / feature,
        "s": feature,
        "angstrom": feature,
    }[units]

    # x is the orbital phase in cycles and W = dtheta/dt in rad s^-1.
    x, _, W = _orbital_time_to_phase(
        t,
        iphase,
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        Rstar,
        Mstar1,
        Mstar2,
        precision=precision_for_phase,
    )

    x = np.asarray(x, dtype=float)
    W = np.asarray(W, dtype=float)

    # ``semimajor`` is the relative semimajor axis a_rel = a1 + a2.
    # Mstar1 is the compact object/emitter, so its barycentric semimajor axis is
    # a1 = a_rel * Mstar2 / (Mstar1 + Mstar2).
    a_emitter = semimajor * Mstar2 / (Mstar1 + Mstar2)

    theta = 2.0 * np.pi * x
    omega = np.deg2rad(periapsis)
    inclination_rad = np.deg2rad(inclination)
    true_anomaly = theta - omega

    # Barycentric orbital radius of the compact object, in units of Rstar.
    p_emitter = a_emitter * (1.0 - eccentricity**2)
    denominator = 1.0 + eccentricity * np.cos(true_anomaly)
    R = p_emitter / denominator

    # dR/dnu, also in units of Rstar per radian. Since omega is constant,
    # dnu/dt = dtheta/dt = W.
    dR_dnu = (
        p_emitter
        * eccentricity
        * np.sin(true_anomaly)
        / denominator**2
    )
    R_dot = dR_dnu * W  # Rstar s^-1

    # Projection consistent with the original expression, corresponding to
    # z = R*cos(theta)*sin(i). Its complete derivative is
    # dz/dt = [R_dot*cos(theta) - R*W*sin(theta)]*sin(i).
    stellar_radius_m = Rstar * rsun_m
    v_orbital_los = (
        R_dot * np.cos(theta) - R * W * np.sin(theta)
    ) * stellar_radius_m * np.sin(inclination_rad)

    # Preserve the original optional wind/outflow contribution and sign
    # convention.
    v_wind_los = (
        wind_vel
        * 1000.0
        * np.cos(theta)
        * np.sin(inclination_rad)
    )

    v_los = v_orbital_los + v_wind_los
    doppler_factor = (c + v_los) / c

    equation = {
        "keV": kev_ams / (rest_feature * doppler_factor),
        "s": rest_feature * doppler_factor,
        "angstrom": rest_feature * doppler_factor,
    }[units]

    if show_plot:
        fig = plt.figure(figsize=(15, 5))
        ax1 = fig.add_subplot(1, 3, 1)
        ax2 = fig.add_subplot(1, 3, 2)
        ax3 = fig.add_subplot(1, 3, 3, projection="polar")

        ax1.plot(t, equation, label="Model")
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("Expected feature evolution")
        ax1.legend()

        ax2.plot(x, equation, label="Model")
        ax2.set_xlabel("Orbital phase")
        ax2.set_ylabel("Expected feature evolution")
        ax2.legend()

        ax3.plot(theta, R)
        ax3.set_theta_zero_location("N")
        ax3.set_theta_direction(-1)
        ax3.set_rlabel_position(90)
        ax3.set_title("Emitter orbit around the barycentre")

        fig.tight_layout()
        fig.savefig("orbit_doppler_evolution.png", dpi=200, bbox_inches="tight")
        plt.show()

    return t, x, equation

    
# SPIRAL #########################################################################################
def doppler_spiral_theoretical(t, units="keV", show_plot=False, load_directly=False, parameter_list=None, verbose_complete=False):
    """
    Computes the Doppler variation expected from a spiral movement given a time array in seconds.

    A logarithmic spiral is a type of spiral that grows by a constant factor with each turn. The spiral equation
    in polar coordinates is r = a * e^(b * θ), where:
    
    - r is the distance from the origin (radius)
    - θ is the angle from a reference direction (usually the positive x-axis)
    - a is the scale factor determining how quickly the spiral grows
    - b is the rate of rotation controlling the tightness or looseness of the spiral

    Parameters
    ----------
    t : array-like
        Time array in seconds.
    units : str, optional
        Units for the output Doppler variation. Default is "keV". Options include:
        - "keV": Doppler variation in keV.
        - "s": Doppler variation in seconds.
        - "angstrom": Doppler variation in angstroms.
    show_plot : bool, optional
        If True, displays and saves a plot of the spiral and Doppler evolution. Default is False.

    Notes
    -----
    A form will appear to input the necessary orbital parameters. These parameters are saved in a .txt file
    in the current directory and automatically loaded in subsequent runs, avoiding the need to re-enter parameters.
    You can modify only those parameters that require adjustment.

    Returns
    -------
    x : array-like
        Orbital phase array corresponding to the input time array.
    equation : array-like
        Expected Doppler variation computed for the given spiral movement.

    """

    parameter_names=["iphase_spiral", "semimajor_spiral", "b", "omega", "inclination_spiral", "feature"]
    
    fixed_values = _manage_parameters(parameter_names, "spiral", load_directly=load_directly,parameter_list=parameter_list, verbose_complete=verbose_complete )
    iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature = fixed_values
    #...................................................
    feature_ = {
        "keV": kev_ams/feature,
        "s": feature,  # No conversion needed for seconds
        "amstrong": feature  # No conversion needed for lambda
    }
    
    # Raise KeyError if the unit is invalid
    if units not in feature_:
        raise KeyError(f"Invalid unit: {units}. Valid units are 'keV', 's', or 'angstrom'.")
    
    feature = feature_.get(units, 1)
    
    t0 = min(t)
    #...................................................

    x = (t-t0) * omega + iphase_spiral
    R = semimajor_spiral * np.exp(b * 2 * np.pi * x)

    vdop = -R * rsun_m * omega * np.sin(2 * np.pi * x ) * np.sin(2 * np.pi * inclination_spiral/360)
    
    #...................................................
    vdop =  np.array(vdop)
    
    equation_ = {
        "keV": kev_ams/(feature * (c + vdop) / c),
        "s": (feature * (c + vdop) / c),  # No conversion needed for seconds
        "amstrong": (feature * (c + vdop) / c)  # No conversion needed for lambda
    }
    
    equation = equation_.get(units, 1)
    #...................................................
    if show_plot:
    
        fig, (ax1,ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
        
        ax1.errorbar(t, equation,  label='Data', color='b')
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("Expected emission feature evolution")
        ax1.legend()

        ax2.errorbar(x, equation,  label='Data', color='b')
        ax2.set_xlabel("Orbital phase")
        ax2.set_ylabel("Expected emission feature evolution")
        ax2.legend()

        ax3 = plt.subplot(1, 3, 3, projection='polar')
        ax3.plot(x[R>0] * 2 * np.pi, R[R>0],"b")
        ax3.set_theta_zero_location("N")
        ax3.set_theta_direction(-1)
        ax3.set_rlabel_position(90)
        
        plt.savefig("spiral_doppler_evolution.png")

        plt.tight_layout()
        
        del units
    
    return t,x,equation, R
    
# ORBIT IN ORBIT #####################################################################################
def doppler_disc_theoretical(t, units="keV", show_plot=False, load_directly=False, parameter_list=None):
    """
    Computes the Doppler variation expected from orbital movement in a main orbit, assuming a ballistic
    movement of plasma around a compact object or the movement of a mass entering an accretion disc.

    Parameters
    ----------
    t : array-like
        Time array in seconds.
    units : str, optional
        Units for the output Doppler variation. Default is "keV". Options include:
        - "keV": Doppler variation in keV.
        - "s": Doppler variation in seconds.
        - "angstrom": Doppler variation in angstroms.
    show_plot : bool, optional
        If True, displays and saves a plot of the disc and Doppler evolution. Default is False.

    Notes
    -----
    A form will appear to input the necessary orbital parameters. These parameters are saved in a .txt file
    in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter parameters,
    and allows modification of only those parameters that require adjustment.

    Returns
    -------
    t : array-like
        The input time array.
    x : array-like
        Orbital phase array for the first orbit.
    x2 : array-like
        Orbital phase array for the second orbit.
    equation : array-like
        Expected Doppler variation computed for the given orbital and disc movement.
    """

    parameter_names=["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis", "inclination", "Rstar", "Mstar1", "Mstar2", "iphase2", "semimajor2", "orbitalperiod2", "eccentricity2", "periapsis2", "inclination2",  "Mass3","wind_vel", "feature"]
    
    fixed_values = _manage_parameters(parameter_names, "disc",load_directly=load_directly,parameter_list=parameter_list )
    iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2, inclination2,  Mass3, wind_vel, feature = fixed_values

    #...................................................
    feature_ = {
        "keV": kev_ams/feature,
        "s": feature,  # No conversion needed for seconds
        "amstrong": feature  # No conversion needed for lambda
    }
    # Raise KeyError if the unit is invalid
    if units not in feature_:
        raise KeyError(f"Invalid unit: {units}. Valid units are 'keV', 's', or 'angstrom'.")
    
    feature = feature_.get(units, 1)
    
    t0 = min(t)
    #...................................................
    
    abar = semimajor * max(Mstar1,Mstar2)/(Mstar1+Mstar2)
    orbitalperiod_s = orbitalperiod * 24 * 60 * 60

    x,_,W  = _orbital_time_to_phase(t ,iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, max(Mstar1,Mstar2), min(Mstar1,Mstar2)+Mass3, precision=0.01)
    
    R = (abar * (1 - eccentricity ** 2) / (1 + eccentricity * np.cos((x - periapsis ) * 2 * np.pi)))
    vdop1 =  -R *Rstar * rsun_m * W * np.sin(2 * np.pi * x ) * np.sin(2 * np.pi * inclination /360)
    v_rad =   wind_vel*1000 * np.cos(2 * np.pi * x ) * np.sin(2 * np.pi * inclination /360)
    
    #...................................................
    
    abar2 = semimajor2* min(Mstar1,Mstar2)/(min(Mstar1,Mstar2)+Mass3)
    orbitalperiod_s2 = orbitalperiod2 * 24 * 60 * 60

    x2,_,W2  = _orbital_time_to_phase(t ,iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2, Rstar, min(Mstar2,Mstar1), Mass3, precision=0.01)
    #.........................................................
    w2 = periapsis2 #+ (orbitalperiod/orbitalperiod2)*(x2-iphase2)+iphase #CORRECT PERIAPSIS IN 2 ORBIT
     
    R2 = abar2 * (1 - eccentricity2 ** 2) / (1 + eccentricity2 * np.cos((x2 - w2) * 2 * np.pi))
    vdop2 = -R2 * Rstar * rsun_m * W2 * np.sin(2 * np.pi * x2 ) * np.sin(2 * np.pi * inclination2/360)
    
    #...................................................
    
    vdop =  vdop2 + vdop1 +v_rad
    
    equation_ = {
        "keV": kev_ams/(feature * (c + vdop) / c),
        "s": (feature * (c + vdop) / c),  # No conversion needed for seconds
        "amstrong": (feature * (c + vdop) / c)  # No conversion needed for lambda
    }

    equation = equation_.get(units, 1)
    #...................................................
    if show_plot:
    
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

        ax1.errorbar(t, equation,  label='Data', color='b')
        ax1.set_xlabel("Time")
        ax1.set_ylabel("Expected emission feature evolution")
        ax1.legend()

        ax2.errorbar(x, equation,  label='Data', color='b')
        ax2.set_xlabel("Orbital phase")
        ax2.set_ylabel("Expected emission feature evolution")
        ax2.legend()
        
        ax3 = plt.subplot(1, 3, 3, projection='polar')

        ax3.plot(x * 2 * np.pi, R+R2,"b")
        ax3.set_theta_zero_location("N")
        ax3.set_theta_direction(-1)
        ax3.set_rlabel_position(90)
        
        plt.savefig("disc_doppler_evolution.png")

        plt.tight_layout()
    #...................................................
        del units
    
    return t,x,x2,equation
    
# SPIRAL IN ORBIT ####################################################################################
def doppler_spiral_in_orbit_theoretical(t, units="keV", show_plot=False, load_directly=False, parameter_list=None):
    """
    Computes the Doppler variation expected from an orbital movement with a logarithmic spiral component,
    given a time array in seconds.

    Parameters
    ----------
    t : array-like
        Time array in seconds.
    units : str, optional
        Units for the output Doppler variation. Default is "keV". Options include:
        - "keV": Doppler variation in keV.
        - "s": Doppler variation in seconds.
        - "angstrom": Doppler variation in angstroms.
    show_plot : bool, optional
        If True, displays and saves a plot of the spiral and Doppler evolution. Default is False.

    Notes
    -----
    A form will appear to input the necessary orbital parameters. These parameters are saved in a .txt file
    in the current directory and automatically loaded in subsequent runs, avoiding the need to re-enter
    all parameters. Only those parameters that require adjustment need to be modified.

    Returns
    -------
    x : array-like
        Orbital phase array corresponding to the input time array.
    equation : array-like
        Expected Doppler variation computed for the given orbital movement with the spiral component.

    Description
    -----------
    The logarithmic spiral grows in size by a constant factor with each turn. The spiral is described by
    the polar coordinates equation:
    r = a * e^(b * θ), where:
    - r is the distance from the origin (radius)
    - θ is the angle from a reference direction (usually the positive x-axis)
    - a is the scale factor that determines how quickly the spiral grows.
    - b is the rate of rotation, controlling the tightness or looseness of the spiral.
    """


    parameter_names=["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis", "inclination", "iphase_spiral", "semimajor_spiral", "b", "omega", "inclination_spiral", "Rstar", "Mstar1", "Mstar2", "wind_vel", "feature"]
    
    fixed_values = _manage_parameters(parameter_names, "spiral_in_orbit",load_directly=load_directly,parameter_list=parameter_list )
    iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, Rstar, Mstar1, Mstar2, wind_vel, feature = fixed_values

    #...................................................
    feature_ = {
        "keV": kev_ams/feature,
        "s": feature,  # No conversion needed for seconds
        "amstrong": feature  # No conversion needed for lambda
    }
    # Raise KeyError if the unit is invalid
    if units not in feature_:
        raise KeyError(f"Invalid unit: {units}. Valid units are 'keV', 's', or 'angstrom'.")
        
    feature = feature_.get(units, 1)
    t0 = min(t)
    #...................................................
    
    abar = semimajor* max(Mstar1,Mstar2)/(Mstar1+Mstar2)
    orbitalperiod_s = orbitalperiod * 24 * 60 * 60
    
    x,_,W = _orbital_time_to_phase(t , iphase,semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)

    R = abar * (1-eccentricity ** 2)/(1+eccentricity * np.cos((x-periapsis/360) * 2 * np.pi)) #R of ellipse
       
    vdop =  -R * Rstar * rsun_m * W * np.sin(2 * np.pi * x ) * np.sin(2 * np.pi * inclination /360)
    v_rad =   wind_vel*1000 * np.cos(2 * np.pi * x ) * np.sin(2 * np.pi * inclination /360)
    #...................................................

    x2 = (t-t0) * omega + iphase_spiral
    R2 = semimajor_spiral * np.exp(b * 2 * np.pi * x2)

    vdop2 = -R2 * Rstar * rsun_m * omega * np.sin(2 * np.pi * x2 ) * np.sin(2 * np.pi * inclination_spiral/360)
    
    #...................................................
    vdop1 = np.array(vdop)
    vdop2 = np.array(vdop2)
    
    vdop =  vdop2 + vdop1 + v_rad

    equation_ = {
        "keV": kev_ams/(feature * (c + vdop) / c),
        "s": (feature * (c + vdop) / c),  # No conversion needed for seconds
        "amstrong": (feature * (c + vdop) / c)  # No conversion needed for lambda
    }
    
    equation = equation_.get(units, 1)
    #...................................................
    if show_plot:
    
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 5))
        
        ax1.errorbar(t, equation,  label='Data', color='b')
        ax1.set_xlabel("Time")
        ax1.set_ylabel("Expected emission feature evolution")
        ax1.legend()
        
        ax2.errorbar(x, equation,  label='Data', color='b')
        ax2.set_xlabel("Orbital phase")
        ax2.set_ylabel("Expected emission feature evolution")
        ax2.legend()

        ax2 = plt.subplot(1, 3, 3, projection='polar')
        ax2.plot(x[R>0] * 2 * np.pi, R2[R>0]+R[R>0],"b")
        ax2.set_theta_zero_location("N")
        ax2.set_theta_direction(-1)
        ax2.set_rlabel_position(90)

        plt.tight_layout()
        
        plt.savefig("spiral_in_orbit_doppler_evolution.png")
    #...................................................
        del units
        
    return t,x,x2,equation
     
