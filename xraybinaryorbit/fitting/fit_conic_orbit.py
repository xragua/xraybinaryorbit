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

#############################################################################
def conic_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, wind_vel, feature, units, method_, extended_binsize):
    """
    Simulate the Doppler and radial velocities in a conic orbital system using either a discrete or extended method. Its a private function called by
    fit_orbit_ps abd fit_orbit_ls.
    
    This function calculates the Doppler shift and radial velocity in an orbital system where the motion follows a conic orbit (usually elliptical).
    The function can calculate these velocities either for discrete time points or extended time bins, depending on the `method_` parameter.
    The output is the Doppler shift or time shift based on the provided units, such as keV, seconds, or wavelength (angstrom).

    Parameters
    ----------
    x_data : array-like
        The input time data (or time bins) to compute the Doppler shift and radial velocities.
    iphase : float
        The initial phase of the orbit, typically in radians or fractions of the orbit period.
    semimajor : float
        The semi-major axis of the orbit, in solar radii or other relevant units.
    orbitalperiod : float
        The orbital period of the system in days.
    eccentricity : float
        The eccentricity of the orbit, ranging from 0 (circular) to 1 (parabolic).
    periapsis : float
        The argument of periapsis, describing the orientation of the orbit in the plane of motion.
    inclination : float
        The inclination of the orbit in degrees, relative to the plane of the sky.
    Rstar : float
        The radius of the star, typically in solar radii.
    Mstar1 : float
        The mass of the primary star (star 1), typically in solar masses.
    Mstar2 : float
        The mass of the secondary object (star 2), typically in solar masses.
    wind_vel : float
        The wind velocity from the star, typically in km/s.
    feature : float
        The spectral feature of interest, either in keV or another unit, that will be Doppler-shifted.
    units : str
        The units for the feature, such as "keV", "s" (seconds), or "amstrong" (angstroms). Determines how the feature is Doppler shifted.
    method_ : str
        The method used for calculation. Can be "extended" (for time bins) or "discrete" (for individual time points).
    extended_binsize : float
        The bin size for extended methods. This is used to average Doppler shifts and velocities over each time bin.

    Returns
    -------
    equation : array-like
        The Doppler-shifted or time-shifted values of the feature, depending on the `units` parameter.
        This can be in keV, seconds, or wavelength (angstroms), depending on the chosen units.

    Notes
    -----
    - The function uses both the Doppler and radial velocity components to calculate the shift in the feature of interest.
    - For the "extended" method, it computes average values over the bins, while for the "discrete" method, it computes values at individual time points.
    - The Doppler shift accounts for both the motion of the star and any wind velocities.
    - The feature is Doppler-shifted according to the motion of the object, and the result is returned in the desired units.

    """
    t = x_data
    
    feature_ = {
        "keV": kev_ams/feature,
        "s": feature,  # No conversion needed for seconds
        "amstrong": feature  # No conversion needed for lambda
    }
    
    feature = feature_.get(units, 1)
    
    orbitalperiod_s = orbitalperiod * 24 * 60 * 60
    abar = semimajor * max(Mstar1,Mstar2) / (Mstar1 + Mstar2)

    vdop_bin=[]
        
    shape_t = t.shape
    t_to_phase = t.reshape(-1)
    
    if method_=="extended":
    
        t_to_phase_puntual = np.mean(t , axis=1)
        #...................................................
        ph_from_t,_,W = _orbital_time_to_phase(t_to_phase ,iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)
        
        phbins = ph_from_t.reshape(shape_t)
        
        size_phase_bin = np.diff(phbins)
        minsizebin = min(size_phase_bin)
        maxph = max(phbins[-1])
        
        phase = np.arange(0,maxph+10,max(minsizebin/10,maxph/100000))
        _,_,W = _orbital_phase_to_time(phase ,iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)
        
        t_to_phase_puntual = np.mean(t , axis=1)
        phase_puntual ,_, W_puntual = _orbital_time_to_phase(t_to_phase_puntual, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)
        
        #...................................................
        R = abar * (1-eccentricity ** 2)/(1+eccentricity * np.cos((phase-periapsis/360) * 2 * np.pi)) #R of ellipse
            
        vdop =  -R * Rstar * rsun_m * W * np.sin(2 * np.pi * phase ) * np.sin(2 * np.pi * inclination /360)
        vrad =  wind_vel * 1000 * np.sin(2 * np.pi * phase) * np.sin(2 * np.pi * inclination /360)

        for i in range(len(phbins)):
    
            if ( size_phase_bin[i] >= extended_binsize):
                vdop_bin.append(np.mean(vdop[(phase >= phbins[i,0]) & (phase <= phbins[i,1])]) + np.mean(vrad[(phase >= phbins[i,0]) & (phase <= phbins[i,1])]))
            
            if ( size_phase_bin[i] < extended_binsize):

                R_puntual = abar * (1-eccentricity ** 2)/(1+eccentricity * np.cos((phase_puntual-periapsis/360) * 2 * np.pi)) #R of ellipse
        
                vdop_puntual =  -R_puntual * Rstar * rsun_m * W_puntual * np.sin(2 * np.pi * phase_puntual) * np.sin(2 * np.pi * inclination /360)
                vrad_puntual = wind_vel * 1000 * np.cos(2 * np.pi * phase_puntual) * np.sin(2 * np.pi * inclination /360)
            
                vdop_bin.append(vdop_puntual[i]+vrad_puntual[i])
                
    #..................................................
    if method_=="discrete":
        
        phase_discrete ,_,W = _orbital_time_to_phase(t ,iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)
            
        R_puntual = abar * (1-eccentricity ** 2)/(1+eccentricity * np.cos((phase_discrete-periapsis/360) * 2 * np.pi)) #R of ellipse

        vdop_discrete =  -R_puntual * Rstar * rsun_m * W * np.sin(2 * np.pi * phase_discrete) * np.sin(2 * np.pi * inclination /360)
        vrad_discrete =  -wind_vel * 1000 * np.cos(2 * np.pi * phase_discrete) * np.sin(2 * np.pi * inclination /360)
            
        vdop_bin = vdop_discrete + vrad_discrete
    #...................................................

    vdop = np.array(vdop_bin)
    
    equation_ = {
        "keV": kev_ams/(feature * (c + vdop) / c),
        "s": (feature * (c + vdop) / c),  # No conversion needed for seconds
        "amstrong": (feature * (c + vdop) / c)  # No conversion needed for lambda
    }
    equation = equation_.get(units, 1)

    return equation

# PS FIT------------------------------------------------------------------------------------------------------

def fit_orbit_ps(x_data, y_data, y_err=0, num_iterations=3, maxiter=1000, swarmsize=100,
                 units="keV", method_="extended", extended_binsize=0.01,load_directly=False, bound_list=None):
    """
    Fits observed orbital modulation data by estimating parameters such as phase, semi-major axis,
    orbital period, eccentricity, inclination, and periapsis using particle swarm optimization (PSO).

    The fitting process employs PSO to iteratively improve parameter estimates by minimizing the
    chi-squared difference between observed and predicted data.

    The function supports two fitting methods:
    
    - "discrete": Suitable for discrete data points (e.g., spectra with small orbital phase ranges; faster).
    - "extended": Suitable for data with varying or extended bin sizes, typical with current instrument
      resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Parameters
    ----------
    x_data : array-like
        Time bins of the observed data.
    y_data : array-like
        Observed data points corresponding to each time bin.
    y_err : array-like, optional
        Error associated with each observed data point. Default is 0.
    num_iterations : int, optional
        Number of iterations for PSO optimization. Default is 3.
    maxiter : int, optional
        Maximum number of iterations for each PSO run. Default is 1000.
    swarmsize : int, optional
        Number of particles in the PSO swarm. Default is 100.
    units : str, optional
        Units of the observed data. Default is "keV".
    method_ : str, optional
        Fitting method to use, either "discrete" or "extended". Default is "extended".
    extended_binsize : float, optional
        Bin size for the extended method. Default is 0.01.

    Notes
    -----
    A form will appear to input the necessary bounds for the orbital parameters. These parameters are saved in
    a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to
    re-enter parameters, allowing modification of only those that require adjustment.

    Returns
    -------
    DataFrame : pd.DataFrame
        Contains the best-fit parameters and their standard deviations.
    ph : array-like
        Array of phases corresponding to the predicted data.
    predicted_data : array-like
        Predicted data based on the best-fit parameters.
    chi_squared : float
        Chi-squared statistic weighted by the error, indicating the quality of the fit.
    """


    #............................................Data prep
    parameter_names = ["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis" ,"inclination", "Rstar", "Mstar1", "Mstar2","wind_vel" ,"feature"]
    #............................................
    t = x_data
    x_data, y_err_weight = _define_x_y_sy(x_data,y_data, y_err)

    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(t, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        print("The number of time points does not allow an extended approach. Changing to discrete")
        method_ = "discrete"
        
    #............................................Objective function
    def objective_function(params):

        iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, wind_vel,feature = params
        predicted_data = conic_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity,
                                                     periapsis ,inclination, Rstar, Mstar1, Mstar2, wind_vel, feature, units, method_,extended_binsize)

        chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)

        return chi_squared

    #............................................PS implementation
    lower_bounds, upper_bounds = _manage_bounds(parameter_names, "bounds_orbit", load_directly=load_directly, bound_list=bound_list)


    best_params_list = []
    chi_list = []
    
    for i in range(num_iterations):
        
        best_params, _ = pso(objective_function, lb=lower_bounds, ub=upper_bounds, maxiter = maxiter, swarmsize = swarmsize)

        best_params_list.append(best_params)
        predicted_data = conic_orbit(x_data, *best_params, units,method_,extended_binsize)

        chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)
        
        chi_list.append(chi_squared)
        
    #.............................Collect results
    mean_params = np.mean (best_params_list, axis=0)
    std_params = np.std (best_params_list, axis=0)

    best_iteration = np.argmin(chi_list)
    best_params = best_params_list[best_iteration]
    best_chi = chi_list[best_iteration]

    (iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, wind_vel, feature) = best_params
    
    #.............................Evaluate results
    predicted_data = conic_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, wind_vel, feature, units, method_,extended_binsize)

    chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)
    #.............................Prepare output
    results = []
    
    for param_name, best_param, std_param in zip(parameter_names, best_params, std_params):
        results.append([best_param, std_param])

    df_results = pd.DataFrame(results, index=parameter_names, columns=['Value', 'Std'])
    df_results.index.name = 'Name of the parameter'
    df_results_transposed = df_results.T

    if method_=="extended":

        t =  np.array([np.mean(i) for i in x_data])
        
        ph,_,_ = _orbital_time_to_phase(t ,
        df_results_transposed.iphase.Value,
        df_results_transposed.semimajor.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
    if method_=="discrete":

        t = x_data
        
        ph,_,_ = _orbital_time_to_phase(t ,
        df_results_transposed.iphase.Value,
        df_results_transposed.semimajor.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
    return df_results_transposed, ph, predicted_data, chi_squared
    
# LS FIT------------------------------------------------------------------------------------------------------
def fit_orbit_ls(x_data, y_data, y_err=0, units="keV", method_="extended", extended_binsize=0.01,load_directly=False, bound_list=None):
    """
    Fits orbital modulation data by estimating parameters such as phase, semi-major axis, orbital period,
    eccentricity, inclination, and periapsis using a traditional Least Squares (LS) method.

    The LS method fits the observed data by minimizing the squared differences between observed and predicted values.

    The function supports two fitting methods:
    
    - "discrete": Suitable for discrete data points (e.g., spectra with small orbital phase ranges; faster).
    - "extended": Suitable for data with varying or extended bin sizes, typical with current instrument
      resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Parameters
    ----------
    x_data : array-like
        Time bins of the observed data.
    y_data : array-like
        Observed data points corresponding to each time bin.
    y_err : array-like, optional
        Error associated with each observed data point. Default is 0.
    units : str, optional
        Units of the observed data. Default is "keV".
    method_ : str, optional
        Fitting method to use, either "discrete" or "extended". Default is "extended".
    extended_binsize : float, optional
        Bin size for the extended method. Default is 0.01.

    Notes
    -----
    A form will appear to input the necessary bounds for the orbital parameters. These parameters are saved in
    a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to
    re-enter parameters, allowing modification of only those that require adjustment.

    Returns
    -------
    DataFrame : pd.DataFrame
        Contains the best-fit parameters and their errors.
    ph : array-like
        Array of phases corresponding to the predicted data.
    predicted_data : array-like
        Predicted data based on the best-fit parameters.
    chi_squared : float
        Chi-squared statistic weighted by the error, indicating the quality of the fit.
    """

    #...........................................data prep
    parameter_names = ["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis" ,"inclination", "Rstar", "Mstar1", "Mstar2","wind_vel" ,"feature"]
    
    t=x_data
    
    x_data, y_err_weight = _define_x_y_sy(x_data,y_data, y_err)

    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(t, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        print("The number of time points does not allow an extended approach. Changing to discrete")
        method_ = "discrete"
        
    #............................................LS implementation
    lower_bounds, upper_bounds = _manage_bounds(parameter_names, "bounds_orbit",load_directly=load_directly, bound_list=bound_list)
    
    model_func = lambda x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2,wind_vel, feature: conic_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, wind_vel,feature, units=units, method_=method_,extended_binsize=extended_binsize)
    
    try:
        fit_params, fit_covariance = curve_fit(model_func, x_data, y_data, bounds=[lower_bounds, upper_bounds], maxfev=1000)
    except RuntimeError:
        raise RuntimeError("Curve fitting did not converge. Try adjusting the bounds.")
        
    errors = np.sqrt(np.diag(fit_covariance))
    
    #............................................Evaluate results
    predicted_data = model_func(x_data, *fit_params)
    residuals = y_data - predicted_data

    rss = np.sum(residuals**2)
    tss = np.sum((y_data - np.mean(y_data))**2)
    r_squared = 1 - (rss / tss)
    
    #..............................................Prepare output

    results = []
    for param_name, best_param, std_param in zip(parameter_names, fit_params, errors):
        results.append([best_param, std_param])
        
    df_results = pd.DataFrame(results, index=parameter_names, columns=['Value', 'Std'])
    df_results.index.name = 'Name of the parameter'
    df_results_transposed = df_results.T

    
    if method_=="extended":

        t =  np.array([np.mean(i) for i in x_data])
        
        ph,_,_ = _orbital_time_to_phase(t ,
        df_results_transposed.iphase.Value,
        df_results_transposed.semimajor.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
    if method_=="discrete":

        t = x_data
        
        ph,_,_ = _orbital_time_to_phase(t ,
        df_results_transposed.iphase.Value,
        df_results_transposed.semimajor.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
    return df_results_transposed,  ph, predicted_data, r_squared
