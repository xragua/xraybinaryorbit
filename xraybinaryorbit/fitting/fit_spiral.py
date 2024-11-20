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

# SPIRAL #############################################################################
def spiral(x_data, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units, method_, extended_binsize):
    """
    Simulate the Doppler shift for a spiral structure in orbit around a central object. This private function is called
    by fit_spiral_ps and fit_spiral_ls.

    This function calculates the Doppler shift in a spiral orbit, where the radial distance of the orbiting material
    expands exponentially as a function of phase. The Doppler shift is computed based on the motion of the material
    and can be returned in different units, such as keV, seconds, or wavelength (angstroms). The function allows
    calculations either over discrete time points or extended time bins, based on the `method_` parameter.

    Parameters
    ----------
    x_data : array-like
        The input time data (or time bins) to compute the Doppler shift values.
    iphase_spiral : float
        The initial phase of the spiral orbit, typically as a fraction of the orbital period or in radians.
    semimajor_spiral : float
        The semi-major axis of the spiral orbit, in solar radii or other relevant units.
    b : float
        The exponential growth factor of the spiral, determining how the radius expands with phase.
    omega : float
        The angular velocity of the spiral, typically in radians per second.
    inclination_spiral : float
        The inclination of the spiral in degrees, relative to the plane of the sky.
    feature : float
        The spectral feature of interest, either in keV or another unit, that will be Doppler-shifted.
    units : str
        The units for the feature, such as "keV", "s" (seconds), or "amstrong" (angstroms). Determines how the feature is Doppler shifted.
    method_ : str
        The method used for calculation. Can be "extended" (for time bins) or "discrete" (for individual time points).
    extended_binsize : float
        The bin size for the extended method. This is used to average Doppler shifts and velocities over each time bin.

    Returns
    -------
    equation : array-like
        The Doppler-shifted values of the feature, depending on the `units` parameter.
        This can be in keV, seconds, or wavelength (angstroms), depending on the chosen units.

    Notes
    -----
    - The function simulates a spiral structure, where the radius grows exponentially as a function of the phase, defined by
      the exponential growth factor `b` and the angular velocity `omega`.
    - For the "extended" method, the function computes average Doppler shift values over time bins, while for the "discrete" method,
      it computes Doppler shifts at individual time points.
    - The Doppler shift is calculated based on the radial velocity of the orbiting material, and the output is returned in the specified units.
    - The Doppler shift accounts for both the motion of the spiral and its inclination with respect to the observer.

    """

    #...................................................
    feature_ = {
        "keV": kev_ams/feature,
        "s": feature,  # No conversion needed for seconds
        "amstrong": feature  # No conversion needed for lambda
    }
    
    feature = feature_.get(units, 1)
    t = x_data
    #...................................................
    
    if method_=="extended":
    
        phbins =  (t - min(t[0])) * omega + iphase_spiral
        size_phase_bin = np.diff(phbins)
        minsizebin = min(size_phase_bin)
        maxph = max(phbins[-1])
        
        phase = np.arange(0,maxph,max(minsizebin/10,maxph/100000))
        
        vdop_bin=[]
        
        R = semimajor_spiral * np.exp(b * 2 * np.pi * phase)
        vdop = -R * rsun_m * omega * np.sin(2 * np.pi * phase ) * np.sin(2 * np.pi * inclination_spiral/360)
        
        
        for i in range(len(phbins)):
        
            if ( size_phase_bin[i] >= extended_binsize):

                vdop_bin.append(np.mean(vdop[(phase >= phbins[i,0]) & (phase <= phbins[i,1])]))
                
            if ( size_phase_bin[i] < extended_binsize):
            
                phase_puntual =  (np.mean(t[i])-min(t[0]))* omega + iphase_spiral
                R_puntual = semimajor_spiral * np.exp(b * 2 * np.pi * phase_puntual)
                vdop_bin.append( -R_puntual * rsun_m * omega * np.sin(2 * np.pi * phase_puntual ) * np.sin(2 * np.pi * inclination_spiral/360))
            
        
    if method_=="discrete":
    #...................................................
        phase_discrete =  (t - min(t)) * omega + iphase_spiral
        R = semimajor_spiral * np.exp(b * 2 * np.pi * phase_discrete)

        vdop_bin = -R * rsun_m * omega * np.sin(2 * np.pi * phase_discrete ) * np.sin(2 * np.pi * inclination_spiral/360)
        #..................................................
    vdop =  np.array(vdop_bin)

    equation_ = {
        "keV": kev_ams/(feature * (c + vdop) / c),
        "s": (feature * (c + vdop) / c),  # No conversion needed for seconds
        "amstrong": (feature * (c + vdop) / c)  # No conversion needed for lambda
    }
    
    # Apply the conversion factor based on the specified units
    equation = equation_.get(units, 1)
    
    return equation
    

# PS FIT------------------------------------------------------------------------------------------------------
def fit_spiral_ps(x_data, y_data, y_err=0, num_iterations=3, maxiter=1000, swarmsize=100,
                  units="keV", method_="extended", extended_binsize=0.01,load_directly=False, bound_list=None):
    """
    Fits orbital modulation data by estimating parameters for a spiral orbit using particle swarm optimization (PSO).

    The fitting process uses a PSO algorithm, which iteratively improves parameter estimates by minimizing the
    chi-squared difference between the observed and predicted data.

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
        Number of PSO iterations to perform. Default is 3.
    maxiter : int, optional
        Maximum number of iterations for PSO. Default is 1000.
    swarmsize : int, optional
        Number of particles in the swarm for PSO. Default is 100.
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
        Chi-squared statistic weighted by the error, indicating the quality of fit.
    """


    #............................................data prep
    parameter_names = ["iphase_spiral", "semimajor_spiral", "b", "omega", "inclination_spiral", "feature"]
    #............................................
    x_data, y_err_weight = _define_x_y_sy(x_data,y_data, y_err)
    
    t=x_data
        
    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(t, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        print("The number of time points does not allow an extended approach. Changing to discrete")
        method_ = "discrete"
        
    #............................................Objective function
    def objective_function(params):
    
        iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature = params

        predicted_data = spiral(x_data, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units, method_,extended_binsize)
      
        chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)

        return chi_squared
        
    #............................................ PS implementation
    lower_bounds, upper_bounds = _manage_bounds(parameter_names, "bounds_spiral",load_directly=load_directly, bound_list=bound_list)
    best_params_list = []
    chi_list = []
    
    for i in range(num_iterations):
        
        best_params, _ = pso(objective_function, lb=lower_bounds, ub=upper_bounds, maxiter = maxiter, swarmsize = swarmsize)
        best_params_list.append(best_params)
        predicted_data = spiral(x_data, *best_params, units,method_,extended_binsize)

        chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)
        
        chi_list.append(chi_squared)
        

    #............................. Collect results
    mean_params = np.mean (best_params_list, axis=0)
    std_params = np.std (best_params_list, axis=0)

    best_iteration = np.argmin(chi_list)
    best_params = best_params_list[best_iteration]
    best_chi = chi_list[best_iteration]

    (iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature) = best_params


    #............................. Evaluate results
    predicted_data = spiral(x_data, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units, method_,extended_binsize)
    chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)
    #............................. Prepare output

    results = []
    for param_name, best_param, std_param in zip(parameter_names, best_params, std_params):
        results.append([best_param, std_param])

    df_results = pd.DataFrame(results, index=parameter_names, columns=['Value', 'Std'])
    df_results.index.name = 'Name of the parameter'

    df_results_transposed = df_results.T

    ph=[]
    
    if method_=="extended":
        t0 = min(x_data[0])
        center_t =  np.array([np.mean(i) for i in x_data])
        ph =  (center_t - t0) * df_results_transposed.omega.Value + df_results_transposed.iphase_spiral.Value
        
    if method_=="discrete":
        t0 = min(x_data)
        t=x_data
        ph = ph =  (t - t0) * df_results_transposed.omega.Value + df_results_transposed.iphase_spiral.Value
    
        
    return df_results_transposed, ph, predicted_data, chi_squared

    
# LS FIT------------------------------------------------------------------------------------------------------
def fit_spiral_ls(x_data, y_data, y_err=0, units="keV", method_="extended", extended_binsize=0.01,load_directly=False, bound_list=None):
    """
    Fits orbital modulation data by estimating parameters for a spiral orbit using a traditional least squares (LS) method.
    Although provided for completeness, the LS method may have limitations due to the complexity of the model.

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
        Chi-squared statistic weighted by the error, indicating the quality of fit.
    """



    parameter_names = ["iphase_spiral", "semimajor_spiral", "b", "omega", "inclination_spiral", "feature"]

    x_data, y_err_weight = _define_x_y_sy(x_data,y_data, y_err)
   
    t=x_data
        
    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(t, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        print("The number of time points does not allow an extended approach. Changing to discrete")
        method_ = "discrete"
    #............................................LS implementation
    lower_bounds, upper_bounds = _manage_bounds(parameter_names, "bounds_spiral",load_directly=load_directly, bound_list=bound_list)
    
    model_func = lambda x_data, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature: spiral(x_data,  iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units=units, method_=method_,extended_binsize=extended_binsize)
    
    try:
        fit_params, fit_covariance = curve_fit(model_func, x_data, np.array(y_data), bounds=[lower_bounds, upper_bounds], maxfev=1000)
    except RuntimeError:
        raise RuntimeError("Curve fitting did not converge. Try adjusting the bounds.")

    errors = np.sqrt(np.diag(fit_covariance))
    
    #............................................. Evaluate results
    predicted_data = model_func(x_data, *fit_params)

    residuals = y_data - predicted_data
    rss = np.sum(residuals**2)

    tss = np.sum((y_data - np.mean(y_data))**2)

    r_squared = 1 - (rss / tss)

    print(f"Coefficient of Determination (R-squared): {r_squared:.4f}")
    
    #.............................Prepare output
    results = []
    for param_name, best_param, std_param in zip(parameter_names, fit_params, errors):
        results.append([best_param, std_param])

    df_results = pd.DataFrame(results, index=parameter_names, columns=['Value', 'Std'])
    df_results.index.name = 'Name of the parameter'

    df_results_transposed = df_results.T
    print(df_results_transposed.to_string())
    
    ph=[]
    
    if method_=="extended":
        t0 = min(x_data[0])
        center_t =  np.array([np.mean(i) for i in x_data])
        ph =  (center_t - t0) * df_results_transposed.omega.Value + df_results_transposed.iphase_spiral.Value
        
    if method_=="discrete":
        t0 = min(x_data)
        t=x_data
        ph = ph =  (t - t0) * df_results_transposed.omega.Value + df_results_transposed.iphase_spiral.Value
    
    return df_results_transposed,  ph, predicted_data, r_squared
