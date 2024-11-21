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

def nh_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, Mass_loss_rate, wind_infinite_velocity, beta, method_, extended_binsize):
    """
    Calculate the hydrogen column density (NH) in an orbital system with a stellar wind. This private function is called by
    fit_nh_ps.

    This function calculates the hydrogen column density (NH) in a binary orbital system with a star that has a stellar wind.
    The function integrates the mass density along the line of sight as the system moves through different orbital phases,
    taking into account the orbital motion and wind velocity profile. It supports both "extended" and "discrete" methods for
    calculation, where the "extended" method computes averaged values over time bins, and the "discrete" method computes
    NH at individual time points.

    Parameters
    ----------
    x_data : array-like
        The input time data (or time bins) to compute NH values.
    iphase : float
        The initial phase of the orbit, typically as a fraction of the orbital period or in radians.
    semimajor : float
        The semi-major axis of the orbit, in solar radii.
    orbitalperiod : float
        The orbital period of the system in days.
    eccentricity : float
        The eccentricity of the orbit, ranging from 0 (circular) to 1 (parabolic).
    periapsis : float
        The argument of periapsis of the orbit, describing the orientation of the orbit in the plane of motion.
    inclination : float
        The inclination of the orbit in degrees, relative to the plane of the sky.
    Rstar : float
        The radius of the star, typically in solar radii.
    Mstar1 : float
        The mass of the primary star, typically in solar masses.
    Mstar2 : float
        The mass of the secondary object, typically in solar masses.
    Mass_loss_rate : float
        The mass loss rate from the star due to stellar wind, typically in solar masses per year.
    wind_infinite_velocity : float
        The terminal wind velocity of the stellar wind, typically in km/s.
    beta : float
        The velocity law exponent for the wind velocity profile (used in the beta-law to describe wind acceleration).
    method_ : str
        The method used for calculation. Can be "extended" (for time bins) or "discrete" (for individual time points).
    extended_binsize : float
        The bin size for the extended method. This is used to average NH values over each time bin.

    Returns
    -------
    nh_bin : np.ndarray
        The hydrogen column density (NH) values for each time bin or each time point, depending on the calculation method.
        The output array is in units of 10^22 cm^-2.

    Notes
    -----
    - The function computes NH based on the integration of mass density along the line of sight, taking into account
      the orbital motion and wind profile. The wind velocity follows a beta-law, where the wind velocity increases as a
      function of distance from the star.
    - For the "extended" method, NH is computed over time bins, and the values are averaged for each bin. For the "discrete"
      method, NH is computed at individual time points.
    - The mass density is calculated using the mass loss rate and wind velocity, and NH is obtained by integrating the density
      along the line of sight.
    - The final NH values are returned in units of 10^22 cm^-2.

    """
    


    def nh_calc(th):

        if not isinstance(th, (list, np.ndarray)):
            th = [th]

        nh = []

        for phase in th:
            Rorb = (abar * (1 - eccentricity**2) / (1 + eccentricity * np.cos((phase - periapsis / 360) * 2 * np.pi))) * Rstar_cm  # in cm

            def integrand(z):
            
                alpha = np.arccos(np.cos(phase * 2 * np.pi) * np.sin(inclination * 2 * np.pi / 360))
                x = np.sqrt(Rorb**2 + z**2 - 2 * Rorb * z * np.cos(alpha))
                v = vinf_cm_s * (1 - Rstar_cm / x)**beta
                rho = M_dot_grams / (4 * np.pi * v * x**2)
                return rho

            ne, _ = quad(integrand, 1, Rstar_cm * 1000)
            nh.append(ne * na / 1e22)

        return np.nan_to_num(nh)


    # Constants and conversions
    vinf_cm_s = wind_infinite_velocity * 100000  # Terminal velocity in cm/s
    Rstar_cm = Rstar * rsun_cm  # Star radius in cm
    
    abar = semimajor * max(Mstar1,Mstar2)/(Mstar1+Mstar2)
    M_dot_grams = Mass_loss_rate * msun/(365*24*60*60) #Mass_loss_rate gr/s
    Rstar_cm = Rstar*rsun_cm #R* in cm

    t = x_data
    shape_t = t.shape
    t_to_phase = t.reshape(-1)

    if method_ == "extended":
        
        ph_from_t, _, _ = _orbital_time_to_phase(t_to_phase, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)
        t_to_phase_punctual = np.mean(t, axis=1)
        phase_punctual, _, _ = _orbital_time_to_phase(t_to_phase_punctual, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)

        phbins = ph_from_t.reshape(shape_t)
        size_phase_bin = np.diff(phbins)
        minsizebin = min(size_phase_bin)
        maxph = max(phbins[-1])
        tphase = np.arange(0, maxph, minsizebin / 3)

        nh_for_bin = nh_calc(tphase)
        nh_bin = []

        for i in range(len(phbins)):
        
            if size_phase_bin[i] >= extended_binsize:
            
                nh_bin.append(np.mean(nh_for_bin[(tphase >= phbins[i, 0]) & (tphase <= phbins[i, 1])]))
                
            else:
                nh_bin.append(nh_calc(np.array([phase_punctual[i]]))[0])

    elif method_ == "discrete":

        
        phase_discrete, _, _ = _orbital_time_to_phase(t, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)
        
        nh_bin = nh_calc(phase_discrete)

    return np.array(nh_bin,dtype=np.float64)
    
# PS FIT------------------------------------------------------------------------------------------------------
def fit_nh_ps(x_data, y_data, y_err=0, num_iterations=3, maxiter=200, swarmsize=20,
              method_="extended", extended_binsize=0.01,load_directly=False, bound_list=None):
    """
    Fits the column density (NH1, x 10^22 cm^-2) encountered by radiation emitted at each orbital phase as
    it travels towards an observer, assuming a spherically distributed, neutral (unionized) stellar wind based
    on the CAK model.

    The fitting process uses a particle swarm optimization (PSO) algorithm to minimize the chi-squared
    difference between observed and predicted data points.

    The function supports two fitting methods:
    
    - "discrete": Suitable for discrete data points (e.g., spectra with small orbital phase ranges; faster).
    - "extended": Suitable for data with varying or extended bin sizes, typical for current instrument
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
        Number of iterations for the PSO algorithm. Default is 3.
    maxiter : int, optional
        Maximum number of iterations for each PSO run. Default is 200.
    swarmsize : int, optional
        Number of particles in the PSO swarm. Default is 20.
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
    df_results_transposed : pd.DataFrame
        Transposed DataFrame of the best-fit parameters and their standard deviations.
    ph : array-like
        Array of phases corresponding to the predicted data.
    predicted_data : array-like
        Predicted data based on the best-fit parameters.
    chi_squared : float
        Chi-squared statistic weighted by the error, indicating the quality of fit.
    """


#............................................Data prep.
    parameter_names=["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis" ,"inclination", "Rstar", "Mstar1", "Mstar2", "Mdot", "v_inf", "beta"]
    
    t = x_data
    x_data, y_err_weight = _define_x_y_sy(x_data,y_data, y_err)
        
    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(t, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        print("The number of time points does not allow an extended approach. Changing to discrete")
        method_ = "discrete"
#............................................Objective function

    def objective_function(params):

        iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, Mdot, v_inf, beta = params
        predicted_data = nh_orbit(x_data,  iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1,
                                     Mstar2, Mdot, v_inf, beta, method_,extended_binsize)

        chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)
         
        return chi_squared
#............................................PS implementation
    lower_bounds, upper_bounds = _manage_bounds(parameter_names, "bounds_nh",load_directly=load_directly, bound_list=bound_list)
    best_params_list = []
    chi_list = []

    for i in range(num_iterations):
    

        best_params, _ = pso(objective_function, lb=lower_bounds, ub=upper_bounds, maxiter = maxiter, swarmsize = swarmsize, phig=2)

        best_params_list.append(best_params)
        predicted_data = nh_orbit(x_data, *best_params, method_,extended_binsize)

        chi_squared = _chi_squared_weighted(y_data, y_err_weight, predicted_data)

        chi_list.append(chi_squared)
        #print(chi_squared)
        
#.............................Get the results
    mean_params = np.mean (best_params_list, axis=0)
    std_params = np.std (best_params_list, axis=0)

    best_iteration = np.argmin(chi_list)
    best_params = best_params_list[best_iteration]
    best_chi = chi_list[best_iteration]

    ( iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, Mdot, v_inf, beta) = best_params

    predicted_data = nh_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, Mdot, v_inf, beta, method_,extended_binsize)
#............................. Final evaluation
    chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)
#............................. Prepare output

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
