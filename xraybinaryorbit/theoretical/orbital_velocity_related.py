import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
import math
import tkinter as tk
from tkinter import messagebox

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
def orbital_phase_to_time(ph, precision=0.01,load_directly=False, parameter_list=None):
    """
    Converts an orbital phase array to a time array for a compact object orbiting a companion star.
    The compact object moves faster at periastron than at apoastro due to the conservation of angular momentum
    and Kepler's laws of planetary motion.

    The increased orbital speed at periastron occurs because, as the compact object moves closer to the central star,
    it must travel faster to maintain the total angular momentum of the system. This behavior follows Kepler's second
    law, which states that objects sweep out equal areas in equal times, and the gravitational force strengthens
    as the bodies approach and weakens as they move apart.

    Parameters
    ----------
    ph : array-like
        Orbital phase array.
    precision : float, optional
        Resolution for the phase array. Default is 0.01.

    Notes
    -----
    A form will appear to input the necessary orbital parameters. These parameters are saved in a .txt file
    in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter parameters,
    allowing modification of only those that require adjustment.
    By setting load_directly=True the data will be authomatically loaded into the function without rising a form. By providing parameter_list=(list of parameters) the parameters will be trated as an imput and saved within the current directory for subsequent runs.

    Returns
    -------
    ph : array-like
        The input orbital phase array.
    time : array-like
        Time array corresponding to the orbital phase.
    W : array-like
        Angular velocity array corresponding to the orbital phase.
    """


    #.............................Load parameters
    parameter_names = ["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis", "Rstar", "Mstar1", "Mstar2"]
    fixed_values = _manage_parameters(parameter_names, "phase_time",load_directly=load_directly,parameter_list=parameter_list )
    
    iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2 = fixed_values
    
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
    
    
def orbital_time_to_phase(t, precision=0.01,load_directly=False, parameter_list=None):
    """
    Converts an orbital time array to a phase array for a compact object orbiting a companion star.
    The compact object moves faster at periastron than at apoastro due to the conservation of angular momentum
    and Kepler's laws of planetary motion.

    The increased orbital speed at periastron occurs because, as the compact object moves closer to the central star,
    it must travel faster to maintain the total angular momentum of the system. This behavior follows Kepler's second
    law, which states that objects sweep out equal areas in equal times, and the gravitational force strengthens
    as the bodies approach and weakens as they move apart.

    Parameters
    ----------
    t : array-like
        Orbital time array.
    precision : float, optional
        Resolution for the phase array. Default is 0.01.

    Notes
    -----
    A form will appear to input the necessary orbital parameters. These parameters are saved in a .txt file
    in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter parameters,
    allowing modification of only those that require adjustment.
    By setting load_directly=True the data will be authomatically loaded into the function without rising a form. By providing parameter_list=(list of parameters) the parameters will be trated as an imput and saved within the current directory for subsequent runs.

    Returns
    -------
    ph : array-like
        Orbital phase array corresponding to the input time array.
    time : array-like
        Time array corresponding to the orbital phase.
    W : array-like
        Angular velocity array corresponding to the orbital phase.

    """


    #.............................Load parameters
    parameter_names = ["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis", "Rstar", "Mstar1", "Mstar2"]
    fixed_values = _manage_parameters(parameter_names, "time_phase",load_directly=load_directly,parameter_list=parameter_list )
    
    iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2 = fixed_values
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
##################################################################################
