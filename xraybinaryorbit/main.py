#PKGS
##########################################################################################
##########################################################################################
##########################################################################################import numpy as np
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

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="divide by zero encountered in divide")

#Constants
##########################################################################################
c = 299792458

msun = (1.98847*10**30)*1000 #gr
rsun_m = 696340*1000 #
rsun_cm = 696340*1000*100 #cm

kev_ams = 1.23984193

na = 6.02214076*10**23/1.00797
mu = 0.5
mp = 1.67E-24


#UNITS MESSAGE
##########################################################################################
print("""

HELLO, nice to see you! :)

PLEASE READ THIS, IT'S VERY IMPORTANT:

These are the units that must be used within this package:

- Rstar: Solar radius
- Mstar: Solar masses
- Inclination: Sexagesimal degrees
- Periapsis: Sexagesimal degrees
- Semimajor: Stellar radius
- Periods: Days (Periods in the case of the period_sliding_window function will support any units)
- Iphase: Radians

A list of the functions contained in this package will be displayed by runing the function list_functions().

As these functions use a lot of parameters, which can sometimes be difficult to handle, we have implemented a user-friendly method for parameter input:
A form will be displayed, and the parameters will be saved in the directory for further interactions. These saved parameters will be used if new parameters are not provided.
For the function to work, the submit button must be pressed.
If the parameters are already saved within the working directory, setting "load_directly=True" no form will be displayed and that parameters will be used within the function.
Alternatively, the input parameters or bounds can be provided as lists, by providing a "parameter_list" or "bound_list" as imputs.

Please, take into account that fits in general will take A LOT of time to complete.

If you need help, contact graciela.sanjurjo@ua.es.
""")

#Constants
##########################################################################################
c = 299792458

msun = (1.98847*10**30)*1000 #gr
rsun_m = 696340*1000 #
rsun_cm = 696340*1000*100 #cm

kev_ams = 1.23984193

na = 6.02214076*10**23/1.00797
mu = 0.5
mp = 1.67E-24

#Helper functions
##########################################################################################
def _advise():
    print("This function may take a significant amount of time to complete. Please consider starting with a smaller swarm size (swarmsize) and fewer iterations (maxiter). You can gradually increase these parameters in subsequent runs if the results are not satisfactory enough.")


def list_functions():
    print("""
    Theoretical:
    - doppler_orbit_theoretical
    - doppler_spiral_theoretical
    - doppler_disc_theoretical
    - doppler_spiral_in_orbit_theoretical
    - density_through_orbit_theoretical
    - absorption_column_through_orbit_theoretical
    - ionization_map_phase
    - orbital_phase_to_time
    - orbital_time_to_phase
    
    Fitting:
    (ps is for particle swarm and ls for least squares. Please, note that ps is preferred as ls does not allways converge)
    - fit_orbit_ps
    - fit_orbit_ls
    - fit_disc_ps
    - fit_disc_ls
    - fit_spiral_ps
    - fit_spiral_ls
    - fit_spiral_in_orbit_ps
    - fit_spiral_in_orbit_ls
    - nh_orbit
    - fit_nh_ps
    
    Timming:
    - hr
    - cr
    - rebin_snr
    - rebin_bins
    - fold_pulse
    - period_sliding_window
    """)
