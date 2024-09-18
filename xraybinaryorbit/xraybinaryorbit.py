#PKGS
##########################################################################################
##########################################################################################
##########################################################################################
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

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="divide by zero encountered in divide")

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

#Helper functions
##########################################################################################
def advise():
    print("Nice! we are working on getting your parameters. It will take a while... please, grab a cup of coffe, go for a walk, go on some vacations... we might finish by Chistmast, unless is Christmast now, if thats the case, by Easter.")



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
##########################################################################################
def gaussian(x, mean, sigma):
    return np.exp(-(x - mean)**2 / (2 * sigma**2)) / (sigma * np.sqrt(2 * np.pi))
#...................................................  Prepare timebins for fitting functions
def time_pairs(edges):

    edge_pairs = [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]
    pairs_=[]

    for pair in edge_pairs:
        pairs_.append(pair)

    return pairs_
#...................................................PCHIP (Piecewise Cubic Hermite Interpolating Polynomial) is a method used for interpolation, particularly for interpolating data points with unevenly spaced independent variable (x) values.
def interpolate_pchip(tmin,tmax,tinterval,t,y, sy):
    
    t_new = np.arange(tmin,tmax,tinterval)
    
    cs = PchipInterpolator(t, y)
    scs = PchipInterpolator(t, sy)
    new_data = cs(t_new)
    snew_data = scs(t_new)
    
    return t_new, new_data, snew_data

#.......................................................
#..................................................... Properly prepare imput for fitting functions
def define_x_y_sy(x_data, y_data, y_err):
    
    x_data = np.array(x_data)
    y_data = np.array(y_data)
    y_err = np.array(y_err)
    
    if np.sum(y_err) !=0 :
    
        if len(y_err.shape) == 2:
        
            y_err_weight = np.array(y_err[1] + y_err[0])
            
        elif len(y_err.shape) == 1:
        
            y_err_weight = np.array(y_err)
            
    if np.sum(y_err) ==0 :
        y_err_weight = np.ones(len(y_data))
            
    if len(y_data) + 1 == len(x_data) and len(x_data.shape) == 1:
        x_data = time_pairs(x_data)
        
    return x_data, y_err_weight
#................................................... Weighted to the error chi squared
def chi_squared_weighted(y_data, y_err, y_pred):

    y_data = np.asarray(y_data)
    y_err = np.asarray(y_err)
    y_pred = np.asarray(y_pred)

    chi_squared = np.sum(((y_data - y_pred) / y_err) ** 2)
    return chi_squared
    
def chi_squared(y_data, y_pred):

    y_data = np.asarray(y_data)
    y_err = np.asarray(y_err)
    y_pred = np.asarray(y_pred)

    chi_squared = np.sum((y_data - y_pred) ** 2)
    return chi_squared

#................................................... Easy way to manage values
def copy_fields(source_file, destination_file, params2=None):
    # Read the content of the source file
    with open(source_file, "r") as source:
        source_lines = source.readlines()
        source_param_names = source_lines[0].strip().split(",")
        source_fixed_values = [float(val) if val != 'nan' else np.nan for val in source_lines[1].strip().split(",")]

    # Check if destination file exists
    try:
        with open(destination_file, "r") as destination:
            dest_lines = destination.readlines()
            dest_param_names = dest_lines[0].strip().split(",")
            dest_fixed_values = [float(val) if val != 'nan' else np.nan for val in dest_lines[1].strip().split(",")]
    except FileNotFoundError:
        # If destination file doesn't exist, create it with params2 and load machine parameters from source_file
        if params2:
            dest_param_names = params2
            dest_fixed_values = [np.nan] * len(params2)
        else:
            raise ValueError("Destination file doesn't exist and params2 is not provided.")

        # Load machine parameters from source_file
        machine_params = []
        with open(source_file, "r") as source:
            for line in source:
                if line.strip():  # Skip empty lines
                    machine_params.append(line.strip())

        # Update dest_fixed_values with machine parameters
        dest_fixed_values[:len(machine_params)] = machine_params

    # Find identical fields and update values
    for i, param_name in enumerate(source_param_names):
        if param_name in dest_param_names:
            index = dest_param_names.index(param_name)
            dest_fixed_values[index] = source_fixed_values[i]

    # Write the updated values to the destination file
    with open(destination_file, "w") as destination:
        destination.write(",".join(map(str, dest_param_names)) + "\n")
        destination.write(",".join(map(str, dest_fixed_values)) + "\n")

    
def load_values_to_interface(param_list, fixed_values, name):
    global fixed_values_list
    fixed_values_list = fixed_values

    root = tk.Tk()
    root.title("Values Input Form")

    # Create entry fields for fixed values
    fixed_entries = {}
    for i, param in enumerate(param_list):
        tk.Label(root, text=f"{param}").grid(row=i, column=0)
        fixed_entries[param] = tk.Entry(root)
        fixed_entries[param].insert(0, str(fixed_values_list[i]))  # Insert previous value
        fixed_entries[param].grid(row=i, column=1)

    # Create submit button
    def submit_form_values():
        global fixed_values_list
        new_values = []
        for param in param_list:
            value = float(fixed_entries[param].get())
            if param == "eccentricity" and not (0 <= value <= 0.999999):
                messagebox.showerror("Invalid Input", "Eccentricity must be between 0 and 0.999999")
                return  # Stop submission if validation fails
            if param == "semimajor" and not (0 < value ):
                messagebox.showerror("Invalid Input", "Semimajor should be higher than 0")
                return  # Stop submission if validation fails
            if param == "Mstar1" and not (0 < value ):
                messagebox.showerror("Invalid Input", "Mstar1 should be higher than 0")
                return  # Stop submission if validation fails
            if param == "Mstar2" and not (0 < value ):
                messagebox.showerror("Invalid Input", "Mstar1 should be higher than 0")
                return  # Stop submission if validation fails
            if param == "Rstar" and not (0 < value ):
                messagebox.showerror("Invalid Input", "Rstar should be higher than 0")
                return  # Stop submission if validation fails
            new_values.append(value)
        
        fixed_values_list = new_values
        root.withdraw()
        root.destroy()  # Close the root window after submitting the form

    submit_button = tk.Button(root, text="Submit", command=submit_form_values)
    submit_button.grid(row=len(param_list), padx=10, pady=10)

    root.mainloop()
    

def manage_parameters(param_list, name):
    global fixed_values_list

    # Load previous fixed values if available
    try:
        with open(f"{name}.txt", "r") as file:
            lines = file.readlines()
            fixed_values_list = [float(val) for val in lines[1].strip().split(",")]
    except FileNotFoundError:
        fixed_values_list = [np.nan] * len(param_list)

    load_values_to_interface(param_list, fixed_values_list, name)

    # Save fixed values to file
    with open(f"{name}.txt", "w") as file:
        file.write(",".join(map(str, param_list)) + "\n")  # Write parameter names
        file.write(",".join(map(str, fixed_values_list)) + "\n")  # Write fixed values

    return fixed_values_list


#................................................... Easy way to manage bounds in fitting functions

def load_bounds_to_interface(param_list, lower_bounds, upper_bounds, name):
    global lower_bounds_list, upper_bounds_list
    
    lower_bounds_list = lower_bounds
    upper_bounds_list = upper_bounds

    root = tk.Tk()
    root.title("Bounds Input Form")

    # Create entry fields for lower bounds
    lower_entries = {}
    for i, param in enumerate(param_list):
        tk.Label(root, text=f"Lower {param}").grid(row=i, column=0)
        lower_entries[param] = tk.Entry(root)
        lower_entries[param].insert(0, str(lower_bounds_list[i]))  # Insert previous value
        lower_entries[param].grid(row=i, column=1)

    # Create entry fields for upper bounds
    upper_entries = {}
    for i, param in enumerate(param_list):
        tk.Label(root, text=f"Upper {param}").grid(row=i, column=2)
        upper_entries[param] = tk.Entry(root)
        upper_entries[param].insert(0, str(upper_bounds_list[i]))  # Insert previous value
        upper_entries[param].grid(row=i, column=3)

    # Create submit button
    def submit_form_bounds():
        global lower_bounds_list, upper_bounds_list
        lower_bounds_list = [float(lower_entries[param].get()) for param in param_list]
        upper_bounds_list = [float(upper_entries[param].get()) for param in param_list]
        root.withdraw()
        root.destroy() # Close the root window after submitting the form

    submit_button = tk.Button(root, text="Submit", command=submit_form_bounds)
    submit_button.grid(row=len(param_list), columnspan=4)

    root.mainloop()




def load_bounds_to_interface(param_list, lower_bounds, upper_bounds, name):
    global lower_bounds_list, upper_bounds_list
    
    lower_bounds_list = lower_bounds
    upper_bounds_list = upper_bounds

    root = tk.Tk()
    root.title("Bounds Input Form")

    # Create entry fields for lower bounds
    lower_entries = {}
    for i, param in enumerate(param_list):
        tk.Label(root, text=f"Lower {param}").grid(row=i, column=0)
        lower_entries[param] = tk.Entry(root)
        lower_entries[param].insert(0, str(lower_bounds_list[i]))  # Insert previous value
        lower_entries[param].grid(row=i, column=1)

    # Create entry fields for upper bounds
    upper_entries = {}
    for i, param in enumerate(param_list):
        tk.Label(root, text=f"Upper {param}").grid(row=i, column=2)
        upper_entries[param] = tk.Entry(root)
        upper_entries[param].insert(0, str(upper_bounds_list[i]))  # Insert previous value
        upper_entries[param].grid(row=i, column=3)

    # Validation function
    def validate_bounds(param, lower_value, upper_value):
        if param == "eccentricity":
            if not (0 <= lower_value < upper_value < 1):
                messagebox.showerror("Invalid Input", "Eccentricity bounds must be 0 <= lower < upper < 1")
                return False
        elif param in ["semimajor","semimajor2","Mstar2","Mstar1", "Rstar"]:
            if not (0 < lower_value < upper_value):
                messagebox.showerror("Invalid Input", f"{param} bounds must be 0 < lower < upper")
                return False
        else:
            if not (lower_value < upper_value):
                messagebox.showerror("Invalid Input", f"{param} bounds must be lower < upper")
                return False
        return True

    # Create submit button
    def submit_form_bounds():
        global lower_bounds_list, upper_bounds_list
        new_lower_bounds = []
        new_upper_bounds = []
        for param in param_list:
            try:
                lower_value = float(lower_entries[param].get())
                upper_value = float(upper_entries[param].get())
            except ValueError:
                messagebox.showerror("Invalid Input", f"{param} bounds must be valid numbers")
                return  # Stop submission if value conversion fails

            if not validate_bounds(param, lower_value, upper_value):
                return  # Stop submission if validation fails

            new_lower_bounds.append(lower_value)
            new_upper_bounds.append(upper_value)
        
        lower_bounds_list = new_lower_bounds
        upper_bounds_list = new_upper_bounds
        root.withdraw()
        root.destroy()  # Close the root window after submitting the form

    submit_button = tk.Button(root, text="Submit", command=submit_form_bounds)
    submit_button.grid(row=len(param_list), columnspan=4, padx=10, pady=10)

    root.mainloop()
    
#........................................................................................
def manage_bounds(param_list, name):
    global lower_bounds_list, upper_bounds_list

    # Load previous bounds if available
    try:
        with open("bounds_"+str(name)+".txt", "r") as file:
            lines = file.readlines()
            param_names = lines[0].strip().split(",")  # Extract parameter names
            lower_bounds_list = [float(val) for val in lines[1].strip().split(",")]
            upper_bounds_list = [float(val) for val in lines[2].strip().split(",")]
    except FileNotFoundError:
        param_names = param_list
        lower_bounds_list = [np.nan] * len(param_list)
        upper_bounds_list = [np.nan] * len(param_list)

    load_bounds_to_interface(param_names, lower_bounds_list, upper_bounds_list, name)

    # Save bounds to file
    with open("bounds_"+str(name)+".txt", "w") as file:
        file.write(",".join(map(str, param_names)) + "\n")  # Write parameter names
        file.write(",".join(map(str, lower_bounds_list)) + "\n")  # Write lower bounds
        file.write(",".join(map(str, upper_bounds_list)) + "\n")  # Write upper bounds

    lower_bounds = [lower_bounds_list[i] for i in range(len(lower_bounds_list)) if not np.isnan(lower_bounds_list[i]) and not np.isnan(upper_bounds_list[i])]
    upper_bounds = [upper_bounds_list[i] for i in range(len(upper_bounds_list)) if not np.isnan(lower_bounds_list[i]) and not np.isnan(upper_bounds_list[i])]
    
    return lower_bounds, upper_bounds
#........................................................................................
# PHASE TO TIME ###########################################################################
def orbital_phase_to_time_(ph, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2 ,precision=0.01):
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

def orbital_time_to_phase_(t,iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01):
    
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
                    
    #.............................Now that we know the relation between time, W and phase respectively, interpolate to obtain phase from our input time
    times_to_interpolate = t-min(t)
    
    phase_interpolator = interp1d(time, th, kind='cubic', fill_value="extrapolate")
    phase = phase_interpolator(times_to_interpolate)
    
    w_interpolator = interp1d(time, W_to_interpolate, kind='cubic', fill_value="extrapolate")
    W = w_interpolator(times_to_interpolate)
    
    return phase, t, W
    
    
def fold_pulse_(t, c, sc, period, snr=None, rebin=None):
    
    phase = (t - min(t)) / period - np.floor((t - min(t)) / period)
    idx_sort = np.argsort(phase)
    
    phase_sorted = phase[idx_sort]
    c_sorted = c[idx_sort]
    sc_sorted = sc[idx_sort]
    
    if snr:
        pulse = rebin_snr_(phase_sorted, c_sorted, sc_sorted, snr)
        return pulse
        
    if rebin:
        pulse = rebin_bins_(phase_sorted, c_sorted, sc_sorted, rebin)
        return pulse
        
    else:
        print("Please provide a Minimum signal-to-noise ratio threshold (snr) or a bin time (rebin) to proceed.")


def rebin_snr_(t, x, sy, snr_threshold):

    
    w=[]
    c_bin=[]
    t_bin=[]
    sc_bin=[]

    c_new=[]
    t_new=[]
    sc_new=[]
    
    mask = np.where(sy > 0)[0]
    
    t=t[mask]
    x=x[mask]
    sy=sy[mask]
    
    for i in range(len(x)-1):

        w.append(pow(1/sy[i],2))
        t_bin.append(t[i])
        c_bin.append(x[i])
        sc_bin.append(sy[i])
        
        sc_weight = pow(1/(sum(np.array(w))),0.5)
        c_weight = sum(np.array(c_bin)*np.array(w))/sum(np.array(w))
        
        snr_now = sc_weight/c_weight

        if (snr_now <= snr_threshold):

            w = np.array(w)
            c_bin = np.array(c_bin)
            sc_bin = np.array(sc_bin)
            t_bin = np.array(t_bin)

            sc_new.append(pow(1/(sum(w)),0.5))
            c_new.append(sum(c_bin*w)/sum(w))
            t_new.append(sum(t_bin*w)/sum(w))

            w=[]
            c_bin=[]
            t_bin=[]
            sc_bin=[]
     
        
    return t_new,c_new,sc_new

#.................................................. Rebin by bins
def rebin_bins_(t, x, sy, nbin):

    c_new=[]
    t_new=[]
    sc_new=[]
    
    t, x, sy = preprocess_data(t, x, sy)
    
    for i in range(len(x)-nbin-1):

        w = (pow(1/sy[i*nbin:(i+1)*nbin],2))
        
        if (sum(w) >0):
        
            t_bin = t[i*nbin:(i+1)*nbin]
            c_bin = x[i*nbin:(i+1)*nbin]
            sc_bin = sy[i*nbin:(i+1)*nbin]

            #...............................

            sc_new.append(pow(1/(sum(w)),0.5))
            c_new.append(sum(c_bin*w)/sum(w))
            t_new.append(sum(t_bin*w)/sum(w))

    return t_new,c_new,sc_new

##########################################################################################
################################ THEORETICAL FUNCTIONS ###################################
##########################################################################################

###################################### DOPPLER ##########################################
# Conic orbit
# Spiral
# Disc in orbit
# Spiral in orbit
##########################################################################################

# CONIC ORBIT #############################################################################
def doppler_orbit_theoretical(t, units="keV", show_plot=False, precision_for_phase=0.01):
    print("""

    Computes the Doppler variation expected from orbital movement given a time array in seconds.

    Parameters:
    - t (array-like): Time array in seconds.
    - units (str, optional): Units for the output. Default is "keV". Can be "keV", "s", or "angstrom".
    - show_plot (bool, optional): If True, a plot of the orbit and Doppler evolution is shown and saved. Default is False.
    - precision_for_phase (float, optional): Precision for phase calculation. Default is 0.01.

    Returns:
    - t (array-like): Input time array.
    - x (array-like): Orbital phase array.
    - equation (array-like): Expected Doppler variation.

    """)
    parameter_names=["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis", "inclination", "Rstar", "Mstar1","Mstar2", "wind_vel", "feature"]
    
    fixed_values = manage_parameters(parameter_names, "orbit")
    iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, wind_vel, feature = fixed_values
    #...................................................
    feature_ = {
        "keV": kev_ams/feature,
        "s": feature,  # No conversion needed for seconds
        "amstrong": feature  # No conversion needed for lambda
    }
    
    feature = feature_.get(units, 1)
        
    t0 = min(t)
    #...................................................
    
    abar = semimajor* max(Mstar1,Mstar2)/(Mstar1+Mstar2)
    orbitalperiod_s = orbitalperiod * 24 * 60 * 60
    
    x,_,W = orbital_time_to_phase_(t ,0, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)

    R = abar * (1-eccentricity ** 2)/(1+eccentricity * np.cos((x-periapsis/360) * 2 * np.pi)) #R of ellipse

    v_dop =  - R * Rstar * rsun_m * W * np.sin(2 * np.pi * x ) * np.sin(2 * np.pi * inclination /360)
    v_rad =  wind_vel*1000 * np.cos(2 * np.pi * x ) * np.sin(2 * np.pi * inclination /360)

    vdop = v_dop + v_rad

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
        ax3.plot(x[R>0] * 2 * np.pi, R[R>0],"b")
        ax3.set_theta_direction(-1)
        ax3.set_rlabel_position(90)

        plt.tight_layout()
        
        plt.savefig("orbit_doppler_evolution.png")
    
        del units
        
    return t,x,equation
    
# SPIRAL #########################################################################################
def doppler_spiral_theoretical(t, units="keV", show_plot=False):
    print("""
    Computes the Doppler variation expected from a spiral movement given a time array in seconds.
    
    A logarithmic spiral is a type of spiral that grows in size by a constant factor with each turn. Its equation in polar coordinates is
    r = a * e^(b * θ), where:

    r is the distance from the origin (radius)
    omega is the angle from a reference direction (usually the positive x-axis)
    a is the scale factor that determines how quickly the spiral grows
    b is the rate of rotation, controlling the tightness or looseness of the spiral

    Computes the Doppler variation expected from a spiral movement given a time array in seconds.

    Parameters:
    - t (array-like): Time array in seconds.
    - units (str, optional): Units for the output. Default is "keV". Can be "keV", "s", or "angstrom".
    - show_plot (bool, optional): If True, a plot of the spiral and Doppler evolution is shown and saved. Default is False.
    - A form will appear to input the necessary orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.
    
    Returns:
    - x (array-like): Orbital phase array.
    - equation (array-like): Expected Doppler variation.
    """)


    parameter_names=["iphase_spiral", "semimajor_spiral", "b", "omega", "inclination_spiral", "feature"]
    
    fixed_values = manage_parameters(parameter_names, "spiral")
    iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature = fixed_values
    #...................................................
    feature_ = {
        "keV": kev_ams/feature,
        "s": feature,  # No conversion needed for seconds
        "amstrong": feature  # No conversion needed for lambda
    }
    
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
    
    return t,x,equation
    
# ORBIT IN ORBIT #####################################################################################
def doppler_disc_theoretical(t, units="keV", show_plot=False):
    print("""
    Computes the Doppler variation expected from an orbital movement in a main orbit,
    assuming a ballistic movement of plasma around a compact object or the movement of a mass
    entering an accretion disc.

    Parameters:
    - t (array-like): Time array in seconds.
    - units (str, optional): Units for the output. Default is "keV". Can be "keV", "s", or "angstrom".
    - show_plot (bool, optional): If True, a plot of the disc and Doppler evolution is shown and saved. Default is False.
    - A form will appear to input the necessary orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - t (array-like): Time array.
    - x (array-like): Orbital phase array for the first orbit.
    - x2 (array-like): Orbital phase array for the second orbit.
    - equation (array-like): Expected Doppler variation.
    """)
    parameter_names=["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis", "inclination", "Rstar", "Mstar1", "Mstar2", "iphase2", "semimajor2", "orbitalperiod2", "eccentricity2", "periapsis2", "inclination2",  "Mass3","wind_vel", "feature"]
    
    fixed_values = manage_parameters(parameter_names, "disc")
    iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2, inclination2,  Mass3, wind_vel, feature = fixed_values

    #...................................................
    feature_ = {
        "keV": kev_ams/feature,
        "s": feature,  # No conversion needed for seconds
        "amstrong": feature  # No conversion needed for lambda
    }
    
    feature = feature_.get(units, 1)
    
    t0 = min(t)
    #...................................................
    
    abar = semimajor * max(Mstar1,Mstar2)/(Mstar1+Mstar2)
    orbitalperiod_s = orbitalperiod * 24 * 60 * 60

    x,_,W  = orbital_time_to_phase_(t ,iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, max(Mstar1,Mstar2), min(Mstar1,Mstar2)+Mass3, precision=0.01)
    
    R = (abar * (1 - eccentricity ** 2) / (1 + eccentricity * np.cos((x - periapsis ) * 2 * np.pi)))
    vdop1 =  -R *Rstar * rsun_m * W * np.sin(2 * np.pi * x ) * np.sin(2 * np.pi * inclination /360)
    v_rad =   wind_vel*1000 * np.cos(2 * np.pi * x ) * np.sin(2 * np.pi * inclination /360)
    
    #...................................................
    
    abar2 = semimajor2* min(Mstar1,Mstar2)/(min(Mstar1,Mstar2)+Mass3)
    orbitalperiod_s2 = orbitalperiod2 * 24 * 60 * 60

    x2,_,W2  = orbital_time_to_phase_(t ,iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2, Rstar, min(Mstar2,Mstar1), Mass3, precision=0.01)
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
def doppler_spiral_in_orbit_theoretical(t, units="keV", show_plot=False):
    print("""
    This function requires a time array in seconds and returns the time, orbital phase, and Doppler variation
    expected under the assumption of an orbital movement with a logarithmic spiral component.

    Parameters:
    - t (array-like): Time array in seconds.
    - units (str, optional): Units for the output. Default is "keV". Can be "keV", "s", or "angstrom".
    - show_plot (bool, optional): If True, a plot of the spiral and Doppler evolution is shown and saved.
                                  Default is False.
    - A form will appear to input the necessary orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - x (array-like): Orbital phase array.
    - equation (array-like): Expected Doppler variation.

    Description:
    A logarithmic spiral grows in size by a constant factor with each turn, with the polar coordinates equation:
    r = a * e^(b * θ), where:
    - r is the distance from the origin (radius)
    - θ is the angle from a reference direction (usually the positive x-axis)
    - a is the scale factor that determines how quickly the spiral grows
    - b is the rate of rotation, controlling the tightness or looseness of the spiral
    """)

    parameter_names=["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis", "inclination", "iphase_spiral", "semimajor_spiral", "b", "omega", "inclination_spiral", "Rstar", "Mstar1", "Mstar2", "wind_vel", "feature"]
    
    fixed_values = manage_parameters(parameter_names, "spiral_in_orbit")
    iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, Rstar, Mstar1, Mstar2, wind_vel, feature = fixed_values

    #...................................................
    feature_ = {
        "keV": kev_ams/feature,
        "s": feature,  # No conversion needed for seconds
        "amstrong": feature  # No conversion needed for lambda
    }
    
    feature = feature_.get(units, 1)
    t0 = min(t)
    #...................................................
    
    abar = semimajor* max(Mstar1,Mstar2)/(Mstar1+Mstar2)
    orbitalperiod_s = orbitalperiod * 24 * 60 * 60
    
    x,_,W = orbital_time_to_phase_(t , iphase,semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)

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
     
###################################### STELLAR WIND DENSITY ##############################
# Density in the orbit
# Absorption colum
##########################################################################################
     
# DENSITY IN THE ORBIT #############################################################################
def density_through_orbit_theoretical(resolution=0.01, show_plot=False):
    print("""
    This function helps visualize the density (gr/cm^2) encountered by a compact object along its orbit.
    It assumes a spherically distributed stellar wind based on the CAK model.

    The function returns a time array, an orbital phase array, and the density experienced by the compact
    object as it moves through its orbit.

    If "show_plot=True," a plot of the density through the orbit will be saved under the name
    "density_through_the_orbit.png."

    Parameters:
    - resolution (float, optional): Resolution for the phase array. Default is 0.01.
    - show_plot (bool, optional): If True, plots the density through the orbit. Default is False.
    - A form will appear to input the necessary orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - time (array-like): Time array.
    - phase (array-like): Orbital phase array.
    - density (array-like): Density through the orbit in gr/cm^2.
    """)
    parameter_names = ["semimajor","orbitalperiod" ,"eccentricity", "periapsis", "Rstar","Mstar1","Mstar2","wind_infinite_velocity","Mass_loss_rate","beta" ]
    
    fixed_values = manage_parameters(parameter_names, "density_through_orbit")
    semimajor,orbitalperiod , eccentricity,periapsis, Rstar, Mstar1, Mstar2, wind_infinite_velocity, Mass_loss_rate, beta = fixed_values

    th = np.arange(0,1,resolution)

    _,time,_ = orbital_phase_to_time_(th,0, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,precision=0.01)

    #........................................................
    abar = semimajor * max(Mstar1,Mstar2)/(Mstar1+Mstar2)
    
    M_dot_grams = Mass_loss_rate * msun/(365*24*60*60) #Mdot gr/s
    vinf_cm_s = wind_infinite_velocity*100000 #cm/s
    Rstar_cm = Rstar*rsun_cm #R* in cm
    
    Rorb = (abar*(1-eccentricity**2)/(1+eccentricity*np.cos((th-periapsis/360)*2*np.pi)))* Rstar * rsun_cm #In cm

    v = vinf_cm_s*(1-Rstar_cm/Rorb)**beta
    
    rho = (M_dot_grams/(4 * np.pi * v[Rorb>Rstar_cm] * Rorb[Rorb>Rstar_cm]**2))
    
    #...................................................
    if show_plot:
    
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

        ax1.errorbar(time[Rorb>Rstar_cm], rho,  label='Data', color='b')
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("Density through the orbit gr/cm$^2$")
        ax1.legend()

        ax2.errorbar(th[Rorb>Rstar_cm], rho,  label='Data', color='b')
        ax2.set_xlabel("Orbital phase")
        ax2.set_ylabel("Density through the orbit gr/cm$^2$")
        ax2.legend()

        ax3 = plt.subplot(1, 3, 3, projection='polar')
        ax3.plot(th[Rorb>0] * 2 * np.pi, Rorb[Rorb>0],"b")
        ax3.set_theta_zero_location("N")
        ax3.set_theta_direction(-1)
        ax3.set_rlabel_position(90)

        plt.tight_layout()
        
        plt.savefig("density_through_the_orbit.png")
    #.....................................................

    return  time[Rorb > Rstar_cm], th[Rorb > Rstar_cm],  rho

# ABSOPTION COLUMN #############################################################################
def absorption_column_through_orbit_theoretical(resolution=0.01, show_plot=True):
    print("""
    This function visualizes the column density (NH1, x 10^22 cm^-2) encountered by radiation emitted at each orbital phase as it travels towards an observer. It assumes a spherically distributed, neutral (unionized) stellar wind based on the CAK model.

    Parameters:
    - resolution (float, optional): Resolution for the phase array. Default is 0.01.
    - show_plot (bool, optional): If True, plots the absorption column density through the orbit. Default is False.
    - A form will appear to input the necessary orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - time (array-like): Time array.
    - phase (array-like): Orbital phase array.
    - NH1 (array-like): Absorption column density (NH1, x 10^22 cm^-2) through the orbit.
    
    (Please, take into account that if the distance to the star is smaller than the stellar radius, the result will be 0).
    """)

    parameter_names = ["semimajor","orbitalperiod" ,"eccentricity", "periapsis" ,"inclination", "Rstar","Mstar1","Mstar2","wind_infinite_velocity","Mass_loss_rate","beta" ]
    
    fixed_values = manage_parameters(parameter_names, "absorption_column_through_orbit")
    semimajor, orbitalperiod,eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, wind_infinite_velocity, Mass_loss_rate, beta = fixed_values

    th = np.arange(0,1,resolution)
    
    _,time,_ = orbital_phase_to_time_(th,0, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,precision=0.01)

    #........................................................
    abar = semimajor * max(Mstar1,Mstar2)/(Mstar1+Mstar2)
    
    M_dot_grams = Mass_loss_rate * msun/(365*24*60*60) #Mdot gr/s
    vinf_cm_s = wind_infinite_velocity*100000 #cm/s
    Rstar_cm = Rstar*rsun_cm

    Rorb_plot = (abar*(1-eccentricity**2)/(1+eccentricity*np.cos((th-periapsis/360)*2*np.pi)))*Rstar*rsun_cm #In cm
    
    nh=[]

    for i in range(len(th)):
        
        Rorb = (abar*(1-eccentricity**2)/(1+eccentricity*np.cos((th[i]-periapsis/360)*2*np.pi)))*Rstar*rsun_cm #In cm

        def integrand(z):
        
            alpha = np.arccos(np.cos(th[i]*2*np.pi)*np.cos(inclination*2*np.pi/360))
            x = np.sqrt(Rorb**2+z**2-2*Rorb*z*np.cos(alpha))
            v = (vinf_cm_s*(1-Rstar_cm/x)**beta)
            rho = (M_dot_grams/(4 * np.pi * v * x**2))
            return rho

        ne, _ = quad(integrand, 1, Rstar_cm * 1000)

        nh.append(ne*na/1e22)
        
    nh = np.nan_to_num(nh)
    #...................................................
    if show_plot:
    
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 5))
        
        ax1.errorbar(time, nh,  label='Data', color='b')
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("NH1 (x 10$^{22}$ cm$^{-2}$)")
        ax1.legend()
        
        ax2.errorbar(th, nh,  label='Data', color='b')
        ax2.set_xlabel("Orbital phase")
        ax2.set_ylabel("NH1 (x 10$^{22}$ cm$^{-2}$)")
        ax2.legend()

        ax3 = plt.subplot(1, 3, 3, projection='polar')
        ax3.plot(th * 2 * np.pi, Rorb_plot,"b")
        ax3.set_theta_zero_location("N")
        ax3.set_theta_direction(-1)
        ax3.set_rlabel_position(90)

        plt.tight_layout()
        
        plt.savefig("NH_through_the_orbir.png")
    #...................................................
    
    return time, th, nh
         
# Ionization parameter map ###############################################################


def ionization_map_phase(size_in_Rstar=0, min_color=None, max_color=None, save_plot=False, name="ionization_map"):
    print("""
    Generates a logarithmic ionization parameter map based on the stellar wind density, the luminosity, and orbital parameters.
    The uncolored area represents the X-ray shadow.

    Parameters:
    - size_in_Rstar (float, optional): Extent of the map from the stellar center in stellar radii. Default is 2 times the semimajor.
    - min_color (float, optional): Minimum color scale value for the ionization parameter. Default is None.
    - max_color (float, optional): Maximum color scale value for the ionization parameter. Default is None.
    - save_plot (bool, optional): Boolean to decide if the plot should be saved. Default is False.
    - name (str, optional): Name of the file to save the plot. Default is "ionization_map".
    - A form will appear to input the necessary orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.
    
    Returns:
    - chi_result (pd.DataFrame): DataFrame containing the ionization parameter map.
    - area (float): The calculated area between bounds in cm^2.
    """)

    parameter_names = [
        "phase", "semimajor", "eccentricity", "periapsis", "Rstar", "Mstar1", "Mstar2",
        "wind_infinite_velocity", "Mass_loss_rate", "beta", "luminosity", "bound1", "bound2"
    ]
    
    # Load fixed values
    fixed_values = manage_parameters(parameter_names, "ionization_map_phase")
    phase, semimajor, eccentricity, periapsis, Rstar, Mstar1, Mstar2, wind_infinite_velocity, Mass_loss_rate, beta, luminosity, bound1, bound2 = fixed_values
    
    # Calculate various parameters
    th = np.arange(0, 1, 0.001)
    abar = semimajor * max(Mstar1, Mstar2) / (Mstar1 + Mstar2)
    M_dot_grams = Mass_loss_rate * msun / (365 * 24 * 60 * 60)
    vinf_cm_s = wind_infinite_velocity * 100000
    Rstar_cm = Rstar * rsun_cm
    
    # Orbital radius calculations
    Rorb = (abar * (1 - eccentricity**2) / (1 + eccentricity * np.cos((phase - periapsis / 360) * 2 * np.pi))) * Rstar_cm
    Rorb_plot = (abar * (1 - eccentricity**2) / (1 + eccentricity * np.cos((th - periapsis / 360) * 2 * np.pi))) * Rstar_cm
    
    if size_in_Rstar==0:
        size_in_Rstar = 2*max(Rorb_plot/Rstar_cm)
    
    
    x = np.arange(1, size_in_Rstar, 0.01) * Rstar_cm
    v = vinf_cm_s * (1 - Rstar_cm / x)**beta
    ro = M_dot_grams / (4 * np.pi * v * x**2)
    
    # Calculate electron number density
    na = 6.02214076e23 / 1.00797
    ne = ro * na
    
    #..............................................Angles for shadow

    phase_ns = phase * 2 * np.pi
    phase_ns_degrees = np.degrees(phase_ns)

    alpha = np.arcsin(1/(Rorb/ Rstar_cm))
    alpha_degrees = np.degrees(alpha)

    alpha2 = np.arcsin(1/size_in_Rstar)
    alpha_degrees2 = np.degrees(alpha2)

    gamma_ = 180-(alpha_degrees2 + alpha_degrees)
    rho_ = gamma_+phase_ns_degrees

    rho = rho_*2*np.pi/360
    gamma = gamma_*2*np.pi/360

    phase_touch = (180-90-alpha_degrees)/360

    #..............................................
    
    # Create a DataFrame to store chi results
    chi_result = pd.DataFrame(index=np.round(x / Rstar_cm, 3))
    cmap = plt.get_cmap('rainbow')

    # Compute chi values
    for i in range(len(th)):
        distance = np.sqrt(Rorb**2 + x**2 - 2 * Rorb * x * np.cos(2 * np.pi * (th[i] - phase)))
        chi = luminosity / (ne * abs(distance)**2)
        chi_result[str(round(th[i], 3))] = chi

    # Determine color scale limits
    if not max_color:
        max_color = np.percentile(np.concatenate(chi_result.values), 90)
        print("max color coefficient is", round(max_color,2))  # Corrected print statement
        
    if not min_color:
        min_color = np.percentile(np.concatenate(chi_result.values), 10)
        print("min color coefficient is", round(min_color,2))

    # Initialize plot
    fig, axs = plt.subplots(1, 1, figsize=(20, 10), subplot_kw={'projection': 'polar'})
    norm = Normalize(vmin=np.log(min_color), vmax=np.log(max_color))
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    
    # Calculate area between bounds and plot
    area_between_bounds = 0
    
    bound_limit_1=[]
    bound_limit_2=[]
    bound_limit_3=[]
    bound_limit_4=[]
    phase_limit=[]
    phase_limit2=[]
    #....................................................................................
    th_ = np.arange(phase_touch+phase ,phase + gamma_/360, 0.001)
    for i in range(len(th_)):

        x_ = (Rorb/Rstar_cm) * np.sin(alpha_degrees * 2 * np.pi / 360) / np.sin((180 - (th_[i] - phase) * 360 - alpha_degrees) * 2 * np.pi / 360)
            
        if (x_ >= 1):
                
            x = np.arange(x_, size_in_Rstar, 0.01)*Rstar_cm
            v = vinf_cm_s * (1 - Rstar_cm / x)**beta
            ro = M_dot_grams / (4 * np.pi * v * x**2)
            ne = ro * na
        
            distance = np.sqrt(Rorb**2 + x**2 - 2 * Rorb * x * np.cos(2 * np.pi * (th_[i] - phase)))
            chi = luminosity / (ne * abs(distance)**2)
            x_chi_select = x[(chi >= bound1) & (chi <= bound2)] / Rstar_cm
            colors = cmap(norm(np.log(chi)))
        

            axs.scatter(np.tile(th_[i] * 2 * np.pi, len(x)), x / Rstar_cm, c=colors, cmap='rainbow',alpha=0.3)

            if len(x_chi_select) > 1:
                
                bound_limit_1.append(min(x_chi_select))
                bound_limit_2.append(max(x_chi_select))
                phase_limit.append(th_[i])
                
                if max(np.diff(x_chi_select))>0.02:
                
                    idx3 = np.argmax(np.diff(x_chi_select))
                    idx4 = np.argmax(np.diff(x_chi_select))+1
                    
                    bound_limit_3.append(x_chi_select[idx3])
                    bound_limit_4.append(x_chi_select[idx4])
                    phase_limit2.append(th_[i])
                    
    #....................................................................................
    th_ = np.arange(phase - gamma_/360 ,-phase_touch + phase, 0.001)
    for i in range(len(th_)):
            
        x_ = (Rorb/Rstar_cm)*np.sin(alpha_degrees*2*np.pi/360)/(np.sin((180-(phase-th_[i])*360-alpha_degrees)*2*np.pi/360))
            
        if (x_ >= 1):
                
            x = np.arange(x_, size_in_Rstar, 0.01)*Rstar_cm
            v = vinf_cm_s * (1 - Rstar_cm / x)**beta
            ro = M_dot_grams / (4 * np.pi * v * x**2)
            ne = ro * na
        
            distance = np.sqrt(Rorb**2 + x**2 - 2 * Rorb * x * np.cos(2 * np.pi * (th_[i] - phase)))
            chi = luminosity / (ne * abs(distance)**2)
            x_chi_select = x[(chi >= bound1) & (chi <= bound2)] / Rstar_cm
            colors = cmap(norm(np.log(chi)))
            distance = luminosity/(ne*chi)**0.5

            axs.scatter(np.tile(th_[i] * 2 * np.pi, len(x)), x / Rstar_cm, c=colors, cmap='rainbow',alpha=0.1)
            
            
            if len(x_chi_select) > 1:
            
                bound_limit_1.append(min(x_chi_select))
                bound_limit_2.append(max(x_chi_select))
                phase_limit.append(th_[i])
                
                if max(np.diff(x_chi_select))>0.02:
                
                    idx3 = np.argmax(np.diff(x_chi_select))
                    idx4 = np.argmax(np.diff(x_chi_select))+1
                    
                    bound_limit_3.append(x_chi_select[idx3])
                    bound_limit_4.append(x_chi_select[idx4])
                    phase_limit2.append(th_[i])

                    
    #....................................................................................
    th_ = np.arange(-phase_touch + phase ,phase_touch+phase, 0.001)
    for i in range(len(th_)):

        x = np.arange(1, size_in_Rstar,0.01)*Rstar_cm
        v = vinf_cm_s * (1 - Rstar_cm / x)**beta
        ro = M_dot_grams / (4 * np.pi * v * x**2)
        ne = ro * na
        
        distance = np.sqrt(Rorb**2 + x**2 - 2 * Rorb * x * np.cos(2 * np.pi * (th_[i] - phase)))
        chi = luminosity / (ne * abs(distance)**2)
        x_chi_select = x[(chi >= bound1) & (chi <= bound2)] / Rstar_cm
        colors = cmap(norm(np.log(chi)))
        

        axs.scatter(np.tile(th_[i] * 2 * np.pi, len(x)), x / Rstar_cm, c=colors, cmap='rainbow',alpha=0.1)

        if len(x_chi_select) > 1:
        
            bound_limit_1.append(min(x_chi_select))
            bound_limit_2.append(max(x_chi_select))
            phase_limit.append(th_[i])

            
            if max(np.diff(x_chi_select))>0.02:
                
                idx3 = np.argmax(np.diff(x_chi_select))
                idx4 = np.argmax(np.diff(x_chi_select))+1
                    
                bound_limit_3.append(x_chi_select[idx3])
                bound_limit_4.append(x_chi_select[idx4])
                phase_limit2.append(th_[i])

    
    # Plot additional elements
    axs.plot(np.linspace(0, 2 * np.pi, 100), np.ones(100), color='black')
    axs.plot(th * 2 * np.pi, Rorb_plot / Rstar_cm, color='black', linestyle='--',alpha=0.1)
    axs.plot(phase * 2 * np.pi, Rorb / Rstar_cm, color='black', marker='.')
    
    #Plot  and area between bounds
    ph_lim=np.array(phase_limit)
    ph_lim2=np.array(phase_limit2)
    x_bound1 = np.array(bound_limit_1)
    x_bound2 = np.array(bound_limit_2)
    x_bound3 = np.array(bound_limit_3)
    x_bound4 = np.array(bound_limit_4)
    
    sorted_indices = np.argsort(ph_lim)

    ph_lim = ph_lim[sorted_indices]*2*np.pi
    x_bound1 = x_bound1[sorted_indices]
    x_bound2 = x_bound2[sorted_indices]
    
    sorted_indices2 = np.argsort(ph_lim2)
    
    ph_lim2 = ph_lim2[sorted_indices2]*2*np.pi
    x_bound3 = x_bound3[sorted_indices2]
    x_bound4 = x_bound4[sorted_indices2]

    axs.plot(ph_lim ,x_bound1 ,"ko", markersize=1)
    axs.plot(ph_lim ,x_bound2,"ko", markersize=1)
    
    axs.plot(ph_lim2 ,x_bound3 ,"ko", markersize=1)
    axs.plot(ph_lim2 ,x_bound4,"ko", markersize=1)
    
    dph_lim = np.diff(ph_lim)
    dph_lim2 = np.diff(ph_lim2)

    area1 = 0.5 * np.sum((x_bound1[:-1] * x_bound1[1:]) * np.sin(dph_lim))
    area2 = 0.5 * np.sum((x_bound2[:-1] * x_bound2[1:]) * np.sin(dph_lim))
    
    area3 = 0.5 * np.sum((x_bound3[:-1] * x_bound3[1:]) * np.sin(dph_lim2))
    area4 = 0.5 * np.sum((x_bound4[:-1] * x_bound4[1:]) * np.sin(dph_lim2))
    
    area_sec_1 = np.abs(area2-area1)*Rstar_cm**2
    area_sec_2 = np.abs(area3-area4)*Rstar_cm**2
    
    area_between_bounds = np.abs(area_sec_1-area_sec_2)
       #.......................
    axs.set_theta_direction(-1)
    axs.set_theta_offset(np.pi / 2)
    
    cbar = plt.colorbar(sm, ax=axs, orientation='vertical')
    cbar.set_label(r'Log $\chi$')

    
    if save_plot:
        plt.savefig(f"{name}.png")
    
    return chi_result, area_between_bounds
###################################### ORBITAL PHASE TO TIME ##############################
# Orbital phase to time aproximation (constant areolar velocity)
# Orbital time to phase (constant areolar velocity and interpolation)
##########################################################################################

# PHASE TO TIME ###########################################################################
def orbital_phase_to_time(ph, precision=0.01):
    print("""
    Converts orbital phase array to time array for a compact object orbiting a companion star.
    The compact object moves faster at periastron than at apoastro. The increased orbital speed at
    periastron is primarily due to the conservation of angular momentum, which dictates that as the
    compact object moves closer to the central star, it must travel faster to maintain the total angular
    momentum of the system. This relationship is further influenced by Kepler’s laws of planetary motion,
    which describe how objects sweep out equal areas in equal times and the gravitational force between two
    bodies, which strengthens as they approach each other and weakens as they move apart.

    Parameters:
    - ph (array-like): Orbital phase array.
    - precision (float, optional): Resolution for the phase array. Default is 0.01.
      A form will appear to input the necessary orbital parameters. These parameters will be saved in a
      .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need
      to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - ph (array-like): Orbital phase array (same as input).
    - time (array-like): Time array corresponding to the orbital phase.
    - W (array-like): Angular velocity array corresponding to the orbital phase.
    """)


    #.............................Load parameters
    parameter_names = ["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis", "Rstar", "Mstar1", "Mstar2"]
    fixed_values = manage_parameters(parameter_names, "phase_time")
    
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
    
    
def orbital_time_to_phase(t, precision=0.01):
    print("""
    Converts orbital time array to phase array for a compact object orbiting a companion star.
    The compact object moves faster at periastron than at apoastro. The increased orbital speed at
    periastron is primarily due to the conservation of angular momentum, which dictates that as the
    compact object moves closer to the central star, it must travel faster to maintain the total angular
    momentum of the system. This relationship is further influenced by Kepler’s laws of planetary motion,
    which describe how objects sweep out equal areas in equal times and the gravitational force between two
    bodies, which strengthens as they approach each other and weakens as they move apart.

    Parameters:
    - t (array-like): Orbital time array.
    - precision (float, optional): Resolution for the phase array. Default is 0.01.
      A form will appear to input the necessary orbital parameters. These parameters will be saved in a
      .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need
      to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - ph (array-like): Orbital phase array corresponding to the input time.
    - time (array-like): Time array corresponding to the orbital phase.
    - W (array-like): Angular velocity array corresponding to the orbital phase.

    The increased orbital speed at periastron is primarily due to the conservation of angular momentum,
    which dictates that as the compact object moves closer to the central star, it must travel faster to
    maintain the total angular momentum of the system. This relationship is further influenced by Kepler’s
    laws of planetary motion, which describe how objects sweep out equal areas in equal times and the
    gravitational force between two bodies, which strengthens as they approach each other and weakens
    as they move apart.
    """)


    #.............................Load parameters
    parameter_names = ["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis", "Rstar", "Mstar1", "Mstar2"]
    fixed_values = manage_parameters(parameter_names, "time_phase")
    
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
                    
    #.............................Now that we know the relation between time, W and phase respectively, interpolate to obtain phase from our input time
    times_to_interpolate = t-min(t)
    
    phase_interpolator = interp1d(time, th, kind='cubic', fill_value="extrapolate")
    phase = phase_interpolator(times_to_interpolate)
    
    w_interpolator = interp1d(time, W_to_interpolate, kind='cubic', fill_value="extrapolate")
    W = w_interpolator(times_to_interpolate)
    
    return phase, t, W
##########################################################################################
################################## FITTING FUNCTIONS #####################################
##########################################################################################

###################################### DOPPLER ##########################################
# Conic orbit
# Spiral
# Orbit in orbit
# Spiral in orbit
##########################################################################################

# CONIC ORBIT #############################################################################
def conic_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, wind_vel, feature, units, method_, extended_binsize):
    
    #..........
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
        ph_from_t,_,W = orbital_time_to_phase_(t_to_phase ,iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)
        
        phbins = ph_from_t.reshape(shape_t)
        
        size_phase_bin = np.diff(phbins)
        minsizebin = min(size_phase_bin)
        maxph = max(phbins[-1])
        
        phase = np.arange(0,maxph+10,minsizebin/10)
        _,_,W = orbital_phase_to_time_(phase ,iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)
        
        t_to_phase_puntual = np.mean(t , axis=1)
        phase_puntual ,_, W_puntual = orbital_time_to_phase_(t_to_phase_puntual, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)
        
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
        
        phase_discrete ,_,W = orbital_time_to_phase_(t ,iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)
            
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
                 units="keV", method_="extended", extended_binsize=0.01):
    print("""
    Fits observed orbital modulation data by estimating parameters such as phase, semi-major axis,
    orbital period, eccentricity, inclination, and periapsis.

    The fitting process utilizes a particle swarm optimization (PSO) algorithm, which iteratively improves
    parameter estimates by minimizing the chi-squared difference between observed and predicted data.

    The function supports two fitting methods:
    - "discrete": Suitable for discrete data points (e.g., spectra with small orbital phase ranges; faster).
    - "extended": Suitable for data with varying or extended bin sizes, typical with current instrument
      resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Parameters:
    - x_data (array-like): Time bins of the observed data.
    - y_data (array-like): Observed data points corresponding to each time bin.
    - y_err (array-like, optional): Error associated with each observed data point (default is 0).
    - num_iterations (int, optional): Number of iterations for PSO optimization (default is 3).
    - maxiter (int, optional): Maximum number of iterations for each PSO run (default is 1000).
    - swarmsize (int, optional): Number of particles in the PSO swarm (default is 100).
    - units (str, optional): Units of the observed data (default is "keV").
    - method_ (str, optional): Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize (float, optional): Bin size for the extended method (default is 0.01).
      A form will appear to input the necessary bounds for the orbital parameters. These parameters will be
      saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids
      the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - DataFrame: Contains the best-fit parameters and their standard deviations.
    - ph (array-like): Array of phases corresponding to the predicted data.
    - predicted_data (array-like): Predicted data based on the best-fit parameters.
    - chi_squared (float): Chi-squared statistic weighted by the error, indicating the quality of the fit.
    """)


    #............................................Data prep
    parameter_names = ["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis" ,"inclination", "Rstar", "Mstar1", "Mstar2","wind_vel" ,"feature"]
    #............................................
    t = x_data
    x_data, y_err_weight = define_x_y_sy(x_data,y_data, y_err)

    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(t, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        method_ == "discrete"
        
    #............................................Objective function
    def objective_function(params):

        iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, wind_vel,feature = params
        predicted_data = conic_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity,
                                                     periapsis ,inclination, Rstar, Mstar1, Mstar2, wind_vel, feature, units, method_,extended_binsize)

        chi_squared = chi_squared_weighted(y_data, y_err_weight,predicted_data)

        return chi_squared

    #............................................PS implementation
    lower_bounds, upper_bounds = manage_bounds(parameter_names, "orbit_bounds")
    best_params_list = []
    chi_list = []
    
    for i in range(num_iterations):
        
        best_params, _ = pso(objective_function, lb=lower_bounds, ub=upper_bounds, maxiter = maxiter, swarmsize = swarmsize)

        best_params_list.append(best_params)
        predicted_data = conic_orbit(x_data, *best_params, units,method_,extended_binsize)

        chi_squared = chi_squared_weighted(y_data, y_err_weight,predicted_data)
        
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

    chi_squared = chi_squared_weighted(y_data, y_err_weight,predicted_data)
    #.............................Prepare output
    results = []
    
    for param_name, best_param, std_param in zip(parameter_names, best_params, std_params):
        results.append([best_param, std_param])

    df_results = pd.DataFrame(results, index=parameter_names, columns=['Value', 'Std'])
    df_results.index.name = 'Name of the parameter'
    df_results_transposed = df_results.T

    if method_=="extended":

        t =  np.array([np.mean(i) for i in x_data])
        
        ph,_,_ = orbital_time_to_phase_(t ,
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
        
        ph,_,_ = orbital_time_to_phase_(t ,
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
def fit_orbit_ls(x_data, y_data, y_err=0, units="keV", method_="extended", extended_binsize=0.01):
    print("""
    Fits orbital modulation data by estimating parameters such as phase, semi-major axis, orbital period,
    eccentricity, inclination, and periapsis.

    The fitting process utilizes a traditional Least Squares (LS) method, which fits the observed data by
    minimizing the squared differences between observed and predicted values.

    The function can handle two fitting methods:
    - "discrete": Suitable for discrete data points (e.g., spectra with small orbital phase ranges; faster).
    - "extended": Suitable for data with varying or extended bin sizes, typical with current instrument
      resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Parameters:
    - x_data (array-like): Time bins of the observed data.
    - y_data (array-like): Observed data points corresponding to each time bin.
    - y_err (array-like, optional): Error associated with each observed data point (default is 0).
    - units (str, optional): Units of the observed data (default is "keV").
    - method_ (str, optional): Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize (float, optional): Bin size for the extended method (default is 0.01).
      A form will appear to input the necessary bounds for the orbital parameters. These parameters will be
      saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids
      the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - DataFrame: Contains the best-fit parameters and their errors.
    - ph (array-like): Array of phases corresponding to the predicted data.
    - predicted_data (array-like): Predicted data based on the best-fit parameters.
    - chi_squared (float): Chi-squared statistic weighted by the error, indicating the quality of the fit.
    """)


    #...........................................data prep
    parameter_names = ["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis" ,"inclination", "Rstar", "Mstar1", "Mstar2","wind_vel" ,"feature"]
    
    t=x_data
    
    x_data, y_err_weight = define_x_y_sy(x_data,y_data, y_err)

    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(t, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        method_ == "discrete"
        
    #............................................LS implementation
    lower_bounds, upper_bounds = manage_bounds(parameter_names, "orbit_bounds")
    
    model_func = lambda x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2,wind_vel, feature: conic_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, wind_vel,feature, units=units, method_=method_,extended_binsize=extended_binsize)
    
    try:
        fit_params, fit_covariance = curve_fit(model_func, x_data, y_data, bounds=[lower_bounds, upper_bounds], maxfev=100000)
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
        
        ph,_,_ = orbital_time_to_phase_(t ,
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
        
        ph,_,_ = orbital_time_to_phase_(t ,
        df_results_transposed.iphase.Value,
        df_results_transposed.semimajor.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
    return df_results_transposed,  ph, predicted_data, r_squared

# ORBIT IN ORBIT #############################################################################
def disc_in_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2, inclination2, Mass3, feature, wind_vel, units, method_, extended_binsize):
    t = x_data
    
    # Feature conversion based on units
    feature_ = {
        "keV": kev_ams / feature,
        "s": feature,  # No conversion needed for seconds
        "amstrong": feature  # No conversion needed for lambda
    }
    
    feature = feature_.get(units, 1)
    
    # Calculate average semimajor axis and convert orbital periods to seconds
    abar = semimajor * max(Mstar1, Mstar2) / (Mstar1 + Mstar2 + Mass3)
    orbitalperiod_s = orbitalperiod * 24 * 60 * 60
    
    abar2 = semimajor2 * min(Mstar1, Mstar2) / (min(Mstar1, Mstar2) + Mass3)
    orbitalperiod_s2 = orbitalperiod2 * 24 * 60 * 60
    
    shape_t = t.shape
    t_to_phase = t.reshape(-1)
    
    
    if method_ == "extended":
        t_to_phase_puntual = np.mean(t, axis=1)
        # Main orbit phase calculations
        ph_from_t, _, W = orbital_time_to_phase_(t_to_phase, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)
        phbins = ph_from_t.reshape(shape_t)
        size_phase_bin = np.diff(phbins)
        minsizebin = min(size_phase_bin)
        maxph = max(phbins[-1])
        
        phase = np.arange(0, maxph + 10, minsizebin / 10)
        _, _, W = orbital_phase_to_time_(phase, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)
        
        t_to_phase_puntual = np.mean(t, axis=1)
        phase_puntual, _, W_puntual = orbital_time_to_phase_(t_to_phase_puntual, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)
        
        # Secondary orbit phase calculations
        ph_from_t2, _, W2 = orbital_time_to_phase_(t_to_phase, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2, Rstar, min(Mstar2, Mstar1), Mass3, precision=0.01)
        phbins2 = ph_from_t2.reshape(shape_t)
        size_phase_bin2 = np.diff(phbins2)
        minsizebin2 = min(size_phase_bin2)
        maxph2 = max(phbins2[-1])
        
        phase2 = np.arange(0, maxph2 + 10, minsizebin2 / 10)
        _, _, W2 = orbital_phase_to_time_(phase2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2, Rstar, min(Mstar2, Mstar1), Mass3, precision=0.01)
        
        phase_puntual2, _, W_puntual2 = orbital_time_to_phase_(t_to_phase_puntual, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2, Rstar, min(Mstar2, Mstar1), Mass3, precision=0.01)
        
        # Main orbit velocity calculations
        R = abar * (1 - eccentricity ** 2) / (1 + eccentricity * np.cos((phase - periapsis / 360) * 2 * np.pi))
        v_dop = -R * Rstar * rsun_m * W * np.sin(2 * np.pi * phase) * np.sin(2 * np.pi * inclination / 360)
        v_rad = wind_vel * 1000 * np.cos(2 * np.pi * phase) * np.sin(2 * np.pi * inclination / 360)
        
        # Secondary orbit velocity calculations
        R2 = abar2 * (1 - eccentricity2 ** 2) / (1 + eccentricity2 * np.cos((phase2 - periapsis2 / 360) * 2 * np.pi))
        v_dop2 = -R2 * Rstar * rsun_m * W2 * np.sin(2 * np.pi * phase2) * np.sin(2 * np.pi * inclination2 / 360)
        
        # Extended main orbit calculations
        vdop_bin1 = []
        
        for i in range(len(phbins)):
            if size_phase_bin[i] >= extended_binsize:
                vdop_bin1.append(np.mean(v_dop[(phase >= phbins[i, 0]) & (phase <= phbins[i, 1])]) +
                                 np.mean(v_rad[(phase >= phbins[i, 0]) & (phase <= phbins[i, 1])]))
            else:
                R_puntual = abar * (1 - eccentricity ** 2) / (1 + eccentricity * np.cos((phase_puntual - periapsis / 360) * 2 * np.pi))
                vdop_puntual = -R_puntual * Rstar * rsun_m * W_puntual * np.sin(2 * np.pi * phase_puntual) * np.sin(2 * np.pi * inclination / 360)
                vrad_puntual = wind_vel * 1000 * np.cos(2 * np.pi * phase_puntual) * np.sin(2 * np.pi * inclination / 360)
                vdop_bin1.append(vdop_puntual[i] + vrad_puntual[i])
        
        # Extended secondary orbit calculations
        vdop_bin2 = []
        
        for i in range(len(phbins2)):
            if size_phase_bin2[i] >= extended_binsize:
                vdop_bin2.append(np.mean(v_dop2[(phase2 >= phbins2[i, 0]) & (phase2 <= phbins2[i, 1])]))
            else:
                R_puntual2 = abar2 * (1 - eccentricity2 ** 2) / (1 + eccentricity2 * np.cos((phase_puntual2 - periapsis2 / 360) * 2 * np.pi))
                vdop_puntual2 = -R_puntual2 * Rstar * rsun_m * W_puntual2 * np.sin(2 * np.pi * phase_puntual2) * np.sin(2 * np.pi * inclination2 / 360)
                vdop_bin2.append(vdop_puntual2[i])
    
    if method_ == "discrete":
        # Discrete main orbit calculations
        phase_discrete, _, W_discrete = orbital_time_to_phase_(t, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)
        R_discrete = abar * (1 - eccentricity ** 2) / (1 + eccentricity * np.cos((phase_discrete - periapsis / 360) * 2 * np.pi))
        vdop_discrete = -R_discrete * Rstar * rsun_m * W_discrete * np.sin(2 * np.pi * phase_discrete) * np.sin(2 * np.pi * inclination / 360)
        vrad_discrete = wind_vel * 1000 * np.cos(2 * np.pi * phase_discrete) * np.sin(2 * np.pi * inclination / 360)
        vdop_bin1 = vdop_discrete + vrad_discrete
        
        # Discrete secondary orbit calculations
        phase_discrete2, _, W_discrete2 = orbital_time_to_phase_(t_to_phase, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2, Rstar, Mstar2, Mass3, precision=0.01)
        R_discrete2 = abar2 * (1 - eccentricity2 ** 2) / (1 + eccentricity2 * np.cos((phase_discrete2 - periapsis2 / 360) * 2 * np.pi))
        vdop_discrete2 = -R_discrete2 * Rstar * rsun_m * W_discrete2 * np.sin(2 * np.pi * phase_discrete2) * np.sin(2 * np.pi * inclination2 / 360)
        vdop_bin2 = vdop_discrete2
    
    # Combine Doppler velocities from both orbits
    vdop1 = np.array(vdop_bin1)
    vdop2 = np.array(vdop_bin2)
    vdop = vdop2 + vdop1
    
    # Calculate final equation based on units
    equation_ = {
        "keV": kev_ams / (feature * (c + vdop) / c),
        "s": (feature * (c + vdop) / c),  # No conversion needed for seconds
        "amstrong": (feature * (c + vdop) / c)  # No conversion needed for lambda
    }
    equation = equation_.get(units, 1)
    
    return equation
    

# PS FIT------------------------------------------------------------------------------------------------------
def fit_disc_ps(x_data, y_data, y_err=0, num_iterations=3, maxiter=1000, swarmsize=100,
                units="keV", method_="extended", extended_binsize=0.01):
    print("""
    Fits orbital modulation data by estimating parameters such as phase, semi-major axis, orbital period,
    eccentricity, and inclination for the main orbit, as well as corresponding parameters for a secondary
    orbit (e.g., ballistic capture of matter around a compact object or an accretion disk).

    The fitting process uses a particle swarm optimization (PSO) algorithm, which iteratively improves
    parameter estimates by minimizing the chi-squared difference between observed and predicted data.

    The function can handle two fitting methods:
    - "discrete": Suitable for discrete data points (e.g., spectra with small orbital phase ranges; faster).
    - "extended": Suitable for data with varying or extended bin sizes, typical with current instrument
      resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Parameters:
    - x_data (array-like): Time bins of the observed data.
    - y_data (array-like): Observed data points corresponding to each time bin.
    - y_err (array-like, optional): Error associated with each observed data point (default is 0).
    - num_iterations (int, optional): Number of iterations for the PSO algorithm (default is 3).
    - maxiter (int, optional): Maximum number of iterations for each PSO run (default is 1000).
    - swarmsize (int, optional): Number of particles in the PSO swarm (default is 100).
    - units (str, optional): Units of the observed data (default is "keV").
    - method_ (str, optional): Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize (float, optional): Bin size for the extended method (default is 0.01).
      A form will appear to input the necessary bounds for the orbital parameters. These parameters will be
      saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids
      the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - DataFrame: Contains the best-fit parameters and their standard deviations.
    - ph (array-like): Array of phases corresponding to the predicted data for the main orbit.
    - ph2 (array-like): Array of phases corresponding to the predicted data for the secondary orbit.
    - predicted_data (array-like): Predicted data based on the best-fit parameters.
    - chi_squared (float): Chi-squared statistic weighted by the error, indicating the quality of fit.
    """)


    #............................................data prep
    parameter_names = ["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis" ,"inclination", "Rstar", "Mstar1", "Mstar2", "iphase2", "semimajor2", "orbitalperiod2", "eccentricity2", "periapsis2" ,"inclination2",  "Mass3", "feature","wind_vel"]

    x_data, y_err_weight = define_x_y_sy(x_data,y_data, y_err)
        
    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
       x_data = np.mean(x_data , axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        method_ == "discrete"
        
    #............................................objective function
    def objective_function(params):

        iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2 ,inclination2, Mass3, feature, wind_vel = params

        predicted_data = disc_in_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2 ,inclination2, Mass3, feature, wind_vel, units, method_,extended_binsize)

        chi_squared = chi_squared_weighted(y_data, y_err_weight,predicted_data)

        return chi_squared
        
    #............................................PS implementation
    lower_bounds, upper_bounds = manage_bounds(parameter_names, "disc_bounds")
    best_params_list = []
    chi_list = []
    
    for i in range(num_iterations):
        
        best_params, _ = pso(objective_function, lb=lower_bounds, ub=upper_bounds, maxiter = maxiter, swarmsize = swarmsize)
        
        best_params_list.append(best_params)
        predicted_data = disc_in_orbit(x_data, *best_params, units,method_,extended_binsize)
        chi_squared = chi_squared_weighted(y_data, y_err_weight,predicted_data)
        
        chi_list.append(chi_squared)

    #.............................Collet results
    mean_params = np.mean (best_params_list, axis=0)
    std_params = np.std (best_params_list, axis=0)

    best_iteration = np.argmin(chi_list)
    best_params = best_params_list[best_iteration]
    best_chi = chi_list[best_iteration]

    (iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2 ,inclination2, Mass3, feature, wind_vel) = best_params

    # ...............................Evaluate reults
    predicted_data = disc_in_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2 ,inclination2, Mass3, feature, wind_vel, units, method_,extended_binsize)

    chi_squared = chi_squared_weighted(y_data, y_err_weight,predicted_data)
    #.............................prepare output

    results = []
    for param_name, best_param, std_param in zip(parameter_names, best_params, std_params):
        results.append([best_param, std_param])

    df_results = pd.DataFrame(results, index=parameter_names, columns=['Value', 'Std'])
    df_results.index.name = 'Name of the parameter'

    df_results_transposed = df_results.T
    
    if method_=="extended":

        t =  np.array([np.mean(i) for i in x_data])
        
        ph,_,_ = orbital_time_to_phase_(t ,
        df_results_transposed.iphase.Value,
        df_results_transposed.semimajor.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        df_results_transposed.Mass3.Value + min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
        ph2,_,_ = orbital_time_to_phase_(t ,
        df_results_transposed.iphase2.Value,
        df_results_transposed.semimajor2.Value,
        df_results_transposed.orbitalperiod2.Value,
        df_results_transposed.eccentricity2.Value,
        df_results_transposed.periapsis2.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        df_results_transposed.Mass3.Value + min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
    if method_=="discrete":

        t = x_data
        
        ph,_,_ = orbital_time_to_phase_(t ,
        df_results_transposed.iphase.Value,
        df_results_transposed.semimajor.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        df_results_transposed.Mass3.Value + min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
        ph2,_,_ = orbital_time_to_phase_(t ,
        df_results_transposed.iphase2.Value,
        df_results_transposed.semimajor2.Value,
        df_results_transposed.orbitalperiod2.Value,
        df_results_transposed.eccentricity2.Value,
        df_results_transposed.periapsis2.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        df_results_transposed.Mass3.Value + min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)

    return df_results_transposed, ph, ph2, predicted_data, chi_squared
    
# lS FIT------------------------------------------------------------------------------------------------------

def fit_disc_ls(x_data, y_data, y_err=0, units="keV", method_="extended", extended_binsize=0.01):
    print("""
    Fits orbital modulation data by estimating parameters such as phase, semi-major axis, orbital period,
    eccentricity, and inclination for the main orbit, as well as corresponding parameters for a secondary
    orbit (e.g., ballistic capture of matter around a compact object or an accretion disk).

    The fitting process uses a traditional least squares (LS) method, provided here for completeness despite
    its potential limitations due to the model's complexity.

    The function can handle two fitting methods:
    - "discrete": Suitable for discrete data points (e.g., spectra with small orbital phase ranges; faster).
    - "extended": Suitable for data with varying or extended bin sizes, typical with current instrument
      resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Parameters:
    - x_data (array-like): Time bins of the observed data.
    - y_data (array-like): Observed data points corresponding to each time bin.
    - y_err (array-like, optional): Error associated with each observed data point (default is 0).
    - units (str, optional): Units of the observed data (default is "keV").
    - method_ (str, optional): Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize (float, optional): Bin size for the extended method (default is 0.01).
      A form will appear to input the necessary bounds for the orbital parameters. These parameters will be
      saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids
      the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - DataFrame: Contains the best-fit parameters and their errors.
    - ph (array-like): Array of phases corresponding to the predicted data.
    - predicted_data (array-like): Predicted data based on the best-fit parameters.
    - chi_squared (float): Chi-squared statistic weighted by the error, indicating the quality of fit.
    """)


    #............................................ data prep
    parameter_names = ["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis" ,"inclination", "Rstar", "Mstar1", "Mstar2", "iphase2", "semimajor2", "orbitalperiod2", "eccentricity2", "periapsis2" ,"inclination2",  "Mass3", "feature","wind_vel"]

    x_data, y_err_weight = define_x_y_sy(x_data,y_data, y_err)
        
    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(x_data, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        method_ == "discrete"
        
    #............................................LS implementation
    lower_bounds, upper_bounds = manage_bounds(parameter_names, "disc_bounds")
    
    model_func = lambda x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2 ,inclination2, Mass3, feature, wind_vel : disc_in_orbit(x_data,  iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2 ,inclination2, Mass3, feature,wind_vel,units=units, method_=method_,extended_binsize=extended_binsize)
    
    try:
        fit_params, fit_covariance = curve_fit(model_func, x_data, y_data,  bounds=[lower_bounds, upper_bounds], maxfev=100000)
    except RuntimeError:
        raise RuntimeError("Curve fitting did not converge. Try adjusting the bounds.")
        
    errors = np.sqrt(np.diag(fit_covariance))
    
    #......................................Evaluate results
    predicted_data = model_func(x_data, *fit_params)
    residuals = y_data - predicted_data
    rss = np.sum(residuals**2)
    tss = np.sum((y_data - np.mean(y_data))**2)
    r_squared = 1 - (rss / tss)

    #.............................Prepare output
    results = []
    for param_name, best_param, std_param in zip(parameter_names, fit_params, errors):
        results.append([best_param, std_param])

    df_results = pd.DataFrame(results, index=parameter_names, columns=['Value', 'Std'])
    df_results.index.name = 'Name of the parameter'

    df_results_transposed = df_results.T

    
    if method_=="extended":

        t =  np.array([np.mean(i) for i in x_data])
        
        ph,_,_ = orbital_time_to_phase_(t ,
        df_results_transposed.iphase.Value,
        df_results_transposed.semimajor.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        df_results_transposed.Mass3.Value + min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
        ph2,_,_ = orbital_time_to_phase_(t ,
        df_results_transposed.iphase2.Value,
        df_results_transposed.semimajor2.Value,
        df_results_transposed.orbitalperiod2.Value,
        df_results_transposed.eccentricity2.Value,
        df_results_transposed.periapsis2.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        df_results_transposed.Mass3.Value + min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
    if method_=="discrete":

        t = x_data
        
        ph,_,_ = orbital_time_to_phase_(t ,
        df_results_transposed.iphase.Value,
        df_results_transposed.semimajor.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        df_results_transposed.Mass3.Value + min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
        ph2,_,_ = orbital_time_to_phase_(t ,
        df_results_transposed.iphase2.Value,
        df_results_transposed.semimajor2.Value,
        df_results_transposed.orbitalperiod2.Value,
        df_results_transposed.eccentricity2.Value,
        df_results_transposed.periapsis2.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        df_results_transposed.Mass3.Value + min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
    return df_results_transposed,  ph,ph2, predicted_data, r_squared

# SPIRAL #############################################################################
def spiral(x_data, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units, method_, extended_binsize):

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
        
        phase = np.arange(0,maxph,minsizebin/10)
        
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
                  units="keV", method_="extended", extended_binsize=0.01):
    print("""
    Fits orbital modulation data by estimating parameters for a spiral orbit.

    The fitting process uses a particle swarm optimization (PSO) algorithm, which iteratively improves the
    parameter estimates by minimizing the chi-squared difference between the observed and predicted data.

    The function can handle two fitting methods:
    - "discrete": Suitable for discrete data points (e.g., spectra with small orbital phase ranges; faster).
    - "extended": Suitable for data with varying or extended bin sizes, typical with current instrument
      resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Parameters:
    - x_data (array-like): Time bins of the observed data.
    - y_data (array-like): Observed data points corresponding to each time bin.
    - y_err (array-like, optional): Error associated with each observed data point (default is 0).
    - num_iterations (int, optional): Number of PSO iterations to perform (default is 3).
    - maxiter (int, optional): Maximum number of iterations for PSO (default is 1000).
    - swarmsize (int, optional): Number of particles in the swarm for PSO (default is 100).
    - units (str, optional): Units of the observed data (default is "keV").
    - method_ (str, optional): Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize (float, optional): Bin size for the extended method (default is 0.01).
      A form will appear to input the necessary bounds for the orbital parameters. These parameters will be
      saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids
      the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - DataFrame: Contains the best-fit parameters and their standard deviations.
    - ph (array-like): Array of phases corresponding to the predicted data.
    - predicted_data (array-like): Predicted data based on the best-fit parameters.
    - chi_squared (float): Chi-squared statistic weighted by the error, indicating the quality of fit.
    """)

    #............................................data prep
    parameter_names = ["iphase_spiral", "semimajor_spiral", "b", "omega", "inclination_spiral", "feature"]
    #............................................
    x_data, y_err_weight = define_x_y_sy(x_data,y_data, y_err)
    
    t=x_data
        
    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(t, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        method_ == "discrete"
        
    #............................................Objective function
    def objective_function(params):
    
        iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature = params

        predicted_data = spiral(x_data, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units, method_,extended_binsize)
      
        chi_squared = chi_squared_weighted(y_data, y_err_weight,predicted_data)

        return chi_squared
        
    #............................................ PS implementation
    lower_bounds, upper_bounds = manage_bounds(parameter_names, "spiral_bounds")
    best_params_list = []
    chi_list = []
    
    for i in range(num_iterations):
        
        best_params, _ = pso(objective_function, lb=lower_bounds, ub=upper_bounds, maxiter = maxiter, swarmsize = swarmsize)
        best_params_list.append(best_params)
        predicted_data = spiral(x_data, *best_params, units,method_,extended_binsize)

        chi_squared = chi_squared_weighted(y_data, y_err_weight,predicted_data)
        
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
    chi_squared = chi_squared_weighted(y_data, y_err_weight,predicted_data)
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
def fit_spiral_ls(x_data, y_data, y_err=0, units="keV", method_="extended", extended_binsize=0.01):
    # ............................................ data prep
    print("""
    Fits orbital modulation data by estimating parameters for a spiral orbit.

    The fitting process uses a traditional least squares (LS) method, provided for completeness due to the
    complexity of the model.

    The function can handle two fitting methods:
    - "discrete": Suitable for discrete data points (e.g., spectra with small orbital phase ranges; faster).
    - "extended": Suitable for data with varying or extended bin sizes, typical with current instrument
      resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Parameters:
    - x_data (array-like): Time bins of the observed data.
    - y_data (array-like): Observed data points corresponding to each time bin.
    - y_err (array-like, optional): Error associated with each observed data point (default is 0).
    - units (str, optional): Units of the observed data (default is "keV").
    - method_ (str, optional): Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize (float, optional): Bin size for the extended method (default is 0.01).
      A form will appear to input the necessary bounds for the orbital parameters. These parameters will be
      saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids
      the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - DataFrame: Contains the best-fit parameters and their errors.
    - ph (array-like): Array of phases corresponding to the predicted data.
    - predicted_data (array-like): Predicted data based on the best-fit parameters.
    - chi_squared (float): Chi-squared statistic weighted by the error, indicating the quality of fit.
    """)


    parameter_names = ["iphase_spiral", "semimajor_spiral", "b", "omega", "inclination_spiral", "feature"]

    x_data, y_err_weight = define_x_y_sy(x_data,y_data, y_err)
   
    t=x_data
        
    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(t, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        method_ == "discrete"
    #............................................LS implementation
    lower_bounds, upper_bounds = manage_bounds(parameter_names, "spiral")
    
    model_func = lambda x_data, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature: spiral(x_data,  iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units=units, method_=method_,extended_binsize=extended_binsize)
    try:
        fit_params, fit_covariance = curve_fit(model_func, x_data, y_data, bounds=[lower_bounds, upper_bounds], maxfev=100000)
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


# SPIRAL IN ORBIT #############################################################################
def spiral_orbit(x_data, iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, inclination_orbit, Rstar, Mstar1, Mstar2, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units, method_,extended_binsize):

    t=x_data

    #...................................................
    feature_ = {
        "keV": kev_ams/feature,
        "s": feature,  # No conversion needed for seconds
        "amstrong": feature  # No conversion needed for lambda
    }
    
    feature = feature_.get(units, 1)
    #...................................................

    abar = semimajor_orbit * max(Mstar1,Mstar2) / (Mstar1 + Mstar2)
    orbitalperiod_s = orbitalperiod * 24 * 60 * 60
    
    shape_t = t.shape
    t_to_phase = t.reshape(-1)
    
    #...................................................
    
    if method_=="extended":
    
        ph_from_t,_,W = orbital_time_to_phase_(t_to_phase ,iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)
        
        phbins = ph_from_t.reshape(shape_t)
            
        size_phase_bin_orbit = np.diff(phbins)
        minsizebin = min(size_phase_bin_orbit)
        maxph = max(ph_from_t)
            
        phase = np.arange(0,maxph+10,minsizebin/10)
        _,_,W = orbital_phase_to_time_(phase ,iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)
            
        t_to_phase_puntual = np.mean(t , axis=1)
        phase_puntual ,_, W_puntual = orbital_time_to_phase_(t_to_phase_puntual, iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)


        R = abar * (1-eccentricity ** 2)/(1+eccentricity * np.cos((phase-periapsis/360) * 2 * np.pi)) #R of ellipse
        vdop =  -R * Rstar * rsun_m * W * np.sin(2 * np.pi * phase ) * np.sin(2 * np.pi * inclination_orbit )
    
        phbins_orbit = ((t - min(t[0])) / orbitalperiod_s) + iphase_orbit
    
        vdop_bin_orbit=[]
        
        for i in range(len(phbins_orbit)):
            
            if ( size_phase_bin_orbit[i] >= extended_binsize):
            
                vdop_bin_orbit.append(np.mean(vdop[(phase >= phbins_orbit[i,0]) & (phase <= phbins_orbit[i,1])]))
            
            if ( size_phase_bin_orbit[i] < extended_binsize):

                R_puntual = abar * (1-eccentricity ** 2)/(1+eccentricity * np.cos((phase_puntual-periapsis/360) * 2 * np.pi)) #R of ellipse
                vdop_puntual =  -R_puntual * Rstar * rsun_m * W_puntual * np.sin(2 * np.pi * phase_puntual) * np.sin(2 * np.pi * inclination_orbit )
                vdop_bin_orbit.append(vdop_puntual)
    
    #..................................................SPIRAL
        phbins_spiral =  (t - min(t[0])) * omega + iphase_spiral
        size_phase_bin_spiral = np.diff(phbins_spiral)
        minsizebin_spiral = min(size_phase_bin_spiral)
        maxph_spiral = max(phbins_spiral[-1])
        
        phase_spiral = np.arange(0,maxph_spiral,minsizebin_spiral/10)
        
        R_spiral = semimajor_spiral * np.exp(b * 2 * np.pi * phase_spiral)
        vdop_spiral = -R_spiral * Rstar* rsun_m * omega * np.sin(2 * np.pi * phase_spiral ) * np.sin(2 * np.pi * inclination_spiral/360)
        
        vdop_bin_spiral=[]
        
        for i in range(len(phbins_spiral)):
        
            if ( size_phase_bin_spiral[i] >= extended_binsize):

                vdop_bin_spiral.append(np.mean(vdop_spiral[(phase_spiral >= phbins_spiral[i,0]) & (phase_spiral <= phbins_spiral[i,1])]))
                
            if ( size_phase_bin_spiral[i] < extended_binsize):
            
                phase_puntual_spiral =  (np.mean(t[i])-min(t[0]))* omega + iphase_spiral
                R_puntual_spiral = semimajor_spiral * np.exp(b * 2 * np.pi * phase_puntual_spiral)
                vdop_bin_spiral.append( -R_puntual_spiral *Rstar* rsun_m * omega * np.sin(2 * np.pi * phase_puntual_spiral ) * np.sin(2 * np.pi * inclination_spiral/360))
            
        
    if method_=="discrete":
            #............................................................ DISCRETE MAIN ORBIT
       
        phase_discrete_orbit,_,W = orbital_time_to_phase_(t ,iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)

        R_discrete_orbit = abar * (1-eccentricity ** 2)/(1+eccentricity * np.cos((phase_discrete_orbit-periapsis/360) * 2 * np.pi)) #R of ellipse
        vdop_discrete_orbit =  -R_discrete_orbit * Rstar * rsun_m * W * np.sin(2 * np.pi * phase_discrete_orbit) * np.sin(2 * np.pi * inclination_orbit)
        vdop_bin_orbit = vdop_discrete_orbit
    
        #................................................... DISCRETE SPIRAL
        phase_discrete_spiral =  (t - min(t)) * omega + iphase_spiral
        R_spiral = semimajor_spiral * np.exp(b * 2 * np.pi * phase_discrete_spiral)

        vdop_bin_spiral = -R_spiral * Rstar* rsun_m * omega * np.sin(2 * np.pi * phase_discrete_spiral ) * np.sin(2 * np.pi * inclination_spiral/360)
        #..................................................
        
        
    vdop_spiral =  np.array(vdop_bin_spiral)
    vdop_orbit =  np.array(vdop_bin_orbit)
    
    vdop = vdop_spiral+vdop_orbit

    equation_ = {
    
        "keV": kev_ams/(feature * (c + vdop) / c),
        "s": (feature * (c + vdop) / c),  # No conversion needed for seconds
        "amstrong": (feature * (c + vdop) / c)  # No conversion needed for lambda
    }
    
    equation = equation_.get(units, 1)
    
    return equation

# PS FIT------------------------------------------------------------------------------------------------------
def fit_spiral_in_orbit_ps(x_data, y_data, y_err=0, num_iterations=3, maxiter=1000, swarmsize=100,
                           units="keV", method_="extended", extended_binsize=0.01):
    print("""
    Fits orbital modulation data by estimating parameters for a spiral orbit contained within a main orbit.

    The fitting process uses a particle swarm optimization (PSO) algorithm, which iteratively improves the
    parameter estimates by minimizing the chi-squared difference between the observed and predicted data.

    The function can handle two fitting methods:
    - "discrete": Suitable for discrete data points (e.g., spectra with small orbital phase ranges; faster).
    - "extended": Suitable for data with varying or extended bin sizes, typical with current instrument
      resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Parameters:
    - x_data (array-like): Time bins of the observed data.
    - y_data (array-like): Observed data points corresponding to each time bin.
    - y_err (array-like, optional): Error associated with each observed data point (default is 0).
    - num_iterations (int, optional): Number of PSO iterations to perform (default is 3).
    - maxiter (int, optional): Maximum number of iterations for PSO (default is 1000).
    - swarmsize (int, optional): Number of particles in the swarm for PSO (default is 100).
    - units (str, optional): Units of the observed data (default is "keV").
    - method_ (str, optional): Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize (float, optional): Bin size for the extended method (default is 0.01).
      A form will appear to input the necessary bounds for the orbital parameters. These parameters will be
      saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids
      the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - DataFrame: Contains the best-fit parameters and their errors.
    - ph (array-like): Array of phases corresponding to the predicted data.
    - predicted_data (array-like): Predicted data based on the best-fit parameters.
    - chi_squared (float): Chi-squared statistic weighted by the error, indicating the quality of fit.
    """)
    #.............................Data prep
    parameter_names = ["iphase_orbit", "semimajor_orbit", "orbitalperiod", "eccentricity", "periapsis", "inclination", "Rstar", "Mstar1", "Mstar2", "iphase_spiral", "semimajor_spiral", "b", "omega", "inclination_spiral","feature"]
    x_data, y_err_weight = define_x_y_sy(x_data,y_data, y_err)
    
    t=x_data
        
    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(t, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        method_ == "discrete"
        
    #............................................Objective function
    def objective_function(params):
        iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, inclination_orbit, Rstar, Mstar1, Mstar2, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature,  = params

        predicted_data = spiral_orbit(x_data, iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, inclination_orbit, Rstar, Mstar1, Mstar2, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units, method_,extended_binsize)

        chi_squared = chi_squared_weighted(y_data, y_err_weight,predicted_data)

        return chi_squared
        
#............................................ PS implementation
    lower_bounds, upper_bounds = manage_bounds(parameter_names, "spiralorbit")
    best_params_list = []
    chi_list = []
    
    for i in range(num_iterations):
        
        best_params, _ = pso(objective_function, lb=lower_bounds, ub=upper_bounds, maxiter = maxiter, swarmsize = swarmsize)

        best_params_list.append(best_params)
        predicted_data = spiral_orbit(x_data, *best_params, units,method_,extended_binsize)

        chi_squared = chi_squared_weighted(y_data, y_err_weight,predicted_data)
        
        chi_list.append(chi_squared)
        

    #............................. Collect results
    mean_params = np.mean (best_params_list, axis=0)
    std_params = np.std (best_params_list, axis=0)

    best_iteration = np.argmin(chi_list)
    best_params = best_params_list[best_iteration]
    best_chi = chi_list[best_iteration]

    (iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, inclination_orbit, Rstar, Mstar1, Mstar2, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature) = best_params
    
    #.............................Evaluate results
    predicted_data = spiral_orbit(x_data, iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, inclination_orbit, Rstar, Mstar1, Mstar2, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units, method_, extended_binsize)
    chi_squared = chi_squared_weighted(y_data, y_err_weight,predicted_data)
    
    #.............................Prepare output
    results = []
    for param_name, best_param, std_param in zip(parameter_names, best_params, std_params):
        results.append([best_param, std_param])

    df_results = pd.DataFrame(results, index=parameter_names, columns=['Value', 'Std'])
    df_results.index.name = 'Name of the parameter'

    df_results_transposed = df_results.T

    if method_=="extended":

        t =  np.array([np.mean(i) for i in x_data])
        
        ph,_,_ = orbital_time_to_phase_(t ,
        df_results_transposed.iphase_orbit.Value,
        df_results_transposed.semimajor_orbit.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
    if method_=="discrete":

        t = x_data
        
        ph,_,_ = orbital_time_to_phase_(t ,
        df_results_transposed.iphase_orbit.Value,
        df_results_transposed.semimajor_orbit.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
       
        
    return df_results_transposed, ph, predicted_data, chi_squared

# PS FIT------------------------------------------------------------------------------------------------------
def fit_spiral_in_orbit_ls(x_data, y_data, y_err=0, units="keV", method_="extended", extended_binsize=0.01):
    print("""
    Fits orbital modulation data by estimating parameters for a spiral orbit contained within a main orbit.

    The fitting process uses a traditional least squares (LS) method, provided for completeness due to the
    complexity of the model.

    The function can handle two fitting methods:
    - "discrete": Suitable for discrete data points (e.g., spectra with small orbital phase ranges; faster).
    - "extended": Suitable for data with varying or extended bin sizes, typical with current instrument
      resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Parameters:
    - x_data (array-like): Time bins of the observed data.
    - y_data (array-like): Observed data points corresponding to each time bin.
    - y_err (array-like, optional): Error associated with each observed data point (default is 0).
    - units (str, optional): Units of the observed data (default is "keV").
    - method_ (str, optional): Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize (float, optional): Bin size for the extended method (default is 0.01).
      A form will appear to input the necessary bounds for the orbital parameters. These parameters will be
      saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids
      the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - DataFrame: Contains the best-fit parameters and their errors.
    - ph (array-like): Array of phases corresponding to the predicted data.
    - predicted_data (array-like): Predicted data based on the best-fit parameters.
    - chi_squared (float): Chi-squared statistic weighted by the error, indicating the quality of fit.
    """)

    #............................................Data prep
    parameter_names = ["iphase_orbit", "semimajor_orbit", "orbitalperiod", "eccentricity", "periapsis", "inclination", "Rstar", "Mstar1", "Mstar2", "iphase_spiral", "semimajor_spiral", "b", "omega", "inclination_spiral"]

    x_data, y_err_weight = define_x_y_sy(x_data,y_data, y_err)
        
    t=x_data
        
    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(t, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        method_ == "discrete"
    #............................................ LS implementation
    lower_bounds, upper_bounds = manage_bounds(parameter_names, "spiralorbit")
    
    model_func = lambda x_data,iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, inclination_orbit, Rstar, Mstar1, Mstar2, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature: spiral_orbit( x_data,iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, inclination_orbit, Rstar, Mstar1, Mstar2, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units=units, method_=method_, extended_binsize=extended_binsize)
    
    try:
        fit_params, fit_covariance = curve_fit(model_func, x_data, y_data, bounds=[lower_bounds, upper_bounds], maxfev=100000)
    except RuntimeError:
        raise RuntimeError("Curve fitting did not converge. Try adjusting the bounds.")
        
    errors = np.sqrt(np.diag(fit_covariance))
    
    #............................................ Evaluate result
    predicted_data = model_func(x_data,*fit_params)
    residuals = y_data - predicted_data

    rss = np.sum(residuals**2)
    tss = np.sum((y_data - np.mean(y_data))**2)
    r_squared = 1 - (rss / tss)
    #............................. Prepare output
        
    results = []
    for param_name, best_param, std_param in zip(parameter_names, fit_params, errors):
        results.append([best_param, std_param])

    df_results = pd.DataFrame(results, index=parameter_names, columns=['Value', 'Std'])
    df_results.index.name = 'Name of the parameter'
    df_results_transposed = df_results.T

    
    if method_=="extended":

        t =  np.array([np.mean(i) for i in x_data])
        
        ph,_,_ = orbital_time_to_phase_(t ,
        df_results_transposed.iphase_orbit.Value,
        df_results_transposed.semimajor_orbit.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
    if method_=="discrete":

        t = x_data
        
        ph,_,_ = orbital_time_to_phase_(t ,
        df_results_transposed.iphase_orbit.Value,
        df_results_transposed.semimajor_orbit.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
       
    
    return df_results_transposed,  ph, predicted_data, r_squared
    
###################################### STELLAR WIND DENSITY ##############################
# Absorption colum
##########################################################################################
     
# Absorption colum #############################################################################
def nh_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, Mass_loss_rate, wind_infinite_velocity, beta, method_, extended_binsize):

    def nh_calc(th):

        if not isinstance(th, (list, np.ndarray)):
            th = [th]

        nh = []

        for phase in th:
            Rorb = (abar * (1 - eccentricity**2) / (1 + eccentricity * np.cos((phase - periapsis / 360) * 2 * np.pi))) * Rstar_cm  # in cm

            def integrand(z):
            
                alpha = np.arccos(np.cos(phase * 2 * np.pi) * np.cos(inclination * 2 * np.pi / 360))
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
        
        ph_from_t, _, _ = orbital_time_to_phase_(t_to_phase, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)
        t_to_phase_punctual = np.mean(t, axis=1)
        phase_punctual, _, _ = orbital_time_to_phase_(t_to_phase_punctual, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)

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

        
        phase_discrete, _, _ = orbital_time_to_phase_(t, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)
        
        nh_bin = nh_calc(phase_discrete)

    return np.array(nh_bin,dtype=np.float64)
    
# PS FIT------------------------------------------------------------------------------------------------------
def fit_nh_ps(x_data, y_data, y_err=0, num_iterations=3, maxiter=200, swarmsize=20,
              method_="extended", extended_binsize=0.01):
    print("""
    Fits the column density (NH1, x 10^22 cm^-2) encountered by radiation emitted at each orbital phase
    as it travels towards an observer. Assumes a spherically distributed, neutral (unionized) stellar wind
    based on the CAK model.

    The fitting process uses a particle swarm optimization (PSO) algorithm to minimize the chi-squared
    difference between observed and predicted data points.

    The function supports two fitting methods:
    - "discrete": Suitable for discrete data points (e.g., spectra with small orbital phase ranges; faster).
    - "extended": Suitable for data with varying or extended bin sizes, typical for current instrument
      resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Parameters:
    - x_data (array-like): Time bins of the observed data.
    - y_data (array-like): Observed data points corresponding to each time bin.
    - y_err (array-like, optional): Error associated with each observed data point (default is 0).
    - num_iterations (int, optional): Number of iterations for the PSO algorithm (default is 3).
    - maxiter (int, optional): Maximum number of iterations for each PSO run (default is 200).
    - swarmsize (int, optional): Number of particles in the PSO swarm (default is 20).
    - method_ (str, optional): Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize (float, optional): Bin size for the extended method (default is 0.01).
      A form will appear to input the necessary bounds for the orbital parameters. These parameters will be
      saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids
      the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - df_results_transposed (DataFrame): Transposed DataFrame of the best-fit parameters and their standard
      deviations.
    - ph (array-like): Array of phases corresponding to the predicted data.
    - predicted_data (array-like): Predicted data based on the best-fit parameters.
    - chi_squared (float): Chi-squared statistic weighted by the error, indicating the quality of fit.
    """)

#............................................Data prep.
    advise()
    parameter_names=["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis" ,"inclination", "Rstar", "Mstar1", "Mstar2", "Mdot", "v_inf", "beta"]
    
    t = x_data
    x_data, y_err_weight = define_x_y_sy(x_data,y_data, y_err)
        
    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(t, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        method_ == "discrete"
#............................................Objective function

    def objective_function(params):

        iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, Mdot, v_inf, beta = params
        predicted_data = nh_orbit(x_data,  iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1,
                                     Mstar2, Mdot, v_inf, beta, method_,extended_binsize)

        chi_squared = chi_squared_weighted(y_data, y_err_weight,predicted_data)
         
        return chi_squared
#............................................PS implementation
    lower_bounds, upper_bounds = manage_bounds(parameter_names, "nh")
    best_params_list = []
    chi_list = []

    for i in range(num_iterations):
    

        best_params, _ = pso(objective_function, lb=lower_bounds, ub=upper_bounds, maxiter = maxiter, swarmsize = swarmsize, phig=2)

        best_params_list.append(best_params)
        predicted_data = nh_orbit(x_data, *best_params, method_,extended_binsize)

        chi_squared = chi_squared_weighted(y_data, y_err_weight, predicted_data)

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
    chi_squared = chi_squared_weighted(y_data, y_err_weight,predicted_data)
#............................. Prepare output

    results = []
    for param_name, best_param, std_param in zip(parameter_names, best_params, std_params):
        results.append([best_param, std_param])

    df_results = pd.DataFrame(results, index=parameter_names, columns=['Value', 'Std'])
    df_results.index.name = 'Name of the parameter'
    df_results_transposed = df_results.T

    
    if method_=="extended":

        t =  np.array([np.mean(i) for i in x_data])
        
        ph,_,_ = orbital_time_to_phase_(t ,
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
        
        ph,_,_ = orbital_time_to_phase_(t ,
        df_results_transposed.iphase.Value,
        df_results_transposed.semimajor.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
       
        
    return df_results_transposed, ph, predicted_data, chi_squared


##########################################################################################
####################         PERIOD RELATED FUNCTIONS                 ####################
##########################################################################################
#HELPER FUNCTIONS--------------------------------------------------------------------------
#.................................................. Hardness ratio
def hr(x, y, ex, ey):
    print("""

    Calculates hardness ratio and errors.
    
    hr = (h - l) / (h + l)

    Parameters:
    -Count rate or flux in a hard band.
    -Count rate or flux in a soft band.
    -Errors in the hard band
    -Errors in the soft band

    Returns:
    -Hardness ratio.
    -Error in the hardness ratio.

    """)

    f_value = (x - y) / (x + y)
    df_dx = 2 * y / (x + y)**2
    df_dy = -2 * x / (x + y)**2
    
    sigma_f = np.sqrt((df_dx * ex)**2 + (df_dy * ey)**2)
    
    return f_value, sigma_f

#.................................................. Color ratio
def cr(x, y, ex, ey):
    print("""

    Calculates color ratio and errors.
    
    cr = h / l

    Parameters:
    -Count rate or flux in a hard band.
    -Count rate or flux in a soft band.
    -Errors in the hard band
    -Errors in the soft band

    Returns:
    -Hardness ratio.
    -Error in the hardness ratio.

    """)

    f_value = x / y
    df_dx = 1 / y
    df_dy = -x / (y**2)
    
    sigma_f = np.sqrt((df_dx * ex)**2 + (df_dy * ey)**2)
    
    return f_value, sigma_f


#.................................................. Rebin signal to noise
def rebin_snr(t, x, sy, snr_threshold):
    print("""

    Calculates a rebinned signal-to-noise lightcurve

    Parameters:
    -Time array
    -Count rate or flux
    -Errors
    -Minumun required signal to noise ratio (tipically 0.2-0.05)

    Returns:
    -Time array
    -Rebined lightcurve
    -Errors
    """)
    
    w=[]
    c_bin=[]
    t_bin=[]
    sc_bin=[]

    c_new=[]
    t_new=[]
    sc_new=[]
    
    mask = np.where(sy > 0)[0]
    
    t=t[mask]
    x=x[mask]
    sy=sy[mask]
    
    for i in range(len(x)-1):

        w.append(pow(1/sy[i],2))
        t_bin.append(t[i])
        c_bin.append(x[i])
        sc_bin.append(sy[i])
        
        sc_weight = pow(1/(sum(np.array(w))),0.5)
        c_weight = sum(np.array(c_bin)*np.array(w))/sum(np.array(w))
        
        snr_now = sc_weight/c_weight

        if (snr_now <= snr_threshold):

            w = np.array(w)
            c_bin = np.array(c_bin)
            sc_bin = np.array(sc_bin)
            t_bin = np.array(t_bin)

            sc_new.append(pow(1/(sum(w)),0.5))
            c_new.append(sum(c_bin*w)/sum(w))
            t_new.append(sum(t_bin*w)/sum(w))

            w=[]
            c_bin=[]
            t_bin=[]
            sc_bin=[]
     
        
    return t_new,c_new,sc_new

#.................................................. Rebin by bins
def rebin_bins(t, x, sy, nbin):

    print("""

    Calculates a rebinned lightcurve

    Parameters:
    -Time array
    -Count rate or flux
    -Errors
    -Time bin (equires a larger time bin compared to the one we currently have).

    Returns:
    -Time array
    -Rebined lightcurve
    -Errors

    """)

    c_new=[]
    t_new=[]
    sc_new=[]
    
    t, x, sy = preprocess_data(t, x, sy)
    
    for i in range(len(x)-nbin-1):

        w = (pow(1/sy[i*nbin:(i+1)*nbin],2))
        
        if (sum(w) >0):
        
            t_bin = t[i*nbin:(i+1)*nbin]
            c_bin = x[i*nbin:(i+1)*nbin]
            sc_bin = sy[i*nbin:(i+1)*nbin]

            #...............................

            sc_new.append(pow(1/(sum(w)),0.5))
            c_new.append(sum(c_bin*w)/sum(w))
            t_new.append(sum(t_bin*w)/sum(w))

    return t_new,c_new,sc_new
#.................................................. Periods sliding window
def fold_pulse(t, c, sc, period, snr=None, rebin=None):
    print("""
    Folds a lightcurve data array based on a given period and optionally rebins it.

    Parameters:
    - t: Time array of the lightcurve.
    - c: Count rate or flux array corresponding to the time array.
    - sc: Errors (standard deviation) associated with the count rate or flux.
    - period: Period to fold the lightcurve, in the same units as 't'.

    Optional Parameters (one at least should be specified):
    - snr: Minimum signal-to-noise ratio threshold. If provided, applies signal-to-noise rebinning using `rebin_snr`.
    - rebin: Number of bins to rebin the folded lightcurve. If provided, applies uniform time binning using `rebin_bins`.

    Returns:
    - Depending on the optional parameters:
      - If `snr` is specified: Time array, rebinned folded lightcurve, rebinned errors after signal-to-noise ratio rebinning.
      - If `rebin` is specified: Time array, rebinned folded lightcurve, rebinned errors after uniform time binning.
      
    Notes:
    - Either `snr` or `rebin` must be specified to proceed with the function.
    """)
    
    phase = (t - min(t)) / period - np.floor((t - min(t)) / period)

    idx_sort = np.argsort(phase)
    
    phase_sorted = phase[idx_sort]
    c_sorted = c[idx_sort]
    sc_sorted = sc[idx_sort]
    
    if snr:
        pulse = rebin_snr(phase_sorted, c_sorted, sc_sorted, snr)
        return pulse
        
    if rebin:
        pulse = rebin_bins(phase_sorted, c_sorted, sc_sorted, rebin)
        return pulse
        
    else:
        print("Please provide a Minimum signal-to-noise ratio threshold (snr) or a bin time (rebin) to proceed.")



def preprocess_data(t, x, sy):
    """ Preprocess the data by removing entries where sc <= 0, NaN, or inf. """
    mask = (sy > 0) & np.isfinite(sy)
    t = t[mask]
    x = x[mask]
    sy = sy[mask]
    
    return t, x, sy
#.................................................. Periods sliding window


def period_sliding_window(t, c, sc, window_sec, step_sec, max_period=None, min_period=None, false_alarm_threshold=0.1, rel_high_for_error=0.9, folded_pulses=False, snr_pulse=0.2, nbin_pulse=None):

    print("""
    Performs period analysis using a sliding window approach on a lightcurve dataset.

    Parameters:
    - t: Time array of the lightcurve.
    - c: Count rate or flux array corresponding to the time array.
    - sc: Errors (standard deviation) associated with the count rate or flux.
    - window_sec: Size of the sliding window in seconds for period analysis.
    - step_sec: Step size in seconds between consecutive windows.
    - max_period: Maximum period to consider in the periodogram analysis (optional).
    - min_period: Minimum period to consider in the periodogram analysis (optional).
    - false_alarm_threshold: Threshold value for false alarm probability in Lomb-Scargle periodogram.
    - rel_high_for_error: Relative height for error estimation in peak_widths function.
    - folded_pulses: If True, folds the lightcurve using `fold_pulse` for each identified period.
    - snr_pulse: Minimum signal-to-noise ratio threshold for folding using `fold_pulse`.
    - nbin_pulse: Number of bins for uniform time binning in folding using `fold_pulse`.

    Returns:
    - result: DataFrame containing the results of the period analysis, including periods, frequencies, powers,
              errors in period, errors in power, false alarm probabilities, and time range.
    - pulses: Dictionary containing folded pulse data for each identified period if `folded_pulses` is True.

    Notes:
    - The function performs Lomb-Scargle periodogram analysis within each sliding window of the specified size.
    - It filters the periods based on false alarm probability and sorts them by power.
    - If `folded_pulses` is True, it folds the lightcurve for each identified period using `fold_pulse` and stores the results.
    """)

    def preprocess_data(t, c, sc):
        """ Preprocess the data by removing entries where sc <= 0, NaN, or inf. """
        mask = (sc > 0) & np.isfinite(sc)
        return t[mask], c[mask], sc[mask]

    def lb_period_freq(t, c, sc, step, window, max_period, min_period, false_alarm_threshold):
        min_freq = 1 / max_period if max_period else None
        max_freq = 1 / min_period if min_period else None

        # Lomb-Scargle periodogram
        freq, power = LombScargle(t, c, sc).autopower(maximum_frequency=max_freq, minimum_frequency=min_freq, samples_per_peak=1000)
        ls = LombScargle(t, c, sc)
        
        pos = find_peaks(power)[0]  # Peak positions
        df_list = []

        if len(pos) > 0:
            # False alarm probability for the peaks
            fa_prob = ls.false_alarm_probability(power[pos])
            index_fa = fa_prob < false_alarm_threshold
            peaks_indices = pos[index_fa]

            if len(peaks_indices) > 0:
                
                # Calculate widths and errors for the peaks
                results_widths = peak_widths(power, peaks_indices, rel_height=rel_high_for_error)
               
                # Loop through each peak and store data
                for i, index in enumerate(peaks_indices):
                    # Approximate frequency error as half the width of the peak in frequency
                    freq_error = max(np.diff(freq[int(results_widths[2][i])-3:int(results_widths[2][i])+3]))*3
        
                    power_error = results_widths[1][i]  # Width corresponds to the peak
                    snr=power[index]/np.median(power)

                    df_list.append({
                        'min_time': min(t),
                        'max_time': max(t),
                        'Frequency': freq[index],
                        'Period': 1 / freq[index],
                        'Power': power[index],
                        'Freq_Error': freq_error,  # Properly indexed freq_error
                        'Period_Error': (freq_error / freq[index]**2)**(1/2),  # Error in period
                        'Power_Error': power_error,  # Corresponding power error
                        'False_alarm': fa_prob[index_fa][i],  # False alarm probability
                        'snr': snr
                    })

        # Return DataFrame with results
        df = pd.DataFrame(df_list).reset_index(drop=True)

        if len(df) > 1:
            df = df.sort_values(by='Power', ascending=False).reset_index(drop=True)

        return df

    # Preprocess the data
    t = np.asarray(t)
    c = np.asarray(c)
    sc = np.asarray(sc)
    
    t, c, sc = preprocess_data(t, c, sc)

    # Collect results using sliding window
    periods = {}
    for i in range(0, len(t) - window_sec - 1, step_sec):
        t_window = t[i:i + window_sec]
        c_window = c[i:i + window_sec]
        sc_window = sc[i:i + window_sec] + 1e-9  # Avoid division by zero

        if len(t_window) == 0:
            continue

        df_periods = lb_period_freq(t_window, c_window, sc_window, step_sec, window_sec, max_period, min_period, false_alarm_threshold).reset_index(drop=True)

        if len(df_periods) > 0:
            periods[i] = df_periods

    # Combine results
    if len(periods) >= 1:
        result = pd.concat(periods).reset_index(drop=True)
    else:
        result = None

    # Fold pulses if requested
    pulses = {}
    if folded_pulses and result is not None:
        for i in range(len(result)):
            idx = (t >= result.min_time[i]) & (t <= result.max_time[i])

            t_ = np.array(t[idx])
            c_ = np.array(c[idx])
            sc_ = np.array(sc[idx])

            ph_pulse, c_pulse, sc_pulse = fold_pulse_(t_, c_, sc_, result.Period[i], snr=snr_pulse, rebin=nbin_pulse)
            pulse_data = {'ph_pulse': ph_pulse, 'c_pulse': c_pulse, 'sc_pulse': sc_pulse}

            pulses[i] = pulse_data

    return result, pulses
