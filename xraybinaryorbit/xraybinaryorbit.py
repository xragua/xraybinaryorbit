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
mp=0.5

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
##########################################################################################
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


#.......................................................
#..................................................... Properly prepare imput for fitting functions
def _define_x_y_sy(x_data, y_data, y_err):
    """
    Prepare and validate x, y, and error data for model fitting.

    This function processes the input data by ensuring that `x_data`, `y_data`,
    and `y_err` are properly formatted as NumPy arrays. It also calculates error weights
    (`y_err_weight`) based on the provided uncertainties. If no uncertainties are provided,
    a default weight of 1 is assigned. Additionally, the function handles time binning
    when `x_data` has one more element than `y_data`.

    Parameters
    ----------
    x_data : array-like
        The independent variable data (e.g., time or position data).
    y_data : array-like
        The dependent variable data (e.g., observed values).
    y_err : array-like
        The uncertainties (errors) associated with `y_data`. Can be a 1D array of errors
        or a 2D array (e.g., lower and upper bounds).

    Returns
    -------
    x_data : array-like
        The processed independent variable data, possibly as time bin pairs.
    y_err_weight : array-like
        The calculated or default error weights for the fitting process.
    """
    x_data = np.array(x_data)
    y_data = np.array(y_data)
    y_err = np.array(y_err)
    
    # If uncertainties are provided, calculate the error weights
    if np.sum(y_err) != 0:
        if len(y_err.shape) == 2:
            y_err_weight = np.array(y_err[1] + y_err[0])  # Sum upper and lower bounds if 2D
        elif len(y_err.shape) == 1:
            y_err_weight = np.array(y_err)  # Use errors directly if 1D
    
    # If no uncertainties are provided, set all weights to 1
    if np.sum(y_err) == 0:
        y_err_weight = np.ones(len(y_data))
    
    # If `x_data` has one more element than `y_data`, assume it's time bin edges and create pairs
    if len(y_data) + 1 == len(x_data) and len(x_data.shape) == 1:
        x_data = _time_pairs(x_data)
        
    return np.array(x_data), np.array(y_err_weight)

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

#................................................... Easy way to manage values
def _copy_fields(source_file, destination_file, params2=None):
    """
    Copy and update fields from a source file to a destination file, matching parameter names and updating values accordingly.

    This function reads parameters and values from a `source_file`, compares them with the parameters in a `destination_file`,
    and updates the destination file with values from the source where the parameter names match. If the destination file
    does not exist, it can be created using `params2`, and machine parameters are loaded from the `source_file`.

    Parameters
    ----------
    source_file : str
        The path to the source file, which contains the original parameter names and values.
    destination_file : str
        The path to the destination file where the updated parameter values will be written.
    params2 : list of str, optional
        A list of parameter names to use when creating the destination file if it does not already exist.
        If not provided, and the destination file is missing, a ValueError is raised.

    Raises
    ------
    ValueError
        If the destination file does not exist and `params2` is not provided, or if `params2` is empty.

    Notes
    -----
    - The first row of the files is expected to contain comma-separated parameter names.
    - The second row is expected to contain corresponding values, with 'nan' used to represent missing values.
    - Only fields that match between the source and destination files will be updated.
    - If the destination file does not exist, and `params2` is provided, the file will be created with
      those parameters and initialized with `nan` values before being updated with source data.

    Returns
    -------
    None
        The function writes the updated parameters and values back to the `destination_file`.
    """
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

    
def _load_values_to_interface(param_list, fixed_values, name):
    """
    Create a graphical user interface (GUI) for inputting and updating a list of parameter values.

    This function opens a Tkinter-based GUI form where users can input or modify values for a given list of parameters.
    The form pre-populates with the previous values (`fixed_values`), and upon submission, performs validation checks
    on specific parameters. If the validation passes, the new values are saved and the form is closed.

    Parameters
    ----------
    param_list : list of str
        A list of parameter names to display in the form as input fields.
    fixed_values : list of float
        A list of initial values corresponding to the parameters in `param_list`, used to pre-fill the input fields.
    name : str
        The name of the form or window, used as the title of the Tkinter window.

    Raises
    ------
    ValueError : If the input values do not meet the validation criteria for specific parameters.
        - 'eccentricity' must be between 0 and 0.999999.
        - 'semimajor', 'Mstar1', 'Mstar2', and 'Rstar' must all be positive.

    Returns
    -------
    None
        The function updates the `fixed_values_list` global variable with the new user-input values.
    """


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
    

def _manage_parameters(param_list, name, load_directly=False):
    """
    Manage parameter values by loading from a file, displaying a user input form for modification,
    and saving the updated values back to the file.

    This function first attempts to load previously saved parameter values from a file named `name.txt`.
    If the file is not found, it initializes the values to `NaN`. After loading the values,
    it displays a graphical user interface (GUI) using the `_load_values_to_interface` function
    to allow the user to update the parameter values. Finally, it saves the updated values back to the same file.

    Parameters
    ----------
    param_list : list of str
        A list of parameter names to be managed and saved. This list is used both in the GUI form
        and for saving the names to the output file.
    name : str
        The base name for the text file where the parameters will be saved and loaded from.
        The file is expected to be named `{name}.txt`.

    Returns
    -------
    fixed_values_list : list of float
        A list of the final parameter values, either loaded from the file or updated through user input.
    """

    global fixed_values_list

    # Try to load previous fixed values if available
    file_exists = True
    try:
        with open(f"{name}.txt", "r") as file:
            lines = file.readlines()
            fixed_values_list = [float(val) for val in lines[1].strip().split(",")]
    except FileNotFoundError:
        file_exists = False
        fixed_values_list = [np.nan] * len(param_list)

    # If load_directly is True but file does not exist, still display the GUI
    if not load_directly or not file_exists:
        _load_values_to_interface(param_list, fixed_values_list, name)

    # Save fixed values to file
    with open(f"{name}.txt", "w") as file:
        file.write(",".join(map(str, param_list)) + "\n")  # Write parameter names
        file.write(",".join(map(str, fixed_values_list)) + "\n")  # Write fixed values

    return fixed_values_list


#................................................... Easy way to manage bounds in fitting functions

    
def _load_bounds_to_interface(param_list, lower_bounds, upper_bounds, name):
    """
    Create a graphical user interface (GUI) for inputting and validating lower and upper bounds for parameters.

    This function opens a Tkinter-based GUI form where users can enter or modify both the lower and upper bounds
    for a given list of parameters. The form pre-populates with the previous bounds (`lower_bounds`, `upper_bounds`).

    Parameters
    ----------
    param_list : list of str
        A list of parameter names to display in the form as input fields for setting lower and upper bounds.
    lower_bounds : list of float
        A list of initial lower bound values corresponding to the parameters in `param_list`.
    upper_bounds : list of float
        A list of initial upper bound values corresponding to the parameters in `param_list`.
    name : str
        The name of the form or window, used as the title of the Tkinter window.

    """

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
    
def _manage_bounds(param_list, name, load_directly=False):
    """
    Manage parameter bounds by loading from a file, displaying a user input form for modification, and saving
    the updated bounds back to the file.

    This function first attempts to load previously saved lower and upper bounds for parameters from a file named
    `bounds_{name}.txt`. If the file is not found, it initializes the bounds to `NaN`. The user can modify these bounds
    through a graphical user interface (GUI) provided by `_load_bounds_to_interface`. After modification, the bounds are
    saved back to the same file.

    Parameters
    ----------
    param_list : list of str
        A list of parameter names to be managed and saved. This list is used both in the GUI form
        and for saving the parameter names to the output file.
    name : str
        The base name for the bounds file where the parameters will be saved and loaded from.
        The file is expected to be named `bounds_{name}.txt`.
    load_directly : bool, optional
        If True, the function will load the bounds directly from the file if it exists, without showing the GUI.
        If False (default), the function will display the GUI for the user to modify the bounds, even if the file exists.

    Returns
    -------
    lower_bounds : list of float
        The final list of valid lower bounds (i.e., where both lower and upper bounds are not `NaN`).
    upper_bounds : list of float
        The final list of valid upper bounds (i.e., where both lower and upper bounds are not `NaN`).
    """

    global lower_bounds_list, upper_bounds_list

    file_exists = True

    # Try to load previous bounds if available
    try:
        with open(f"{name}.txt", "r") as file:
            lines = file.readlines()

            param_names = lines[0].strip().split(",")  # Extract parameter names
            lower_bounds_list = [float(val) for val in lines[1].strip().split(",")]
            upper_bounds_list = [float(val) for val in lines[2].strip().split(",")]
            
    except FileNotFoundError:
        file_exists = False
        param_names = param_list
        lower_bounds_list = [np.nan] * len(param_list)
        upper_bounds_list = [np.nan] * len(param_list)

    # If load_directly is False or the file does not exist, display the GUI to allow user modifications
    if not load_directly or not file_exists:
        _load_bounds_to_interface(param_names, lower_bounds_list, upper_bounds_list, name)

    # Save updated bounds to the file
    with open(f"{name}.txt", "w") as file:
        file.write(",".join(map(str, param_names)) + "\n")  # Write parameter names
        file.write(",".join(map(str, lower_bounds_list)) + "\n")  # Write lower bounds
        file.write(",".join(map(str, upper_bounds_list)) + "\n")  # Write upper bounds

    # Filter out NaN values and return valid lower and upper bounds
    lower_bounds = [lower_bounds_list[i] for i in range(len(lower_bounds_list))
                    if not np.isnan(lower_bounds_list[i]) and not np.isnan(upper_bounds_list[i])]
    upper_bounds = [upper_bounds_list[i] for i in range(len(upper_bounds_list))
                    if not np.isnan(lower_bounds_list[i]) and not np.isnan(upper_bounds_list[i])]

    return lower_bounds, upper_bounds

    
    
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
def doppler_orbit_theoretical(t, units="keV", show_plot=False, precision_for_phase=0.01, load_directly=False):
    """
    Computes the Doppler variation expected from orbital movement given a time array in seconds.

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
        If True, displays and saves a plot of the orbit and Doppler evolution. Default is False.
    precision_for_phase : float, optional
        Precision for the phase calculation. Default is 0.01.
        
     Notes
    -----
    A form will appear to input the necessary orbital parameters. These parameters are saved in a .txt file
    in the current directory and automatically loaded in subsequent runs, avoiding the need to re-enter parameters.
    You can modify only those parameters that require adjustment.

    Returns
    -------
    t : array-like
        The input time array.
    x : array-like
        Orbital phase array corresponding to the input time array.
    equation : array-like
        Expected Doppler variation computed for the given orbital movement.
    """

    parameter_names = ["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis", "inclination", "Rstar", "Mstar1", "Mstar2", "wind_vel", "feature"]
    
    fixed_values = _manage_parameters(parameter_names, "orbit",load_directly=load_directly)
    iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, wind_vel, feature = fixed_values
    
    #...................................................
    feature_ = {
        "keV": kev_ams / feature,
        "s": feature,  # No conversion needed for seconds
        "angstrom": feature  # No conversion needed for lambda
    }

    # Raise KeyError if the unit is invalid
    if units not in feature_:
        raise KeyError(f"Invalid unit: {units}. Valid units are 'keV', 's', or 'angstrom'.")
    
    feature = feature_[units]
    
    #...................................................
    
    t0 = min(t)
    
    abar = semimajor * max(Mstar1, Mstar2) / (Mstar1 + Mstar2)
    orbitalperiod_s = orbitalperiod * 24 * 60 * 60
    
    x, _, W = _orbital_time_to_phase(t, 0, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)
    
    R = abar * (1 - eccentricity ** 2) / (1 + eccentricity * np.cos((x - periapsis / 360) * 2 * np.pi))  # R of ellipse

    v_dop = -R * Rstar * rsun_m * W * np.sin(2 * np.pi * x) * np.sin(2 * np.pi * inclination / 360)
    v_rad = wind_vel * 1000 * np.cos(2 * np.pi * x) * np.sin(2 * np.pi * inclination / 360)

    vdop = v_dop + v_rad

    equation_ = {
        "keV": kev_ams / (feature * (c + vdop) / c),
        "s": (feature * (c + vdop) / c),  # No conversion needed for seconds
        "angstrom": (feature * (c + vdop) / c)  # No conversion needed for lambda
    }

    equation = equation_[units]
    
    #...................................................
    
    if show_plot:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
        
        ax1.errorbar(t, equation, label='Data', color='b')
        ax1.set_xlabel("Time")
        ax1.set_ylabel("Expected emission feature evolution")
        ax1.legend()

        ax2.errorbar(x, equation, label='Data', color='b')
        ax2.set_xlabel("Orbital phase")
        ax2.set_ylabel("Expected emission feature evolution")
        ax2.legend()

        ax3 = plt.subplot(1, 3, 3, projection='polar')
        ax3.plot(x[R > 0] * 2 * np.pi, R[R > 0], "b")
        ax3.set_theta_direction(-1)
        ax3.set_rlabel_position(90)

        plt.tight_layout()
        plt.savefig("orbit_doppler_evolution.png")
    
    return t, x, equation

    
# SPIRAL #########################################################################################
def doppler_spiral_theoretical(t, units="keV", show_plot=False, load_directly=False):
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
    
    fixed_values = _manage_parameters(parameter_names, "spiral", load_directly=load_directly)
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
    
    return t,x,equation
    
# ORBIT IN ORBIT #####################################################################################
def doppler_disc_theoretical(t, units="keV", show_plot=False, load_directly=False):
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
    
    fixed_values = _manage_parameters(parameter_names, "disc",load_directly=load_directly)
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
def doppler_spiral_in_orbit_theoretical(t, units="keV", show_plot=False, load_directly=False):
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
    
    fixed_values = _manage_parameters(parameter_names, "spiral_in_orbit",load_directly=load_directly)
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
     
###################################### STELLAR WIND DENSITY ##############################
# Density in the orbit
# Absorption colum
##########################################################################################
     
# DENSITY IN THE ORBIT #############################################################################
def density_through_orbit_theoretical(resolution=0.01, show_plot=False, load_directly=False):
    """
    Visualizes the density (gr/cm^2) encountered by a compact object along its orbit, assuming a spherically
    distributed stellar wind based on the CAK (Castor-Abbott-Klein) model.

    Parameters
    ----------
    resolution : float, optional
        Resolution for the phase array. Default is 0.01.
    show_plot : bool, optional
        If True, displays and saves a plot of the density through the orbit. The plot is saved as
        "density_through_the_orbit.png." Default is False.

    Notes
    -----
    A form will appear to input the necessary orbital parameters. These parameters are saved in a .txt file
    in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter parameters,
    allowing modification of only those that require adjustment.

    Returns
    -------
    time : array-like
        The time array corresponding to the orbital movement.
    phase : array-like
        The orbital phase array.
    density : array-like
        The density encountered by the compact object through the orbit, measured in gr/cm^2.

    """

    parameter_names = ["semimajor","orbitalperiod" ,"eccentricity", "periapsis", "Rstar","Mstar1","Mstar2","wind_infinite_velocity","Mass_loss_rate","beta" ]
    
    fixed_values = _manage_parameters(parameter_names, "density_through_orbit",load_directly=load_directly)
    semimajor,orbitalperiod , eccentricity,periapsis, Rstar, Mstar1, Mstar2, wind_infinite_velocity, Mass_loss_rate, beta = fixed_values

    th = np.arange(0,1,resolution)

    _,time,_ = _orbital_phase_to_time(th,0, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,precision=0.01)

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
def absorption_column_through_orbit_theoretical(resolution=0.01, show_plot=True, load_directly=False):
    """
    Visualizes the column density (NH1, x 10^22 cm^-2) encountered by radiation emitted at each orbital
    phase as it travels towards an observer. Assumes a spherically distributed, neutral (unionized)
    stellar wind based on the CAK (Castor-Abbott-Klein) model.

    Parameters
    ----------
    resolution : float, optional
        Resolution for the phase array. Default is 0.01.
    show_plot : bool, optional
        If True, displays and saves a plot of the absorption column density through the orbit.
        Default is False.

    Notes
    -----
    A form will appear to input the necessary orbital parameters. These parameters are saved in a .txt file
    in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter
    parameters, allowing modification of only those that require adjustment.
    If the distance to the star is smaller than the stellar radius, the result will be 0.

    Returns
    -------
    time : array-like
        Time array corresponding to the orbital movement.
    phase : array-like
        Orbital phase array.
    NH1 : array-like
        Absorption column density (NH1, x 10^22 cm^-2) through the orbit.
    
    """


    parameter_names = ["semimajor","orbitalperiod" ,"eccentricity", "periapsis" ,"inclination", "Rstar","Mstar1","Mstar2","wind_infinite_velocity","Mass_loss_rate","beta" ]
    
    fixed_values = _manage_parameters(parameter_names, "absorption_column_through_orbit",load_directly=load_directly)
    semimajor, orbitalperiod,eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, wind_infinite_velocity, Mass_loss_rate, beta = fixed_values

    th = np.arange(0,1,resolution)
    
    _,time,_ = _orbital_phase_to_time(th,0, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,precision=0.01)

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
        
        plt.savefig("NH_through_the_orbit.png")
    #...................................................
    
    return time, th, nh
         
         

# DENSITY AND LOGCHI #############################################################################

def density_and_ionization_orbital_phase_theoretical(resolution=0.01, size=10, show_plot=True, load_directly=False):
    """
    Calculates and visualizes the density and ionization parameter (log(ξ)) encountered by radiation emitted
    at each orbital phase as it travels towards an observer. Assumes a spherically distributed, neutral stellar
    wind based on the CAK (Castor, Abbott, Klein) model. The density profile and ionization parameter are calculated
    along the path from a neutron star (NS) through the stellar wind of its companion.

    Parameters
    ----------
    resolution : float, optional
        The resolution for the phase array (orbital phase). Default is 0.01.
    size : float, optional
        The scaling factor for the size of the path from the NS through the stellar wind. Determines how far into
        the wind the path is extended. Default is 10.
    show_plot : bool, optional
        If True, displays and saves plots of the density, ionization parameter, and orbital path. Default is True.
    load_directly : bool, optional
        If True, attempts to load previously saved orbital parameters from a file. If False, prompts the user to
        input the parameters. Default is False.

    Returns
    -------
    z : array-like
        Array of distances along the path from the neutron star to the observer (in cm).
    density : array-like
        Density profile of the stellar wind (in cm$^{-3}$) along the path.
    chi : array-like
        Ionization parameter (log(ξ)) calculated at each point along the path.

    """

    parameter_names = [ "orb_phase", "luminosity","semimajor", "eccentricity", "periapsis", "inclination","Rstar", "Mstar1", "Mstar2", "wind_infinite_velocity", "Mass_loss_rate", "beta"]
    
    fixed_values = _manage_parameters(parameter_names, "den_chi_orbphase",load_directly=load_directly)
    orb_phase, luminosity, semimajor, eccentricity, periapsis, inclination ,Rstar, Mstar1, Mstar2, wind_infinite_velocity, Mass_loss_rate, beta = fixed_values

    th = np.arange(0, 1, resolution)
    abar = semimajor * max(Mstar1, Mstar2) / (Mstar1 + Mstar2)#BARICENTER CORRECTION
    
    luminosity_ = luminosity * 1e+32
    M_dot_grams = Mass_loss_rate * msun / (365 * 24 * 60 * 60)  # Mass loss rate in grams/s
    vinf_cm_s = wind_infinite_velocity * 100000  # Wind velocity in cm/s
    Rstar_cm = Rstar * rsun_cm  # Convert stellar radius to cm

    # Calculate orbital radius
    Rorb_plot = (abar * (1 - eccentricity**2) / (1 + eccentricity * np.cos((th - periapsis / 360) * 2 * np.pi))) * Rstar * rsun_cm  # In cm
    Rorb = (abar * (1 - eccentricity**2) / (1 + eccentricity * np.cos((orb_phase - periapsis / 360) * 2 * np.pi))) * Rstar * rsun_cm  # In cm
    
    # Create path from the NS to some distance towards the observer (z=distrance travelled by the emitted radiation fromm the NS towards the observer)
    z = np.arange(0, Rorb * size, Rorb * size / 10000)  # Range for distance in the wind

    # Calculate the angle `alpha`
    alpha = np.arccos(np.cos(orb_phase * 2 * np.pi) * np.cos(inclination * 2 * np.pi / 360))
    
    # Calculate x (distance from z points to donnor)
    cosalpha = np.round(np.cos(alpha),10)
    x = np.sqrt(Rorb**2 + z**2 - 2 * Rorb * z * cosalpha)
    
    # Velocity and density in the wind depending on distance to the donnor (depending on distance travelled from NS)
    v = (vinf_cm_s * (1 - Rstar_cm / x)**beta)
    density = (M_dot_grams / (4 * np.pi * v * x**2  * 1.67E-24* 0.5))  # Renamed from rho to density (mu, mp)
    
    # Calculate the chi parameter for each z
    chi = np.log(luminosity_ / (density * 4 * np.pi * abs(z)**2))
    #...................................................
    if show_plot:
    
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
        
        ax1.errorbar(z/Rstar_cm, density,  label='Data', color='b')
        ax1.set_xlabel("Path from NS (R*)")
        ax1.set_ylabel("Density cm$^{-3}$)")
        ax1.legend()
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        
        ax2.errorbar(z/Rstar_cm, chi,  label='Data', color='b')
        ax2.set_xlabel("Path from NS (R*)")
        ax2.set_ylabel("Ionization Parameter (log(ξ))")
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.legend()

        ax3 = plt.subplot(1, 3, 3, projection='polar')
        ax3.plot(th * 2 * np.pi, Rorb_plot,"b")
        ax3.set_theta_zero_location("N")
        ax3.set_theta_direction(-1)
        ax3.set_rlabel_position(90)

        plt.tight_layout()
        
        plt.savefig("density_and_chi_orbital_phase.png")
    #...................................................
    
    return z, density, chi
         
         
         
# Ionization parameter map ###############################################################


def ionization_map_phase(size_in_Rstar=0, min_color=None, max_color=None, save_plot=False, name="ionization_map", load_directly=False):
    """
    Generates a logarithmic ionization parameter map based on the stellar wind density, luminosity, and
    orbital parameters. The uncolored area in the map represents the X-ray shadow.

    Parameters
    ----------
    size_in_Rstar : float, optional
        Extent of the map from the stellar center in stellar radii. Default is 2 times the semimajor axis.
    min_color : float, optional
        Minimum value for the color scale of the ionization parameter. Default is None.
    max_color : float, optional
        Maximum value for the color scale of the ionization parameter. Default is None.
    save_plot : bool, optional
        If True, saves the generated plot. Default is False.
    name : str, optional
        Name of the file to save the plot. Default is "ionization_map".

    Notes
    -----
    A form will appear to input the necessary orbital parameters. These parameters are saved in a .txt file
    in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter parameters,
    allowing modification of only those that require adjustment.

    Returns
    -------
    chi_result : pd.DataFrame
        DataFrame containing the ionization parameter map.
    area : float
        The calculated area between bounds in cm^2.
    """


    parameter_names = [
        "phase", "semimajor", "eccentricity", "periapsis", "Rstar", "Mstar1", "Mstar2",
        "wind_infinite_velocity", "Mass_loss_rate", "beta", "luminosity", "bound1", "bound2"
    ]
    
    # Load fixed values
    fixed_values = _manage_parameters(parameter_names, "ionization_map_phase",load_directly=load_directly)
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
def orbital_phase_to_time(ph, precision=0.01,load_directly=False):
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
    fixed_values = _manage_parameters(parameter_names, "phase_time",load_directly=load_directly)
    
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
    
    
def orbital_time_to_phase(t, precision=0.01,load_directly=False):
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
    fixed_values = _manage_parameters(parameter_names, "time_phase",load_directly=load_directly)
    
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
def _conic_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, wind_vel, feature, units, method_, extended_binsize):
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
                 units="keV", method_="extended", extended_binsize=0.01,load_directly=False):
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
        predicted_data = _conic_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity,
                                                     periapsis ,inclination, Rstar, Mstar1, Mstar2, wind_vel, feature, units, method_,extended_binsize)

        chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)

        return chi_squared

    #............................................PS implementation
    lower_bounds, upper_bounds = _manage_bounds(parameter_names, "bounds_orbit", load_directly=load_directly)


    best_params_list = []
    chi_list = []
    
    for i in range(num_iterations):
        
        best_params, _ = pso(objective_function, lb=lower_bounds, ub=upper_bounds, maxiter = maxiter, swarmsize = swarmsize)

        best_params_list.append(best_params)
        predicted_data = _conic_orbit(x_data, *best_params, units,method_,extended_binsize)

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
    predicted_data = _conic_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, wind_vel, feature, units, method_,extended_binsize)

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
def fit_orbit_ls(x_data, y_data, y_err=0, units="keV", method_="extended", extended_binsize=0.01,load_directly=False):
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
    lower_bounds, upper_bounds = _manage_bounds(parameter_names, "bounds_orbit",load_directly=load_directly)
    
    model_func = lambda x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2,wind_vel, feature: _conic_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, wind_vel,feature, units=units, method_=method_,extended_binsize=extended_binsize)
    
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

# ORBIT IN ORBIT #############################################################################
def _disc_in_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2, inclination2, Mass3, feature, wind_vel, units, method_, extended_binsize):
    """
    Simulate the Doppler and radial velocities for a binary system with a third mass object (disc in orbit). This private function is
    called by fit_disc_ps and fit_disc_ls.

    This function calculates the Doppler shifts and radial velocities in a system with a binary star (primary and secondary objects)
    and a third body (disc or additional mass). The method can calculate these velocities either for discrete time points or extended
    time bins, depending on the value of the `method_` parameter. The output is the Doppler shift or time shift, based on the
    provided units, such as keV, seconds, or wavelength (angstroms).

    Parameters
    ----------
    x_data : array-like
        The input time data (or time bins) for which Doppler shifts and radial velocities are computed.
    iphase : float
        The initial phase of the main orbit, typically in radians or as a fraction of the orbit period.
    semimajor : float
        The semi-major axis of the main orbit, typically in solar radii or other relevant units.
    orbitalperiod : float
        The orbital period of the main system in days.
    eccentricity : float
        The eccentricity of the main orbit, ranging from 0 (circular) to 1 (parabolic).
    periapsis : float
        The argument of periapsis of the main orbit, describing the orientation of the orbit in the plane of motion.
    inclination : float
        The inclination of the main orbit in degrees, relative to the plane of the sky.
    Rstar : float
        The radius of the star, typically in solar radii.
    Mstar1 : float
        The mass of the primary object (main star), typically in solar masses.
    Mstar2 : float
        The mass of the secondary object, typically in solar masses.
    iphase2 : float
        The initial phase of the secondary orbit, typically in radians or as a fraction of the orbit period.
    semimajor2 : float
        The semi-major axis of the secondary orbit, typically in solar radii.
    orbitalperiod2 : float
        The orbital period of the secondary system in days.
    eccentricity2 : float
        The eccentricity of the secondary orbit, ranging from 0 (circular) to 1 (parabolic).
    periapsis2 : float
        The argument of periapsis of the secondary orbit, describing the orientation of the orbit in the plane of motion.
    inclination2 : float
        The inclination of the secondary orbit in degrees.
    Mass3 : float
        The mass of the third object (disc or additional mass) in solar masses.
    feature : float
        The spectral feature of interest, either in keV or another unit, that will be Doppler-shifted.
    wind_vel : float
        The wind velocity from the star, typically in km/s.
    units : str
        The units for the feature, such as "keV", "s" (seconds), or "amstrong" (angstroms). Determines how the feature is Doppler-shifted.
    method_ : str
        The method used for calculation. Can be "extended" (for time bins) or "discrete" (for individual time points).
    extended_binsize : float
        The bin size for the extended method. This is used to average Doppler shifts and velocities over each time bin.

    Returns
    -------
    equation : array-like
        The Doppler-shifted or time-shifted values of the feature, depending on the `units` parameter.
        This can be in keV, seconds, or wavelength (angstroms), depending on the chosen units.

    Notes
    -----
    - The function performs Doppler shift and radial velocity calculations for both the main binary orbit and the secondary orbit
      associated with the third mass object.
    - For the "extended" method, the function computes average values over time bins, while for the "discrete" method,
      it computes values at individual time points.
    - The Doppler shift accounts for both the motion of the stars and any wind velocities from the star.
    - The final result returns the Doppler-shifted feature, adjusted for both the main and secondary orbits, in the specified units.

    """

    # Feature conversion based on units
    feature_ = {
        "keV": kev_ams / feature,
        "s": feature,  # No conversion needed for seconds
        "amstrong": feature  # No conversion needed for lambda
    }
    
    t=x_data
    
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
        ph_from_t, _, W = _orbital_time_to_phase(t_to_phase, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)
        phbins = ph_from_t.reshape(shape_t)
        size_phase_bin = np.diff(phbins)
        minsizebin = min(size_phase_bin)
        maxph = max(phbins[-1])
        
        phase = np.arange(0, maxph + 10, minsizebin / 10)
        _, _, W = _orbital_phase_to_time(phase, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)
        
        t_to_phase_puntual = np.mean(t, axis=1)
        phase_puntual, _, W_puntual = _orbital_time_to_phase(t_to_phase_puntual, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)
        
        # Secondary orbit phase calculations
        ph_from_t2, _, W2 = _orbital_time_to_phase(t_to_phase, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2, Rstar, min(Mstar2, Mstar1), Mass3, precision=0.01)
        phbins2 = ph_from_t2.reshape(shape_t)
        size_phase_bin2 = np.diff(phbins2)
        minsizebin2 = min(size_phase_bin2)
        maxph2 = max(phbins2[-1])
        
        phase2 = np.arange(0, maxph2 + 10, minsizebin2 / 10)
        _, _, W2 = _orbital_phase_to_time(phase2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2, Rstar, min(Mstar2, Mstar1), Mass3, precision=0.01)
        
        phase_puntual2, _, W_puntual2 = _orbital_time_to_phase(t_to_phase_puntual, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2, Rstar, min(Mstar2, Mstar1), Mass3, precision=0.01)
        
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
        phase_discrete, _, W_discrete = _orbital_time_to_phase(t, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01)
        R_discrete = abar * (1 - eccentricity ** 2) / (1 + eccentricity * np.cos((phase_discrete - periapsis / 360) * 2 * np.pi))
        vdop_discrete = -R_discrete * Rstar * rsun_m * W_discrete * np.sin(2 * np.pi * phase_discrete) * np.sin(2 * np.pi * inclination / 360)
        vrad_discrete = wind_vel * 1000 * np.cos(2 * np.pi * phase_discrete) * np.sin(2 * np.pi * inclination / 360)
        vdop_bin1 = vdop_discrete + vrad_discrete
        
        # Discrete secondary orbit calculations
        phase_discrete2, _, W_discrete2 = _orbital_time_to_phase(t_to_phase, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2, Rstar, Mstar2, Mass3, precision=0.01)
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
                units="keV", method_="extended", extended_binsize=0.01,load_directly=False):
    """
    Fits orbital modulation data by estimating parameters such as phase, semi-major axis, orbital period,
    eccentricity, and inclination for the main orbit, as well as corresponding parameters for a secondary
    orbit (e.g., ballistic capture of matter around a compact object or an accretion disk).

    The fitting process uses a particle swarm optimization (PSO) algorithm, which iteratively improves
    parameter estimates by minimizing the chi-squared difference between observed and predicted data.

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
        Number of iterations for the PSO algorithm. Default is 3.
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
        Array of phases corresponding to the predicted data for the main orbit.
    ph2 : array-like
        Array of phases corresponding to the predicted data for the secondary orbit.
    predicted_data : array-like
        Predicted data based on the best-fit parameters.
    chi_squared : float
        Chi-squared statistic weighted by the error, indicating the quality of fit.
    """



    #............................................data prep
    parameter_names = ["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis" ,"inclination", "Rstar", "Mstar1", "Mstar2", "iphase2", "semimajor2", "orbitalperiod2", "eccentricity2", "periapsis2" ,"inclination2",  "Mass3", "feature","wind_vel"]

    x_data, y_err_weight = _define_x_y_sy(x_data,y_data, y_err)
        
    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
       x_data = np.mean(x_data , axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        print("The number of time points does not allow an extended approach. Changing to discrete")
        method_ = "discrete"
        
    #............................................objective function
    def objective_function(params):

        iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2 ,inclination2, Mass3, feature, wind_vel = params

        predicted_data = _disc_in_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2 ,inclination2, Mass3, feature, wind_vel, units, method_,extended_binsize)

        chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)

        return chi_squared
        
    #............................................PS implementation
    lower_bounds, upper_bounds = _manage_bounds(parameter_names, "bounds_disc",load_directly=load_directly)
    best_params_list = []
    chi_list = []
    
    for i in range(num_iterations):
        
        best_params, _ = pso(objective_function, lb=lower_bounds, ub=upper_bounds, maxiter = maxiter, swarmsize = swarmsize)
        
        best_params_list.append(best_params)
        predicted_data = _disc_in_orbit(x_data, *best_params, units,method_,extended_binsize)
        chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)
        
        chi_list.append(chi_squared)

    #.............................Collet results
    mean_params = np.mean (best_params_list, axis=0)
    std_params = np.std (best_params_list, axis=0)

    best_iteration = np.argmin(chi_list)
    best_params = best_params_list[best_iteration]
    best_chi = chi_list[best_iteration]

    (iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2 ,inclination2, Mass3, feature, wind_vel) = best_params

    # ...............................Evaluate reults
    predicted_data = _disc_in_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2 ,inclination2, Mass3, feature, wind_vel, units, method_,extended_binsize)

    chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)
    #.............................prepare output

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
        df_results_transposed.Mass3.Value + min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
        ph2,_,_ = _orbital_time_to_phase(t ,
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
        
        ph,_,_ = _orbital_time_to_phase(t ,
        df_results_transposed.iphase.Value,
        df_results_transposed.semimajor.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        df_results_transposed.Mass3.Value + min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
        ph2,_,_ = _orbital_time_to_phase(t ,
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

def fit_disc_ls(x_data, y_data, y_err=0, units="keV", method_="extended", extended_binsize=0.01,load_directly=False):
    """
    Fits orbital modulation data by estimating parameters such as phase, semi-major axis, orbital period,
    eccentricity, and inclination for the main orbit, as well as corresponding parameters for a secondary
    orbit (e.g., ballistic capture of matter around a compact object or an accretion disk).

    The fitting process uses a traditional least squares (LS) method, provided here for completeness despite
    its potential limitations due to the model's complexity.

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



    #............................................ data prep
    parameter_names = ["iphase", "semimajor", "orbitalperiod", "eccentricity", "periapsis" ,"inclination", "Rstar", "Mstar1", "Mstar2", "iphase2", "semimajor2", "orbitalperiod2", "eccentricity2", "periapsis2" ,"inclination2",  "Mass3", "feature","wind_vel"]

    x_data, y_err_weight = _define_x_y_sy(x_data,y_data, y_err)
        
    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(x_data, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        print("The number of time points does not allow an extended approach. Changing to discrete")
        method_ = "discrete"
        
    #............................................LS implementation
    lower_bounds, upper_bounds = _manage_bounds(parameter_names, "bounds_disc", load_directly=load_directly)
    
    model_func = lambda x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2 ,inclination2, Mass3, feature, wind_vel : _disc_in_orbit(x_data,  iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, iphase2, semimajor2, orbitalperiod2, eccentricity2, periapsis2 ,inclination2, Mass3, feature,wind_vel,units=units, method_=method_,extended_binsize=extended_binsize)
    
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
        
        ph,_,_ = _orbital_time_to_phase(t ,
        df_results_transposed.iphase.Value,
        df_results_transposed.semimajor.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        df_results_transposed.Mass3.Value + min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
        ph2,_,_ = _orbital_time_to_phase(t ,
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
        
        ph,_,_ = _orbital_time_to_phase(t ,
        df_results_transposed.iphase.Value,
        df_results_transposed.semimajor.Value,
        df_results_transposed.orbitalperiod.Value,
        df_results_transposed.eccentricity.Value,
        df_results_transposed.periapsis.Value,
        df_results_transposed.Rstar.Value,
        max(df_results_transposed.Mstar2.Value, df_results_transposed.Mstar1.Value),
        df_results_transposed.Mass3.Value + min(df_results_transposed.Mstar2.Value,df_results_transposed.Mstar1.Value),  precision=0.01)
        
        ph2,_,_ = _orbital_time_to_phase(t ,
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
def _spiral(x_data, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units, method_, extended_binsize):
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
                  units="keV", method_="extended", extended_binsize=0.01,load_directly=False):
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

        predicted_data = _spiral(x_data, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units, method_,extended_binsize)
      
        chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)

        return chi_squared
        
    #............................................ PS implementation
    lower_bounds, upper_bounds = _manage_bounds(parameter_names, "bounds_spiral",load_directly=load_directly)
    best_params_list = []
    chi_list = []
    
    for i in range(num_iterations):
        
        best_params, _ = pso(objective_function, lb=lower_bounds, ub=upper_bounds, maxiter = maxiter, swarmsize = swarmsize)
        best_params_list.append(best_params)
        predicted_data = _spiral(x_data, *best_params, units,method_,extended_binsize)

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
    predicted_data = _spiral(x_data, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units, method_,extended_binsize)
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
def fit_spiral_ls(x_data, y_data, y_err=0, units="keV", method_="extended", extended_binsize=0.01,load_directly=False):
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
    lower_bounds, upper_bounds = _manage_bounds(parameter_names, "bounds_spiral",load_directly=load_directly)
    
    model_func = lambda x_data, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature: _spiral(x_data,  iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units=units, method_=method_,extended_binsize=extended_binsize)
    
    try:
        fit_params, fit_covariance = curve_fit(model_func, x_data, np.array(y_data), bounds=[lower_bounds, upper_bounds], maxfev=100000)
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
def _spiral_orbit(x_data, iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, inclination_orbit, Rstar, Mstar1, Mstar2, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units, method_, extended_binsize):
    """
    Simulate the Doppler shifts for a system that includes both a main orbital motion and a spiral structure. This
    private function is used in fit_spiral_in_orbit_ps and fit_spiral_in_orbit_ls.

    This function calculates the Doppler shifts in a system that has both an elliptical orbit (main orbit)
    and an expanding spiral structure. The Doppler shifts are computed for both the main orbital motion
    and the spiral structure, and the results are combined. The Doppler shifts can be returned in various
    units, such as keV, seconds, or wavelength (angstroms). The function supports both "extended" and "discrete"
    calculation methods.

    Parameters
    ----------
    x_data : array-like
        The input time data (or time bins) to compute the Doppler shifts for both the main orbit and spiral structure.
    iphase_orbit : float
        The initial phase of the main orbit, typically as a fraction of the orbital period or in radians.
    semimajor_orbit : float
        The semi-major axis of the main orbit, in solar radii or other relevant units.
    orbitalperiod : float
        The orbital period of the main system in days.
    eccentricity : float
        The eccentricity of the main orbit, ranging from 0 (circular) to 1 (parabolic).
    periapsis : float
        The argument of periapsis of the main orbit, describing the orientation of the orbit in the plane of motion.
    inclination_orbit : float
        The inclination of the main orbit in degrees, relative to the plane of the sky.
    Rstar : float
        The radius of the star, typically in solar radii.
    Mstar1 : float
        The mass of the primary object (main star), typically in solar masses.
    Mstar2 : float
        The mass of the secondary object, typically in solar masses.
    iphase_spiral : float
        The initial phase of the spiral structure, typically as a fraction of the orbital period or in radians.
    semimajor_spiral : float
        The semi-major axis of the spiral, representing the initial distance from the center.
    b : float
        The exponential growth factor of the spiral, determining how the radius expands with phase.
    omega : float
        The angular velocity of the spiral, typically in radians per second.
    inclination_spiral : float
        The inclination of the spiral in degrees, relative to the plane of the sky.
    feature : float
        The spectral feature of interest, either in keV or another unit, that will be Doppler-shifted.
    units : str
        The units for the feature, such as "keV", "s" (seconds), or "amstrong" (angstroms). Determines how the feature is Doppler-shifted.
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
    - The function models two distinct components: the main elliptical orbit and a spiral structure that expands with time.
    - For the "extended" method, the function computes average Doppler shift values over time bins, while for the "discrete" method,
      it computes Doppler shifts at individual time points.
    - The Doppler shifts are calculated for both the main orbit and the spiral, and their contributions are combined.
    - The Doppler shifts account for the motion of both the spiral and the orbital system, as well as their respective inclinations.
    - The final result is returned in the specified units (keV, seconds, or angstroms).

    """
    

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
    
        ph_from_t,_,W = _orbital_time_to_phase(t_to_phase ,iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)
        
        phbins = ph_from_t.reshape(shape_t)
            
        size_phase_bin_orbit = np.diff(phbins)
        minsizebin = min(size_phase_bin_orbit)
        maxph = max(ph_from_t)
            
        phase = np.arange(0,maxph+10,max(minsizebin/10,maxph/100000))
        _,_,W = _orbital_phase_to_time(phase ,iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)
            
        t_to_phase_puntual = np.mean(t , axis=1)
        phase_puntual ,_, W_puntual = _orbital_time_to_phase(t_to_phase_puntual, iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)


        R = abar * (1-eccentricity ** 2)/(1+eccentricity * np.cos((phase-periapsis/360) * 2 * np.pi)) #R of ellipse
        vdop =  -R * Rstar * rsun_m * W * np.sin(2 * np.pi * phase ) * np.sin(2 * np.pi * inclination_orbit )
    
        phbins_orbit = ((t - min(t[0])) / orbitalperiod_s) + iphase_orbit
    
        vdop_bin_orbit=[]
        
        for i in range(len(phbins_orbit)):
            
            if ( size_phase_bin_orbit[i] >= extended_binsize):
            
                vdop_bin_orbit.append(np.mean(vdop[(phase >= phbins_orbit[i,0]) & (phase <= phbins_orbit[i,1])]))
            
            if ( size_phase_bin_orbit[i] < extended_binsize):

                R_puntual = abar * (1-eccentricity ** 2)/(1+eccentricity * np.cos((phase_puntual[i]-periapsis/360) * 2 * np.pi)) #R of ellipse
                vdop_puntual =  -R_puntual * Rstar * rsun_m * W_puntual[i] * np.sin(2 * np.pi * phase_puntual[i]) * np.sin(2 * np.pi * inclination_orbit )
                vdop_bin_orbit.append(vdop_puntual)
    
    #..................................................SPIRAL
        phbins_spiral =  (t - min(t[0])) * omega + iphase_spiral
        size_phase_bin_spiral = np.diff(phbins_spiral)
        minsizebin_spiral = min(size_phase_bin_spiral)
        maxph_spiral = max(phbins_spiral[-1])
        
        phase_spiral = np.arange(0,maxph_spiral,max(minsizebin_spiral/10,maxph_spiral/100000))
        
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
       
        phase_discrete_orbit,_,W = _orbital_time_to_phase(t ,iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,  precision=0.01)

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
                           units="keV", method_="extended", extended_binsize=0.01,load_directly=False):
    """
    Fits orbital modulation data by estimating parameters for a spiral orbit contained within a main orbit using
    particle swarm optimization (PSO).

    The fitting process uses a PSO algorithm to iteratively improve parameter estimates by minimizing the
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
        Contains the best-fit parameters and their errors.
    ph : array-like
        Array of phases corresponding to the predicted data.
    predicted_data : array-like
        Predicted data based on the best-fit parameters.
    chi_squared : float
        Chi-squared statistic weighted by the error, indicating the quality of fit.
    """

    #.............................Data prep
    parameter_names = ["iphase_orbit", "semimajor_orbit", "orbitalperiod", "eccentricity", "periapsis", "inclination", "Rstar", "Mstar1", "Mstar2", "iphase_spiral", "semimajor_spiral", "b", "omega", "inclination_spiral","feature"]
    x_data, y_err_weight = _define_x_y_sy(x_data,y_data, y_err)
    
    t=x_data
        
    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(t, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        print("The number of time points does not allow an extended approach. Changing to discrete")
        method_ = "discrete"
        
    #............................................Objective function
    def objective_function(params):
        iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, inclination_orbit, Rstar, Mstar1, Mstar2, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature,  = params

        predicted_data = _spiral_orbit(x_data, iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, inclination_orbit, Rstar, Mstar1, Mstar2, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units, method_,extended_binsize)

        chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)

        return chi_squared
        
#............................................ PS implementation
    lower_bounds, upper_bounds = _manage_bounds(parameter_names, "bounds_spiral_orbit",load_directly=load_directly)
    best_params_list = []
    chi_list = []
    
    for i in range(num_iterations):
        
        best_params, _ = pso(objective_function, lb=lower_bounds, ub=upper_bounds, maxiter = maxiter, swarmsize = swarmsize)

        best_params_list.append(best_params)
        predicted_data = _spiral_orbit(x_data, *best_params, units,method_,extended_binsize)

        chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)
        
        chi_list.append(chi_squared)
        

    #............................. Collect results
    mean_params = np.mean (best_params_list, axis=0)
    std_params = np.std (best_params_list, axis=0)

    best_iteration = np.argmin(chi_list)
    best_params = best_params_list[best_iteration]
    best_chi = chi_list[best_iteration]

    (iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, inclination_orbit, Rstar, Mstar1, Mstar2, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature) = best_params
    
    #.............................Evaluate results
    predicted_data = _spiral_orbit(x_data, iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, inclination_orbit, Rstar, Mstar1, Mstar2, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units, method_, extended_binsize)
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
        
        ph,_,_ = _orbital_time_to_phase(t ,
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
def fit_spiral_in_orbit_ls(x_data, y_data, y_err=0, units="keV", method_="extended", extended_binsize=0.01,load_directly=False):
    """
    Fits orbital modulation data by estimating parameters for a spiral orbit contained within a main orbit using
    a traditional least squares (LS) method. Although the LS method may be limited due to the complexity of the model,
    it is provided for completeness.

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


    #............................................Data prep
    parameter_names = ["iphase_orbit", "semimajor_orbit", "orbitalperiod", "eccentricity", "periapsis", "inclination", "Rstar", "Mstar1", "Mstar2", "iphase_spiral", "semimajor_spiral", "b", "omega", "inclination_spiral"]

    x_data, y_err_weight = _define_x_y_sy(x_data,y_data, y_err)
        
    t=x_data
        
    if method_ == "discrete" and len(np.shape(x_data)) == 2 and len(x_data) == len(y_data):
        x_data = np.mean(t, axis=1)
    
    if len(np.shape(x_data)) == 1 and len(x_data) == len(y_data):
        print("The number of time points does not allow an extended approach. Changing to discrete")
        method_ = "discrete"
    #............................................ LS implementation
    lower_bounds, upper_bounds = _manage_bounds(parameter_names, "bounds_spiral_orbit",load_directly=load_directly)
    
    model_func = lambda x_data,iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, inclination_orbit, Rstar, Mstar1, Mstar2, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature: _spiral_orbit( x_data,iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity, periapsis, inclination_orbit, Rstar, Mstar1, Mstar2, iphase_spiral, semimajor_spiral, b, omega, inclination_spiral, feature, units=units, method_=method_, extended_binsize=extended_binsize)
    
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
        
        ph,_,_ = _orbital_time_to_phase(t ,
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
        
        ph,_,_ = _orbital_time_to_phase(t ,
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
def _nh_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, Mass_loss_rate, wind_infinite_velocity, beta, method_, extended_binsize):
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
              method_="extended", extended_binsize=0.01,load_directly=False):
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
    _advise()
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
        predicted_data = _nh_orbit(x_data,  iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1,
                                     Mstar2, Mdot, v_inf, beta, method_,extended_binsize)

        chi_squared = _chi_squared_weighted(y_data, y_err_weight,predicted_data)
         
        return chi_squared
#............................................PS implementation
    lower_bounds, upper_bounds = _manage_bounds(parameter_names, "bounds_nh",load_directly=load_directly)
    best_params_list = []
    chi_list = []

    for i in range(num_iterations):
    

        best_params, _ = pso(objective_function, lb=lower_bounds, ub=upper_bounds, maxiter = maxiter, swarmsize = swarmsize, phig=2)

        best_params_list.append(best_params)
        predicted_data = _nh_orbit(x_data, *best_params, method_,extended_binsize)

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

    predicted_data = _nh_orbit(x_data, iphase, semimajor, orbitalperiod, eccentricity, periapsis ,inclination, Rstar, Mstar1, Mstar2, Mdot, v_inf, beta, method_,extended_binsize)
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


##########################################################################################
####################         PERIOD RELATED FUNCTIONS                 ####################
##########################################################################################
#HELPER FUNCTIONS--------------------------------------------------------------------------
#.................................................. Hardness ratio
def hr(x, y, ex, ey):
    """
    Calculates the hardness ratio and its associated errors.

    The hardness ratio (HR) is calculated using the formula:
    
    HR = (h - l) / (h + l)

    where:
    - h is the count rate or flux in the hard band.
    - l is the count rate or flux in the soft band.

    Parameters
    ----------
    x : float or array-like
        Count rate or flux in the hard band.
    y : float or array-like
        Count rate or flux in the soft band.
    ex : float or array-like
        Errors associated with the hard band.
    ey : float or array-like
        Errors associated with the soft band.

    Returns
    -------
    hr : float or array-like
        The calculated hardness ratio.
    hr_error : float or array-like
        The error in the hardness ratio.
    """


    f_value = (x - y) / (x + y)
    df_dx = 2 * y / (x + y)**2
    df_dy = -2 * x / (x + y)**2
    
    sigma_f = np.sqrt((df_dx * ex)**2 + (df_dy * ey)**2)
    
    return f_value, sigma_f

#.................................................. Color ratio
def cr(x, y, ex, ey):
    """
    Calculates the color ratio and its associated errors.

    The color ratio (CR) is calculated using the formula:
    
    CR = h / l

    where:
    - h is the count rate or flux in the hard band.
    - l is the count rate or flux in the soft band.

    Parameters
    ----------
    x : float or array-like
        Count rate or flux in the hard band.
    y : float or array-like
        Count rate or flux in the soft band.
    ex : float or array-like
        Errors associated with the hard band.
    ey : float or array-like
        Errors associated with the soft band.

    Returns
    -------
    cr : float or array-like
        The calculated color ratio.
    cr_error : float or array-like
        The error in the color ratio.
    """


    f_value = x / y
    df_dx = 1 / y
    df_dy = -x / (y**2)
    
    sigma_f = np.sqrt((df_dx * ex)**2 + (df_dy * ey)**2)
    
    return f_value, sigma_f


#.................................................. Rebin signal to noise
def rebin_snr(t, x, sy, snr_threshold):
    """
    Calculates a rebinned signal-to-noise ratio (SNR) lightcurve.

    Parameters
    ----------
    t : array-like
        Time array.
    x : array-like
        Count rate or flux array.
    sy : array-like
        Error array associated with the count rate or flux.
    snr_threshold : float
        Minimum required signal-to-noise ratio (typically between 0.2 and 0.05).

    Returns
    -------
    t_rebinned : array-like
        Rebinned time array.
    x_rebinned : array-like
        Rebinned lightcurve.
    sy_rebinned : array-like
        Rebinned errors.
    """

    
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
def rebin_bins(t, c, sc, nbin):
    """
    Calculates a rebinned lightcurve using a specified number of time bins.

    Parameters
    ----------
    t : array-like
        Time array.
    x : array-like
        Count rate or flux array.
    sy : array-like
        Error array associated with the count rate or flux.
    nbin : int
        Number of time bins for rebinning (requires a larger time bin compared to the original one).

    Returns
    -------
    t_rebinned : array-like
        Rebinned time array.
    x_rebinned : array-like
        Rebinned lightcurve.
    sy_rebinned : array-like
        Rebinned errors.
    """


    c_new=[]
    t_new=[]
    sc_new=[]
    
    t, c, sc = preprocess_data(t, c, sc)
    
    for i in range(len(c)-nbin-1):

        w = (pow(1/sc[i*nbin:(i+1)*nbin],2))
        
        if (sum(w) >0):
        
            t_bin = t[i*nbin:(i+1)*nbin]
            c_bin = c[i*nbin:(i+1)*nbin]
            sc_bin = sc[i*nbin:(i+1)*nbin]

            #...............................

            sc_new.append(pow(1/(sum(w)),0.5))
            c_new.append(sum(c_bin*w)/sum(w))
            t_new.append(sum(t_bin*w)/sum(w))

    return t_new,c_new,sc_new
#.................................................. Periods sliding window
def fold_pulse(t, c, sc, period, snr=None, rebin=None):
    """
    Folds a lightcurve data array based on a given period and optionally rebins it using either signal-to-noise
    ratio (SNR) or uniform time binning.

    Parameters
    ----------
    t : array-like
        Time array of the lightcurve.
    c : array-like
        Count rate or flux array corresponding to the time array.
    sc : array-like
        Errors (standard deviation) associated with the count rate or flux.
    period : float
        Period to fold the lightcurve, in the same units as the time array `t`.

    Optional Parameters
    -------------------
    snr : float, optional
        Minimum signal-to-noise ratio threshold. If provided, applies signal-to-noise rebinning using `rebin_snr`.
    rebin : int, optional
        Number of bins to rebin the folded lightcurve. If provided, applies uniform time binning using `rebin_bins`.

    Returns
    -------
    If `snr` is specified:
        t_rebinned : array-like
            Rebinned time array after signal-to-noise ratio rebinning.
        c_rebinned : array-like
            Rebinned folded lightcurve.
        sc_rebinned : array-like
            Rebinned errors after signal-to-noise ratio rebinning.

    If `rebin` is specified:
        t_rebinned : array-like
            Rebinned time array after uniform time binning.
        c_rebinned : array-like
            Rebinned folded lightcurve.
        sc_rebinned : array-like
            Rebinned errors after uniform time binning.

    Notes
    -----
    Either `snr` or `rebin` must be specified to proceed with the function.
    """

    
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


def period_sliding_window(t, c, sc, window_sec, step_sec, max_period=None, min_period=None, false_alarm_threshold=0.1,
                          rel_high_for_error=0.9, folded_pulses=False, snr_pulse=0.2, nbin_pulse=None):
    """
    Performs period analysis using a sliding window approach on a lightcurve dataset.

    Parameters
    ----------
    t : array-like
        Time array of the lightcurve.
    c : array-like
        Count rate or flux array corresponding to the time array.
    sc : array-like
        Errors (standard deviation) associated with the count rate or flux.
    window_sec : float
        Size of the sliding window in seconds for period analysis.
    step_sec : float
        Step size in seconds between consecutive windows.
    max_period : float, optional
        Maximum period to consider in the periodogram analysis. Default is None.
    min_period : float, optional
        Minimum period to consider in the periodogram analysis. Default is None.
    false_alarm_threshold : float, optional
        Threshold value for false alarm probability in Lomb-Scargle periodogram. Default is 0.1.
    rel_high_for_error : float, optional
        Relative height for error estimation in the `peak_widths` function. Default is 0.9.
    folded_pulses : bool, optional
        If True, folds the lightcurve using `fold_pulse` for each identified period. Default is False.
    snr_pulse : float, optional
        Minimum signal-to-noise ratio threshold for folding using `fold_pulse`. Default is 0.2.
    nbin_pulse : int, optional
        Number of bins for uniform time binning in folding using `fold_pulse`. Default is None.

    Returns
    -------
    result : DataFrame
        DataFrame containing the results of the period analysis, including:
        - periods,
        - frequencies,
        - powers,
        - errors in period,
        - errors in power,
        - false alarm probabilities,
        - time range.
    pulses : dict, optional
        Dictionary containing folded pulse data for each identified period if `folded_pulses` is True.

    Notes
    -----
    - For each window, the function returns a list of all peaks found in the periodogram.
    - The function performs Lomb-Scargle periodogram analysis within each sliding window of the specified size.
    - It filters periods based on false alarm probability and sorts them by power.
    - If `folded_pulses` is True, it folds the lightcurve for each identified period using `fold_pulse`
      and stores the results.
    """


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
                    
                    # Ensure that the frequency slice contains enough elements
                    start_idx = max(int(results_widths[2][i]) - 3, 0)
                    end_idx = min(int(results_widths[2][i]) + 3, len(freq) - 1)
                    
                    if end_idx > start_idx:  # Ensure there are at least two elements
                        freq_error = max(np.diff(freq[start_idx:end_idx])) * 3
                    else:
                        freq_error = 0  # Assign 0 if there's not enough data to calculate frequency error
                
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
    
    
    min_period_ = 2 * np.mean(np.diff(t))
    
    if (min_period_ >= min_period):

        print(f"\033[91mWarning: The minimum period given is too low to be accurately captured with the current sampling rate. According to Nyquist-Shannon Sampling Theorem, the minimum period allowed to be identified is {min_period_:.4f} seconds.\033[0m")

    
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

            ph_pulse, c_pulse, sc_pulse = fold_pulse(t_, c_, sc_, result.Period[i], snr=snr_pulse, rebin=nbin_pulse)
            pulse_data = {'ph_pulse': ph_pulse, 'c_pulse': c_pulse, 'sc_pulse': sc_pulse}

            pulses[i] = pulse_data

    return result, pulses
