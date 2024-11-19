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
    
def _manage_parameters(param_list, name, load_directly=False, parameter_list=None):
    """
    Manage parameter values by loading from a file, displaying a user input form for modification,
    and saving the updated values back to the file. This version also allows loading parameters
    from a given array.

    Parameters
    ----------
    param_list : list of str
        A list of parameter names to be managed and saved. This list is used both in the GUI form
        and for saving the names to the output file.
    name : str
        The base name for the text file where the parameters will be saved and loaded from.
        The file is expected to be named `{name}.txt`.
    load_directly : bool, optional
        If True, loads parameters directly from the file without displaying the GUI.
    parameter_list : list of float, optional
        If provided, this array of values will be used as the parameter values instead of loading from a file,
        and no form will be displayed.

    Returns
    -------
    fixed_values_list : list of float
        A list of the final parameter values, either loaded from the file, provided as an array,
        or updated through user input.
    """

    global fixed_values_list

    # Initialize the flag to ensure it is always defined
    should_print_after_modification = False

    # If parameter_list is provided, load values directly and skip the form
    if parameter_list is not None:
        fixed_values_list = parameter_list
        print("Loaded parameters from the provided array:")
        for param, value in zip(param_list, fixed_values_list):
            print(f"{param}: {value}")

        # Save the provided parameters directly to the file
        with open(f"{name}.txt", "w") as file:
            file.write(",".join(map(str, param_list)) + "\n")  # Write parameter names
            file.write(",".join(map(str, fixed_values_list)) + "\n")  # Write fixed values

        # Skip displaying the form
        return fixed_values_list

    # If parameter_list is not provided, proceed with loading from the file or showing the form
    file_exists = True
    try:
        with open(f"{name}.txt", "r") as file:
            lines = file.readlines()
            fixed_values_list = [float(val) for val in lines[1].strip().split(",")]
            if load_directly:
                print("Loaded parameters from file:")
                for param, value in zip(param_list, fixed_values_list):
                    print(f"{param}: {value}")
            else:
                should_print_after_modification = True
    except FileNotFoundError:
        file_exists = False
        fixed_values_list = [np.nan] * len(param_list)
        print("File not found. Initializing parameters to NaN:")
        for param in param_list:
            print(f"{param}: NaN")

    # If load_directly is False or the file does not exist, display the GUI
    if not load_directly or not file_exists:
        _load_values_to_interface(param_list, fixed_values_list, name)
        should_print_after_modification = True

    # Print the parameters after modification if needed
    if should_print_after_modification:
        print("Parameters after modification or user input:")
        for param, value in zip(param_list, fixed_values_list):
            print(f"{param}: {value}")

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
    
def _manage_bounds(param_list, name, load_directly=False, bound_list=None):
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


def _manage_bounds(param_list, name, load_directly=False, bound_list=None):
    """
    Manage parameter bounds by loading from a file, displaying a user input form for modification,
    and saving the updated bounds back to the file. This version also allows loading bounds from a
    provided array.

    Parameters
    ----------
    param_list : list of str
        A list of parameter names to be managed and saved. This list is used both in the GUI form
        and for saving the parameter names to the output file.
    name : str
        The base name for the bounds file where the parameters will be saved and loaded from.
        The file is expected to be named `bounds_{name}.txt`.
    load_directly : bool, optional
        If True, loads bounds directly from the file without displaying the GUI.
    bound_list : list of tuple, optional
        If provided, this array of tuples (lower, upper) will be used as the bounds instead of
        loading from a file, and no form will be displayed.

    Returns
    -------
    lower_bounds : list of float
        The final list of valid lower bounds.
    upper_bounds : list of float
        The final list of valid upper bounds.
    """

    global lower_bounds_list, upper_bounds_list

    # If bound_list is provided, load values directly and skip the form
    if bound_list is not None:
    
        lower_bounds_list, upper_bounds_list = bound_list
        print("Loaded bounds from the provided array:")
        for param, lower, upper in zip(param_list, lower_bounds_list, upper_bounds_list):
            print(f"{param}: Lower = {lower}, Upper = {upper}")

        # Save the provided bounds directly to the file
        with open(f"{name}.txt", "w") as file:
            file.write(",".join(map(str, param_list)) + "\n")  # Write parameter names
            file.write(",".join(map(str, lower_bounds_list)) + "\n")  # Write lower bounds
            file.write(",".join(map(str, upper_bounds_list)) + "\n")  # Write upper bounds


    if bound_list is None:
        # Initialize the flag to ensure it is always defined
        file_exists = True
        should_print_after_modification = False

        # Try to load previous bounds if available
        try:
            with open(f"{name}.txt", "r") as file:
                lines = file.readlines()
                lower_bounds_list = [float(val) for val in lines[1].strip().split(",")]
                upper_bounds_list = [float(val) for val in lines[2].strip().split(",")]
                if load_directly:
                    print("Loaded bounds from file:")
                    for param, lower, upper in zip(param_list, lower_bounds_list, upper_bounds_list):
                        print(f"{param}: Lower = {lower}, Upper = {upper}")
                else:
                    should_print_after_modification = True
        except FileNotFoundError:
            file_exists = False
            lower_bounds_list = [np.nan] * len(param_list)
            upper_bounds_list = [np.nan] * len(param_list)
            print("File not found. Initializing bounds to NaN:")
            for param in param_list:
                print(f"{param}: Lower = NaN, Upper = NaN")

        # If load_directly is False or the file does not exist, display the GUI
        if not load_directly or not file_exists:
            _load_bounds_to_interface(param_list, lower_bounds_list, upper_bounds_list, name)
            should_print_after_modification = True

        # Print the bounds after modification if needed
        if should_print_after_modification:
            print("Bounds after modification or user input:")
            for param, lower, upper in zip(param_list, lower_bounds_list, upper_bounds_list):
                print(f"{param}--- Lower = {lower}, Upper = {upper}")

        # Save updated bounds to the file
        with open(f"{name}.txt", "w") as file:
            file.write(",".join(map(str, param_list)) + "\n")  # Write parameter names
            file.write(",".join(map(str, lower_bounds_list)) + "\n")  # Write lower bounds
            file.write(",".join(map(str, upper_bounds_list)) + "\n")  # Write upper bounds

    return lower_bounds_list, upper_bounds_list
