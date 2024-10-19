
## FUNCTIONS

###Theoretical Functions

These functions will enable us to model the theoretical behavior of various astrophysical features influenced by orbital modulations with the aim of comparing our observational results with the expected.

####doppler_orbit_theoretical
-------------------------------------------------------------------------------
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
----------------------------------------------------------------------------------
####doppler_spiral_theoretical
----------------------------------------------------------------------------------
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
----------------------------------------------------------------------------------
####doppler_disc_theoretical
----------------------------------------------------------------------------------
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
----------------------------------------------------------------------------------
####doppler_spiral_in_orbit_theoretical
----------------------------------------------------------------------------------
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

----------------------------------------------------------------------------------
####density_through_orbit_theoretical
----------------------------------------------------------------------------------
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

----------------------------------------------------------------------------------
####absorption_column_through_orbit_theoretical
----------------------------------------------------------------------------------
    This function visualizes the column density (NH1, x 10^22 cm^-2) encountered by radiation emitted at each orbital phase as it travels towards an observer. It assumes a spherically distributed, neutral (unionized) stellar wind based on the CAK model.

    Parameters:
    - resolution (float, optional): Resolution for the phase array. Default is 0.01.
    - show_plot (bool, optional): If True, plots the absorption column density through the orbit. Default is False.
    - A form will appear to input the necessary orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - time (array-like): Time array.
    - phase (array-like): Orbital phase array.
    - NH1 (array-like): Absorption column density (NH1, x 10^22 cm^-2) through the orbit.
----------------------------------------------------------------------------------
####density_and_ionization_orbital_phase_theoretical
----------------------------------------------------------------------------------
    Calculates and visualizes the density and ionization parameter (log(ξ)) encountered by radiation emitted
    at each orbital phase as it travels towards an observer. Assumes a spherically distributed, neutral stellar
    wind based on the CAK (Castor, Abbott, Klein) model. The density profile and ionization parameter are calculated
    along the path from a neutron star (NS) through the stellar wind of its companion.

    Parameters:
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

    Returns:
    z : array-like
        Array of distances along the path from the neutron star to the observer (in cm).
    density : array-like
        Density profile of the stellar wind (in cm$^{-3}$) along the path.
    chi : array-like
        Ionization parameter (log(ξ)) calculated at each point along the path.
    
----------------------------------------------------------------------------------
####ionization_map_phase
----------------------------------------------------------------------------------
    Generates a logarithmic ionization parameter map based on the stellar wind density
    and the ionization parameters influenced by the source luminosity and distance to the compact object.

    Parameters:
    - size_in_Rstar (float, optional): Extent of the map from the stellar center in stellar radii. Default is 3.
    - min_color (float, optional): Minimum color scale value for the ionization parameter. Default is None.
    - max_color (float, optional): Maximum color scale value for the ionization parameter. Default is None.
    - save_plot (bool, optional): Boolean to decide if the plot should be saved. Default is False.
    - name (str, optional): Name of the file to save the plot. Default is "ionization_map".
    - A form will appear to input the necessary orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.
    

    Returns:
    - chi_result (pd.DataFrame): DataFrame containing the ionization parameter map.
    - area (float): The calculated area between bounds.
----------------------------------------------------------------------------------
####orbital_phase_to_time
----------------------------------------------------------------------------------
    Converts orbital phase array to time array for a compact object orbiting a companion star.
    The compact object moves faster at periastron than at apoastro. The increased orbital speed at periastron is primarily due to the conservation of angular momentum, which dictates that as the compact object moves closer to the central star, it must travel faster to maintain the total angular momentum of the system. This relationship is further influenced by Kepler’s laws of planetary motion, which describe how objects sweep out equal areas in equal times and the gravitational force between two bodies, which strengthens as they approach each other and weakens as they move apart.

    Parameters:
    - ph (array-like): Orbital phase array.
    - precision (float, optional): Resolution for the phase array. Default is 0.01.
    - A form will appear to input the necessary orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - ph (array-like): Orbital phase array (same as input).
    - time (array-like): Time array corresponding to the orbital phase.
    - W (array-like): Angular velocity array corresponding to the orbital phase.
    
----------------------------------------------------------------------------------
####orbital_time_to_phase
----------------------------------------------------------------------------------
    Converts orbital time array to phase array for a compact object orbiting a companion star.
    The compact object moves faster at periastron than at apoastro. The increased orbital speed at periastron is primarily due to the conservation of angular momentum, which dictates that as the compact object moves closer to the central star, it must travel faster to maintain the total angular momentum of the system. This relationship is further influenced by Kepler’s laws of planetary motion, which describe how objects sweep out equal areas in equal times and the gravitational force between two bodies, which strengthens as they approach each other and weakens as they move apart.

    Parameters:
    - ph (array-like): Orbital phase array.
    - precision (float, optional): Resolution for the phase array. Default is 0.01.
    - A form will appear to input the necessary orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.

    Returns:
    - ph (array-like): Orbital phase array (same as input).
    - time (array-like): Time array corresponding to the orbital phase.
    - W (array-like): Angular velocity array corresponding to the orbital phase.
    
    
    The increased orbital speed at periastron is primarily due to the conservation of angular momentum, which dictates that as the compact object moves closer to the central star, it must travel faster to maintain the total angular momentum of the system. This relationship is further influenced by Kepler’s laws of planetary motion, which describe how objects sweep out equal areas in equal times and the gravitational force between two bodies, which strengthens as they approach each other and weakens as they move apart.
----------------------------------------------------------------------------------
###Fitting 

If we have real data from a feature (NS spin period, emission line...) affected by Doppler or absorption or stellar density affected by orbital modulations we can recover the orbital parameters with the following functions. Keep in mind that "_ps" referes to particle swarm while "_ls" to least squares fitting approaches. "_ls" might not allways converge and "_ps" might take a lot of time to complete.

####fit_orbit_ps
----------------------------------------------------------------------------------
    fit_orbit_ps fits observed orbital modulation data by estimating parameters such as phase, semi-major axis, orbital period, eccentricity, inclination, and periapsis.

    The fitting process utilizes a particle swarm optimization (PSO) algorithm, which iteratively improves parameter estimates by minimizing the chi-squared difference between observed and predicted data.

    The function can handle two fitting methods:

    - Discrete: Suitable for discrete data points (e.g., spectra with small orbital phase ranges, faster).
    - Extended: Suitable for data with varying or extended bin sizes, typical with current instrument resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Inputs:
    - x_data: Time bins of the observed data.
    - y_data: Observed data points corresponding to each time bin.
    - y_err: Error associated with each observed data point (default is 0).
    - num_iterations: Number of iterations for PSO optimization (default is 3).
    - maxiter: Maximum number of iterations for each PSO run (default is 1000).
    - swarmsize: Number of particles in the PSO swarm (default is 100).
    - units: Units of the observed data (default is "keV").
    - method_: Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize: Bin size for extended method (default is 0.01).
    - A form will appear to input the necessary bounds for the orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.


    Outputs:
    - DataFrame of the best-fit parameters and their standard deviations.
    - ph: Array of phases corresponding to the predicted data.
    - predicted_data: Predicted data based on the best-fit parameters.
    - r_squared: Coefficient of determination (R-squared) indicating the quality of fit.
    
----------------------------------------------------------------------------------
####fit_orbit_ls
----------------------------------------------------------------------------------

    fit_orbit_ls fits orbital modulation data by estimating parameters such as phase, semi-major axis, orbital period, eccentricity, inclination, and periapsis.

    The fitting process utilizes a traditional Least Squares (LS) method, which fits the observed data by minimizing the squared differences between observed and predicted values.

    The function can handle two fitting methods:

    - Discrete: Suitable for discrete data points (e.g., spectra with small orbital phase ranges, faster).
    - Extended: Suitable for data with varying or extended bin sizes, typical with current instrument resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Inputs:
    - x_data: Time bins of the observed data.
    - y_data: Observed data points corresponding to each time bin.
    - y_err: Error associated with each observed data point (default is 0).
    - units: Units of the observed data (default is "keV").
    - method_: Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize: Bin size for extended method (default is 0.01).
    - A form will appear to input the necessary bounds for the orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.


    Outputs:
    - DataFrame of the best-fit parameters and their errors.
    - ph: Array of phases corresponding to the predicted data.
    - predicted_data: Predicted data based on the best-fit parameters.
    - r_squared: Coefficient of determination (R-squared) indicating the quality of fit.

----------------------------------------------------------------------------------
####fit_disc_ps
----------------------------------------------------------------------------------
    fit_disc_ps fits orbital modulation data by estimating parameters such as phase, semi-major axis, orbital period, eccentricity, inclination for the main orbit, and corresponding parameters for a secondary orbit (e.g., ballistic capture of matter around a compact object or an accretion disk).

    The fitting process uses a particle swarm optimization (PSO) algorithm, which iteratively improves parameter estimates by minimizing the chi-squared difference between observed and predicted data.

    The function can handle two fitting methods:

    - Discrete: Suitable for discrete data points (e.g., spectra with small orbital phase ranges, faster).
    - Extended: Suitable for data with varying or extended bin sizes, typical with current instrument resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Inputs:
    - x_data: Time bins of the observed data.
    - y_data: Observed data points corresponding to each time bin.
    - y_err: Error associated with each observed data point (default is 0).
    - num_iterations: Number of iterations for the PSO algorithm (default is 3).
    - maxiter: Maximum number of iterations for each PSO run (default is 1000).
    - swarmsize: Number of particles in the PSO swarm (default is 100).
    - units: Units of the observed data (default is "keV").
    - method_: Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize: Bin size for extended method (default is 0.01).
    - A form will appear to input the necessary bounds for the orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.


    Outputs:
    - DataFrame of the best-fit parameters and their standard deviations.
    - ph: Array of phases corresponding to the predicted data for the main orbit.
    - ph2: Array of phases corresponding to the predicted data for the secondary orbit.
    - predicted_data: Predicted data based on the best-fit parameters.
    - r_squared: Coefficient of determination (R-squared) indicating the quality of fit.

----------------------------------------------------------------------------------
####fit_disc_ls
----------------------------------------------------------------------------------
    fit_disc_ls fits orbital modulation data by estimating parameters such as phase, semi-major axis, orbital period, eccentricity, inclination for the main orbit, and corresponding parameters for a secondary orbit (e.g., ballistic capture of matter around a compact object or an accretion disk).

    The fitting process uses a traditional least squares (LS) method, provided here for completeness despite its potential limitations due to the model's complexity.

    The function can handle two fitting methods:

    - Discrete: Suitable for discrete data points (e.g., spectra with small orbital phase ranges, faster).
    - Extended: Suitable for data with varying or extended bin sizes, typical with current instrument resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Inputs:
    - x_data: Time bins of the observed data.
    - y_data: Observed data points corresponding to each time bin.
    - y_err: Error associated with each observed data point (default is 0).
    - units: Units of the observed data (default is "keV").
    - method_: Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize: Bin size for extended method (default is 0.01).
    - A form will appear to input the necessary bounds for the orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.


    Outputs:
    - DataFrame of the best-fit parameters and their errors.
    - ph: Array of phases corresponding to the predicted data.
    - predicted_data: Predicted data based on the best-fit parameters.
    - r_squared: Coefficient of determination (R-squared) indicating the quality of fit.

----------------------------------------------------------------------------------
####fit_spiral_ps
----------------------------------------------------------------------------------
    fit_spiral_ps fits orbital modulation data by estimating parameters for a spiral orbit.

    The fitting process uses a particle swarm optimization (PSO) algorithm, which iteratively improves the parameter estimates by minimizing the chi-squared difference between the observed and predicted data.

    The function can handle two fitting methods:

    - Discrete: Suitable for discrete data points (e.g., spectra with small orbital phase ranges, faster).
    - Extended: Suitable for data with varying or extended bin sizes, typical with current instrument resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Inputs:
    - x_data: Time bins of the observed data.
    - y_data: Observed data points corresponding to each time bin.
    - y_err: Error associated with each observed data point (default is 0).
    - num_iterations: Number of PSO iterations to perform (default is 3).
    - maxiter: Maximum number of iterations for PSO (default is 1000).
    - swarmsize: Number of particles in the swarm for PSO (default is 100).
    - units: Units of the observed data (default is "keV").
    - method_: Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize: Bin size for extended method (default is 0.01).
    - A form will appear to input the necessary bounds for the orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.


    Outputs:
    - DataFrame of the best-fit parameters and their standard deviations.
    - ph: Array of phases corresponding to the predicted data.
    - predicted_data: Predicted data based on the best-fit parameters.
    - r_squared: Coefficient of determination (R-squared) indicating the quality of fit.

----------------------------------------------------------------------------------
####fit_spiral_ls
----------------------------------------------------------------------------------
    fit_spiral_ls fits orbital modulation data by estimating parameters for a spiral orbit.

    The fitting process uses a traditional least squares (LS) method, provided for completeness due to the complexity of the model.

    The function can handle two fitting methods:

    - Discrete: Suitable for discrete data points (e.g., spectra with small orbital phase ranges, faster).
    - Extended: Suitable for data with varying or extended bin sizes, typical with current instrument resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Inputs:
    - x_data: Time bins of the observed data.
    - y_data: Observed data points corresponding to each time bin.
    - y_err: Error associated with each observed data point (default is 0).
    - units: Units of the observed data (default is "keV").
    - method_: Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize: Bin size for extended method (default is 0.01).
    - A form will appear to input the necessary bounds for the orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.


    Outputs:
    - DataFrame of the best-fit parameters and their errors.
    - ph: Array of phases corresponding to the predicted data.
    - predicted_data: Predicted data based on the best-fit parameters.
    - r_squared: Coefficient of determination (R-squared) indicating the quality of fit.

----------------------------------------------------------------------------------
####fit_spiral_in_orbit_ps
----------------------------------------------------------------------------------
    fit_spiral_in_orbit_ps fits orbital modulation data by estimating parameters for a spiral orbit contained in a main orbit.

    The fitting process uses a particle swarm optimization (PSO) algorithm, which iteratively improves the parameter estimates by minimizing the chi-squared difference between the observed and predicted data.

    The function can handle two fitting methods:

    - Discrete: Suitable for discrete data points (e.g., spectra with small orbital phase ranges, faster).
    - Extended: Suitable for data with varying or extended bin sizes, typical with current instrument resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Inputs:
    - x_data: Time bins of the observed data.
    - y_data: Observed data points corresponding to each time bin.
    - y_err: Error associated with each observed data point (default is 0).
    - units: Units of the observed data (default is "keV").
    - method_: Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize: Bin size for extended method (default is 0.01).
    - A form will appear to input the necessary bounds for the orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.


    Outputs:
    - DataFrame of the best-fit parameters and their errors.
    - ph: Array of phases corresponding to the predicted data.
    - predicted_data: Predicted data based on the best-fit parameters.
    - r_squared: Coefficient of determination (R-squared) indicating the quality of fit.

----------------------------------------------------------------------------------
####fit_spiral_in_orbit_ls
----------------------------------------------------------------------------------
    fit_spiral_in_orbit_ps fits orbital modulation data by estimating parameters for a spiral orbit contained in a main orbit.

    The fitting process uses a traditional least squares (LS) method, provided for completeness due to the complexity of the model.

    The function can handle two fitting methods:

    - Discrete: Suitable for discrete data points (e.g., spectra with small orbital phase ranges, faster).
    - Extended: Suitable for data with varying or extended bin sizes, typical with current instrument resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Inputs:
    - x_data: Time bins of the observed data.
    - y_data: Observed data points corresponding to each time bin.
    - y_err: Error associated with each observed data point (default is 0).
    - units: Units of the observed data (default is "keV").
    - method_: Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize: Bin size for extended method (default is 0.01).
    - A form will appear to input the necessary bounds for the orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.


    Outputs:
    - DataFrame of the best-fit parameters and their errors.
    - ph: Array of phases corresponding to the predicted data.
    - predicted_data: Predicted data based on the best-fit parameters.
    - r_squared: Coefficient of determination (R-squared) indicating the quality of fit.
----------------------------------------------------------------------------------
####fit_nh_ps
----------------------------------------------------------------------------------
    fit_nh_ps fits the column density (NH1, x 10^22 cm^-2) encountered by radiation emitted at each orbital phase as it travels towards an observer. It assumes a spherically distributed, neutral (unionized) stellar wind based on the CAK model.

    The fitting process uses a particle swarm optimization (PSO) algorithm to minimize the chi-squared difference between observed and predicted data points.

    The function supports two fitting methods:

    - Discrete: Suitable for discrete data points (e.g., spectra with small orbital phase ranges, faster).
    - Extended: Suitable for data with varying or extended bin sizes, typical for current instrument resolutions (e.g., XMM-Newton and Chandra) and short X-ray binary (XRB) orbits.

    Inputs:
    - x_data: Time bins of the observed data.
    - y_data: Observed data points corresponding to each time bin.
    - y_err: Error associated with each observed data point (default is 0).
    - num_iterations: Number of iterations for the PSO algorithm (default is 3).
    - maxiter: Maximum number of iterations for each PSO run (default is 200).
    - swarmsize: Number of particles in the PSO swarm (default is 20).
    - method_: Fitting method to use, either "discrete" or "extended" (default is "extended").
    - extended_binsize: Bin size for the extended method (default is 0.01).
    - A form will appear to input the necessary bounds for the orbital parameters. These parameters will be saved in a .txt file in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter all parameters; you can modify only those that require adjustment.


    Outputs:
    - df_results_transposed: Transposed DataFrame of the best-fit parameters and their standard deviations.
    - ph: Array of phases corresponding to the predicted data.
    - predicted_data: Predicted data based on the best-fit parameters.
    - r_squared: Coefficient of determination (R-squared) indicating the quality of fit.

----------------------------------------------------------------------------------
### Timing related 

In this section we collect some usefull timing functions. From them we can obtain NS pulses in a sliding window approach allowing us to observe if Doppler modulations are present and recover the orbital parameters using fitting functions.

####hr
----------------------------------------------------------------------------------
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
----------------------------------------------------------------------------------
####cr
----------------------------------------------------------------------------------
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
----------------------------------------------------------------------------------
#### rebin_snr 
----------------------------------------------------------------------------------
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
----------------------------------------------------------------------------------
#### rebin_bins
----------------------------------------------------------------------------------
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
----------------------------------------------------------------------------------
####  fold_pulse
----------------------------------------------------------------------------------
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
----------------------------------------------------------------------------------
####   period_sliding_window
----------------------------------------------------------------------------------
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
----------------------------------------------------------------------------------

