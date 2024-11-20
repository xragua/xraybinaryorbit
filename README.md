![Codecov](https://img.shields.io/codecov/c/github/xragua/xraybinaryorbit)

[![status](https://joss.theoj.org/papers/7f4fbc6a63adc4ba5a07bb340a4d7246/status.svg)](https://joss.theoj.org/papers/7f4fbc6a63adc4ba5a07bb340a4d7246)


# XRAYBINARYORBIT: A Python Package for Analyzing Orbital Modulations in X-ray Binaries


**X-ray binaries are truly fascinating!** In these extreme environments, a compact object—either a neutron star or black hole—draws in matter from a companion star, producing intense X-ray emissions. These systems offer a unique window into extreme physics, from the effects of strong gravity and relativistic jets to the presence of intense magnetic fields.

Orbital modulations are observed in nearly all X-ray binary systems. These variations arise from the orbital motions of the system, driven by the relative velocities of the two stars and their changing configurations with respect to each other and the observer.

We can observe how the center of an emission line slightly changes during a phase resolved analysis, or how the NS spin period slightly varies following a trend. These phenomena can be caused by Doppler effect. Our code will help you turn these observations into the orbital parameters that cause those Doppler shifts.

But this simple idea can get tricky when you consider all the factors involved. Inclination, eccentricity, periapsis, distance to the barycenter (which depends on the masses of the stars involved), and orbital phases really matter in this analysis. Plus, if there’s eccentricity different than 0, the velocity around the orbit isn’t constant. 

If we think about stellar wind (the matter accreted by our compact objects) there are many combinations that can result into orbital modulations. With the eccentricity the density around the orbit changes, and thus, the accreted matter. With eccentricity and inclination the absorption column faced by the emitted radiation varies depending on the orbital phase, so does the ionization of the wind.

These orbital modulations are easy to grasp but not so straightforward to analyze—yet that’s where our tools come in.

The primary challenge in this type of analysis has long been the lack of sufficient resolution for detailed phase-resolved observations. However, upcoming high-resolution missions, such as **XRISM** and **New Athena**, promise to take these analyses to the next level. But it’s not just about improved resolution—advances in computational power are equally crucial. Many of these tools have already been successfully applied in studies using **XMM-Newton** and **Chandra** data, enabling analyses that would have been impossible just a few years ago.


Our code comprises three different king of functions: 

- **Theoretical: To predice which modulations should be expected.
- **Fitting: To retrieve orbital parameters and stellar wind characteristics from our observed data.
- **Timming: Some simple functions to manage lightcurves, calculating periods and foldin pulses. 


---

## Installation

You can install the package directly from PyPI using pip: **pip install xraybinaryorbit**.

Or download the code from https://github.com/xragua/xraybinaryorbit/releases/tag/0.2.9.

Some examples of their usage are presented https://github.com/xragua/xraybinaryorbit/tree/main/examples.

---

## Dependencies

The following Python libraries are required:

-- numpy
-- pandas
-- scipy
-- pyswarm
-- matplotlib
-- tk
-- astropy
-- pytest>=6.0
-- xraybinaryorbit

---

## Usage

### General Usage

The software provides a user-friendly interface for managing the various parameters that influence orbital modulations. Upon first use, the user inputs parameters through a form, which are then saved in a file within the running directory. This file is automatically loaded in future runs, eliminating the need to re-enter all parameters.
 If the file is absent, the form will reappear for new inputs. Alternatively, setting `load_directly=True` will bypass the form and run the code using previously saved parameters (only if the file exists).
Alternatively, the input parameters or bounds can be provided as lists, by providing a "parameter_list" or "bound_list" as imputs. This values will be saved for future runs within the working directory as a .txt for subsequent runs.

### Fitting Functions

For fitting orbital parameters, the software offers two approaches: least squares (LS) and particle swarm optimization (PSO) denoted by `_ls` and `_ps`, within the functions name respectively. The LS methods is faster but may fail to converge in certain cases, whereas the PS0 is more robust but computationally intensive. Key parameters for PSO include `num_iterations`, `maxiter`, and `swarmsize`. It is recommended to start with smaller values for these parameters (e.g., `num_iterations=3`, `maxiter=100`, and `swarmsize=20`) to evaluate computational demand and adjust accordingly. 

### Extended vs. Discrete Fitting

The software provides two fitting methods: `extended` and `discrete`. When fitting data to a certain orbital modulation, a list of `values` (either energies, pulses or wave lengths) corresponding to `time sections` (e.g., phase-resolved spectra or lightcurve segments, wich have a defined time lenght) will be our data to fit. For short-period orbital modulations, which are often sinusoidal, the center of the time section may not properly represent the modulation (that would be a discrete approach, i.e. an discrete time/orbital phase for a discrete value). Instead, it’s necessary to consider the average modulation across the entire `time section` and fit that to the observed `values`.

To handle this complexity, the `extended_binsize` parameter is used. If the size of the time section covers more than `extended_binsize` orbital phases (e.g., `extended_binsize=0.01`), the point is treated as extended; otherwise, it is treated as discrete authomatically, within the `extended` mode, which is the default. 

The `extended` approach can be utilized by providing a list of `time pairs`, or alternatively, a single list of `times` with one extra element compared to the `values` list, from which `time pairs` will be automatically created. If the list of `times` is the same length as the `values` list, the `discrete` approach will be used instead, even if the `extended` method was selected (a warnig will appear indicating it so).

---

## Documentation

Detailed documentation, including function references and tutorials, can be found at: [XRAYBINARYORBIT Documentation](https://xragua.github.io/xraybinaryorbit/).

---

## Contributing

Contributions are welcome! Please check the [CONTRIBUTING FILE](https://github.com/xragua/xraybinaryorbit/blob/main/contributing.md) for details on how to contribute.

---

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/xragua/xraybinaryorbit/blob/main/LICENSE) file for details.

---
