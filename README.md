
# XRAYBINARYORBIT: A Python Package for analyzing X-ray binaries orbital modulations


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

- Python >= 3.7
- NumPy >= 1.19
- SciPy >= 1.5
- Matplotlib >= 3.2 (optional, for plotting)
- Astropy >= 4.0 (for astronomical time and coordinate handling)
- PySwarms >= 1.1 (for particle swarm optimization)
- Tinker 
---

## Usage

### For all Functions:
In this code, a user-friendly form is used to handle the various parameters that influence orbital modulations. The user inputs are saved in a file in the running directory, which is automatically loaded during future runs to avoid re-entering parameters. If the file doesn’t exist, the form will load for new inputs. However, if "load_directly=True" and the file exists, the code will run using the saved parameters.


### For Fitting Functions:
For fitting the orbital parameters, two approaches are available: least squares (LS) and particle swarm (PSO), named *_ls and *_ps respectively. The least squares method is faster but doesn’t always converge, while the particle swarm method is more reliable but computationally expensive. Key parameters for particle swarm include num_iterations, maxiter, and swarmsize. It’s recommended to start with smaller values for these parameters (e.g., num_iterations=3, maxiter=250, and swarmsize=50) to gauge the computational demand and adjust as needed.

For each of these fiting approaches we have two different methods: extended and discrete. In order to fit our data to some orbital modulation we have a list of "values" corresponding to a list of "time sections" (the duration of our phase resolved spectra, lightcurve section etc...). 
When we are triying to fit an orbital modulation with for example, a short period, taking into account that orbital modulations are in general sinusoidal, is not possible to just use the center of the "time section" but we have to consider the average of the orbital modulation during the extent of the "time sections" and fit that to our "values".

As this highly complicates the calculus, the parameter "extended_binsize" gives us a value from which we will consider the point to be discrete or extended (i.e. if extended_binsize=0.01, if the size of our time section covers more than 0.01 orbital phases, it will be trated as extended, and if less as discrete). 
To use the extended approach we can provide our list of times as pairs (see Fe XXV Doppler Shifts Cen X-3 example) or we can provide a list of "times", one component longer than the "values" list and the pairs will be created authomatically. If the list of "times" is the same length of the list of "values" the discrete approach will be used.

---

## Documentation

Detailed documentation, including function references and tutorials, can be found at: [XRAYBINARYORBIT Documentation](https://xragua.github.io/xraybinaryorbit/).

---

## Contributing

Contributions are welcome! Please check the [CONTRIBUTING.md](link-to-contributing) for details on how to contribute.

---

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/xragua/xraybinaryorbit/blob/main/LICENSE) file for details.

---
