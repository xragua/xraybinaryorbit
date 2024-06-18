# xraybinaryorb

xraybinaryorb is a user-friendly Python package designed to facilitate the analysis of orbital modulations in X-ray binaries. Orbital modulations in these systems are always present, potentially providing valuable information about XRB systems, but due to the lack of resolution, they are generally overlooked. These modulations can be complex, with numerous parameters influencing the observed X-ray flux.

### Features
To address these challenges, xraybinaryorb offers an input system based on interactive forms. These forms, once loaded, are saved for future interactions. Each function is briefly described when loaded, helping to guide its execution.

With the advent of new, high-resolution telescopes like Athena, the importance of accurately analyzing these modulations has increased significantly. xraybinaryorb provides a suite of functions that simplify this process, making it more accessible to researchers. Whether you're dealing with eclipses, Doppler shifts, or absorption variations, xraybinaryorb offers the tools you need to effectively interpret and understand the intricate behaviors of X-ray binaries.

The primary physics used for these packages rely on the Doppler effect, Kepler's laws, and the CAK model (Castor, Abbott, and Klein) as described in:

Castor, J. I., Abbott, D. C., & Klein, R. I. (1975). Radiation-driven winds in Of stars. *Astrophysical Journal*, 195, 157-174.

### Available Functions

#### Theoretical Functions
- **doppler_orbit_theoretical**: Calculates the Doppler effect in the orbital motion.
- **doppler_spiral_theoretical**: Models the Doppler effect in spiral structures.
- **doppler_disc_theoretical**: Models the Doppler effect in accretion discs.
- **doppler_spiral_in_orbit_theoretical**: Combines orbital and spiral Doppler effects.
- **density_through_orbit_theoretical**: Computes wind density through the orbit.
- **absorption_column_through_orbit_theoretical**: Calculates absorption column variations through the orbit.
- **ionization_map_phase**: Maps ionization levels across orbital phases.
- **orbital_phase_to_time**: Converts orbital phase to time.
- **orbital_time_to_phase**: Converts orbital time to phase.

#### Fitting Functions
(Note: 'ps' stands for particle swarm optimization and 'ls' stands for least squares. Particle swarm optimization is preferred as least squares does not always converge.)
- **fit_orbit_ps**: Fits orbital parameters using particle swarm optimization.
- **fit_orbit_ls**: Fits orbital parameters using least squares.
- **fit_disc_ps**: Fits disc parameters using particle swarm optimization.
- **fit_disc_ls**: Fits disc parameters using least squares.
- **fit_spiral_ps**: Fits spiral parameters using particle swarm optimization.
- **fit_spiral_ls**: Fits spiral parameters using least squares.
- **fit_spiral_in_orbit_ps**: Fits combined spiral and orbit parameters using particle swarm optimization.
- **fit_spiral_in_orbit_ls**: Fits combined spiral and orbit parameters using least squares.
- **nh_orbit**: Calculates the column density (N_H) variations through the orbit.
- **fit_nh_ps**: Fits N_H variations using particle swarm optimization.

#### Timing Functions
- **hr**: Calculates hardness ratios.
- **cr**: Computes count rates.
- **rebin_snr**: Rebins data to achieve a specific signal-to-noise ratio.
- **rebin_bins**: Rebins data into a specified number of bins.
- **fold_pulse**: Folds pulse profiles over the orbital period.
- **period_sliding_window**:(based in Lomb Scargle periodogram) Analyzes period changes using a sliding window method.

### Examples
Some examples require additional files, which are contained in each folder. If used, all files should be kept in the running directory.

#### Examples based on theoretical data:
- **example_doppler_theoretical**: Demonstrates Doppler effect calculations in theoretical models.
- **example_density_NH**: Shows how to calculate density and column density (N_H) variations.
- **example_time_phase_relation**: Illustrates the conversion between orbital time and phase.

#### Examples based on real data:
- **example_ionization_map**: Provides an example of mapping ionization levels using real observational data.
- **example_NS_doppler_shifts**: Demonstrates the analysis of neutron star Doppler shifts using real data.
