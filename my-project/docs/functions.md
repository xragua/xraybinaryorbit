# Function reference

This page documents the current `xraybinaryorbit` development API. The scientific equations, assumptions, and interpretation of the optimization methods are described in [Scientific background and numerical methods](science.md).

!!! note "Development documentation"
    This reference follows the corrected development modules reviewed in July 2026. When using a released PyPI version, check its version number because older releases may retain different conventions, defaults, or return values.

!!! warning "Theoretical Doppler wrappers"
    `doppler_orbit_theoretical` is aligned with the corrected eccentric-orbit kernel. The currently inspected `doppler_spiral_theoretical`, `doppler_disc_theoretical`, and `doppler_spiral_in_orbit_theoretical` wrappers still retain parts of the historical velocity implementation. Their signatures and current return values are documented here, but their internal physics should be synchronized with the corrected fitting kernels before the next release.

---

## Common conventions

### Time input

Fitting functions support two evaluation modes:

- `method_="discrete"`: `x_data` is a one-dimensional array containing one time per measurement.
- `method_="extended"`: `x_data` has shape `(N, 2)` and contains the start and end time of every measurement bin. A one-dimensional array of consecutive bin edges can also be converted internally into pairs.

Times are expressed in seconds. `extended_binsize` is a threshold in **orbital-phase cycles**. Bins narrower than this value are evaluated at their midpoint; wider bins are sampled and averaged.

### Doppler units

The canonical unit names are:

```python
units="keV"
units="s"
units="angstrom"
```

`feature` is the unshifted energy, period, or wavelength in the selected units.

!!! warning "Historical spelling"
    Older Doppler code used the misspelling `"amstrong"`. The corrected fitting kernels use `"angstrom"` and may retain the old spelling only as a compatibility alias. The remaining theoretical Doppler wrappers should be checked before release so that they all follow the same canonical spelling.

### Parameter and bound input

The theoretical functions support `load_directly` and `parameter_list`:

- `parameter_list=[...]` supplies parameters programmatically in the documented order.
- Without `parameter_list`, the parameter manager can open the interactive form.
- `load_directly=True` loads previously stored values without reopening the form, when those values exist.

Fitting functions similarly support `load_directly` and `bound_list`. The recommended explicit structure is

```python
bound_list = [lower_bounds, upper_bounds]
```

where both arrays follow the parameter order documented for that fit.

### Result tables

All fitting functions return a transposed `pandas.DataFrame` whose columns are parameter names and whose rows are:

- `Value`: selected best-fit value;
- `Std`: uncertainty or stability estimate.

The meaning of `Std` depends on the optimizer:

- `_ps`: standard deviation of the best parameter values across independent PSO runs;
- `_ls`: square root of the corresponding diagonal element of the local covariance matrix.

These are not equivalent uncertainty definitions.

---

# Theoretical functions

## `doppler_orbit_theoretical`

```python
doppler_orbit_theoretical(
    t,
    units="keV",
    show_plot=False,
    precision_for_phase=0.01,
    load_directly=False,
    parameter_list=None,
)
```

Calculate the Doppler modulation produced by the barycentric motion of an emitting compact object in an eccentric binary orbit. The corrected velocity includes both radial and tangential orbital terms.

**Parameter order**

```text
iphase, semimajor, orbitalperiod, eccentricity, periapsis,
inclination, Rstar, Mstar1, Mstar2, wind_vel, feature
```

Current conventions:

- `semimajor`: relative semi-major axis in donor-radius units;
- `Mstar1`: emitting compact object;
- `Mstar2`: donor or companion;
- `wind_vel`: optional projected phenomenological flow velocity in km s^-1;
- `feature`: rest energy, period, or wavelength.

**Returns**

```python
time, orbital_phase, shifted_feature
```

---

## `doppler_spiral_theoretical`

```python
doppler_spiral_theoretical(
    t,
    units="keV",
    show_plot=False,
    load_directly=False,
    parameter_list=None,
    verbose_complete=False,
)
```

Calculate the Doppler modulation associated with an isolated logarithmic spiral.

**Parameter order**

```text
iphase_spiral, semimajor_spiral, b, omega,
inclination_spiral, feature
```

`omega` is interpreted as cycles per second. `semimajor_spiral` is the initial spiral scale in solar radii in this standalone theoretical function.

**Returns**

```python
time, spiral_phase, shifted_feature, spiral_radius
```

---

## `doppler_disc_theoretical`

```python
doppler_disc_theoretical(
    t,
    units="keV",
    show_plot=False,
    load_directly=False,
    parameter_list=None,
)
```

Calculate a two-orbit Doppler model intended to represent an emitting region moving in an inner orbit while the inner system follows an outer binary orbit.

**Parameter order**

```text
iphase, semimajor, orbitalperiod, eccentricity, periapsis,
inclination, Rstar, Mstar1, Mstar2,
iphase2, semimajor2, orbitalperiod2, eccentricity2,
periapsis2, inclination2, Mass3, wind_vel, feature
```

Recommended physical interpretation:

- `Mstar1`: compact object at the centre of the inner orbit;
- `Mstar2`: donor or outer companion;
- `Mass3`: body or emitting disc element carrying the feature;
- `semimajor`: relative outer semi-major axis;
- `semimajor2`: relative inner semi-major axis.

**Returns**

```python
time, outer_phase, inner_phase, shifted_feature
```

---

## `doppler_spiral_in_orbit_theoretical`

```python
doppler_spiral_in_orbit_theoretical(
    t,
    units="keV",
    show_plot=False,
    load_directly=False,
    parameter_list=None,
)
```

Calculate the combined Doppler modulation from an eccentric binary orbit and a logarithmic spiral component.

**Parameter order**

```text
iphase, semimajor, orbitalperiod, eccentricity, periapsis,
inclination, iphase_spiral, semimajor_spiral, b, omega,
inclination_spiral, Rstar, Mstar1, Mstar2, wind_vel, feature
```

**Returns**

```python
time, orbital_phase, spiral_phase, shifted_feature
```

---

## `density_through_orbit_theoretical`

```python
density_through_orbit_theoretical(
    resolution=0.01,
    show_plot=False,
    load_directly=False,
    parameter_list=None,
)
```

Calculate the smooth-wind mass density at the compact object's relative orbital position.

**Parameter order**

```text
semimajor, orbitalperiod, eccentricity, periapsis, Rstar,
wind_infinite_velocity, Mass_loss_rate, beta
```

**Returns**

```python
time, orbital_phase, mass_density
```

`mass_density` is in g cm^-3.

---

## `absorption_column_through_orbit_theoretical`

```python
absorption_column_through_orbit_theoretical(
    resolution=0.01,
    show_plot=True,
    load_directly=False,
    parameter_list=None,
)
```

Calculate the equivalent hydrogen column between the compact object and the observer through a smooth beta-law wind.

**Parameter order**

```text
semimajor, orbitalperiod, eccentricity, periapsis, inclination,
Rstar, hydrogen_mass_fraction, wind_infinite_velocity,
Mass_loss_rate, beta
```

**Returns**

```python
time, orbital_phase, NH
```

`NH` is in units of 10^22 cm^-2. Under the current package convention, sightlines that intersect the stellar photosphere return `NH=0` as an eclipse flag.

---

## `density_and_ionization_orbital_phase_theoretical`

```python
density_and_ionization_orbital_phase_theoretical(
    resolution=0.01,
    size=10,
    show_plot=True,
    load_directly=False,
    parameter_list=None,
)
```

Calculate hydrogen number density, ionization parameter, and remaining column along a ray from the compact object towards the observer at one selected orbital phase.

**Parameter order**

```text
orb_phase, luminosity, semimajor, eccentricity, periapsis,
inclination, Rstar, hydrogen_mass_fraction,
wind_infinite_velocity, Mass_loss_rate, beta
```

`luminosity` retains the historical scale of 10^32 erg s^-1.

**Returns**

```python
distance_in_Rstar, n_H, log10_xi, NH_remaining
```

- `n_H`: cm^-3;
- `log10_xi`: logarithm of the ionization parameter in erg cm s^-1;
- `NH_remaining`: 10^22 cm^-2.

---

## `ionization_map_phase`

```python
ionization_map_phase(
    size_in_Rstar=0,
    min_color=None,
    max_color=None,
    save_plot=False,
    name="ionization_map",
    load_directly=False,
    parameter_list=None,
)
```

Generate a two-dimensional polar map of the ionization parameter around the donor, including the geometrical stellar shadow.

**Parameter order**

```text
phase, semimajor, eccentricity, periapsis, Rstar,
hydrogen_mass_fraction, wind_infinite_velocity, Mass_loss_rate,
beta, luminosity, bound1, bound2
```

Here `luminosity` is supplied directly in erg s^-1.

**Returns**

```python
xi_map, electron_density_map, area_between_bounds
```

The first two outputs are `pandas.DataFrame` objects. The area is returned in cm^2. If `size_in_Rstar=0`, the map extent is selected automatically from the orbit.

---

## `orbital_phase_to_time`

```python
orbital_phase_to_time(
    ph,
    precision=0.01,
    load_directly=False,
    parameter_list=None,
)
```

Convert orbital phase to time using Kepler's second law.

**Parameter order**

```text
iphase, orbitalperiod, eccentricity, periapsis
```

**Returns**

```python
phase, time, angular_velocity
```

`time` is in seconds and `angular_velocity` is in rad s^-1.

---

## `orbital_time_to_phase`

```python
orbital_time_to_phase(
    t,
    precision=0.01,
    load_directly=False,
    parameter_list=None,
)
```

Convert time to orbital phase using the inverse swept-area relation.

The parameter order and returned quantities are the same as for `orbital_phase_to_time`; the second returned array is the original time input.

---

# Fitting functions

## Common fitting arguments

- `x_data`: times or time bins.
- `y_data`: observed values.
- `y_err`: one-sigma errors. When omitted, the package constructs working weights; the absolute chi-square should then be interpreted cautiously.
- `units`: output units for Doppler models.
- `method_`: `"discrete"` or `"extended"`.
- `extended_binsize`: midpoint/averaging threshold in orbital-phase cycles.
- `bound_list`: lower and upper parameter bounds.
- `load_directly`: load saved bounds through the parameter interface.
- `num_iterations`: number of complete independent PSO runs.
- `swarmsize`: number of particles in each PSO run.
- `maxiter`: maximum internal swarm updates per run.

### PSO outputs

PSO functions select the independent run with the lowest weighted chi-square and return that run's parameters. Their `Std` row is the run-to-run standard deviation of the best parameter values.

### Least-squares outputs

Least-squares functions use bounded `scipy.optimize.curve_fit`. When non-zero measurement errors are supplied, they are passed as absolute `sigma` values. The initial point is the midpoint of the bounds. Their `Std` row is calculated from the local covariance matrix.

---

## `fit_orbit_ps`

```python
fit_orbit_ps(
    x_data,
    y_data,
    y_err=0,
    num_iterations=3,
    maxiter=1000,
    swarmsize=100,
    units="keV",
    method_="extended",
    extended_binsize=0.01,
    load_directly=False,
    bound_list=None,
)
```

Fit an eccentric binary Doppler model with Particle Swarm Optimization.

**Bound order**

```text
iphase, semimajor, orbitalperiod, eccentricity, periapsis,
inclination, Rstar, Mstar1, Mstar2, wind_vel, feature
```

**Returns**

```python
results, orbital_phase, predicted_data, chi_squared
```

---

## `fit_orbit_ls`

```python
fit_orbit_ls(
    x_data,
    y_data,
    y_err=0,
    units="keV",
    method_="extended",
    extended_binsize=0.01,
    load_directly=False,
    bound_list=None,
)
```

Fit the same eccentric binary Doppler model with bounded non-linear least squares.

The bound order is identical to `fit_orbit_ps`.

**Returns**

```python
results, orbital_phase, predicted_data, r_squared
```

---

## `fit_disc_ps`

```python
fit_disc_ps(
    x_data,
    y_data,
    y_err=0,
    num_iterations=3,
    maxiter=1000,
    swarmsize=100,
    units="keV",
    method_="extended",
    extended_binsize=0.01,
    load_directly=False,
    bound_list=None,
)
```

Fit the hierarchical two-orbit Doppler model with PSO.

**Bound order**

```text
iphase, semimajor, orbitalperiod, eccentricity, periapsis,
inclination, Rstar, Mstar1, Mstar2,
iphase2, semimajor2, orbitalperiod2, eccentricity2,
periapsis2, inclination2, Mass3, feature, wind_vel
```

**Returns**

```python
results, outer_phase, inner_phase, predicted_data, chi_squared
```

---

## `fit_disc_ls`

```python
fit_disc_ls(
    x_data,
    y_data,
    y_err=0,
    units="keV",
    method_="extended",
    extended_binsize=0.01,
    load_directly=False,
    bound_list=None,
)
```

Fit the same hierarchical model with bounded non-linear least squares.

**Returns**

```python
results, outer_phase, inner_phase, predicted_data, r_squared
```

---

## `fit_spiral_ps`

```python
fit_spiral_ps(
    x_data,
    y_data,
    y_err=0,
    num_iterations=3,
    maxiter=1000,
    swarmsize=100,
    units="keV",
    method_="extended",
    extended_binsize=0.01,
    load_directly=False,
    bound_list=None,
)
```

Fit an isolated logarithmic-spiral Doppler model with PSO.

**Bound order**

```text
iphase_spiral, semimajor_spiral, b, omega,
inclination_spiral, feature
```

**Returns**

```python
results, spiral_phase, predicted_data, chi_squared
```

---

## `fit_spiral_ls`

```python
fit_spiral_ls(
    x_data,
    y_data,
    y_err=0,
    units="keV",
    method_="extended",
    extended_binsize=0.01,
    load_directly=False,
    bound_list=None,
)
```

Fit the isolated spiral model with bounded non-linear least squares.

**Returns**

```python
results, spiral_phase, predicted_data, r_squared
```

---

## `fit_spiral_in_orbit_ps`

```python
fit_spiral_in_orbit_ps(
    x_data,
    y_data,
    y_err=0,
    num_iterations=3,
    maxiter=1000,
    swarmsize=100,
    units="keV",
    method_="extended",
    extended_binsize=0.01,
    load_directly=False,
    bound_list=None,
)
```

Fit a combined eccentric-orbit and logarithmic-spiral Doppler model with PSO.

**Bound order**

```text
iphase_orbit, semimajor_orbit, orbitalperiod, eccentricity,
periapsis, inclination, Rstar, Mstar1, Mstar2,
iphase_spiral, semimajor_spiral, b, omega,
inclination_spiral, feature
```

**Returns**

```python
results, orbital_phase, predicted_data, chi_squared
```

The current fitting function returns the main orbital phase. The spiral phase can be reconstructed from the fitted `iphase_spiral` and `omega` values.

---

## `fit_spiral_in_orbit_ls`

```python
fit_spiral_in_orbit_ls(
    x_data,
    y_data,
    y_err=0,
    units="keV",
    method_="extended",
    extended_binsize=0.01,
    load_directly=False,
    bound_list=None,
)
```

Fit the combined orbit--spiral model with bounded non-linear least squares.

**Returns**

```python
results, orbital_phase, predicted_data, r_squared
```

---

## `fit_nh_ps`

```python
fit_nh_ps(
    x_data,
    y_data,
    y_err=0,
    num_iterations=3,
    maxiter=200,
    swarmsize=20,
    method_="extended",
    extended_binsize=0.01,
    load_directly=False,
    bound_list=None,
)
```

Fit the smooth-wind equivalent hydrogen column with PSO.

**Bound order**

```text
iphase, semimajor, orbitalperiod, eccentricity, periapsis,
inclination, Rstar, Mdot, v_inf, beta
```

**Returns**

```python
results, orbital_phase, predicted_data, chi_squared
```

Masses are not fitted because the orbital period fixes the phase--time relation and the absorption geometry depends on the relative separation from the donor. The current fitting kernel uses a hydrogen mass fraction of 0.70 and integrates the line of sight to 1000 donor radii.

---

## `fit_nh_old`

`fit_nh_old.py` is a deprecated compatibility layer. It accepts the historical mass arguments but forwards the calculation to the corrected mass-independent `N_H` implementation. It should not be treated as a second physical model.

---

# Timing functions

## `preprocess_data`

```python
preprocess_data(t, x, sy)
```

Convert aligned arrays to NumPy arrays, remove rows containing non-finite values or non-positive uncertainties, and sort the retained data by increasing time.

**Returns**

```python
clean_time, clean_values, clean_errors
```

---

## `hr`

```python
hr(hard, soft, hard_error, soft_error)
```

Return the hardness ratio

$$
HR=\frac{H-S}{H+S}
$$

and its propagated one-sigma uncertainty. Inputs can be scalars or broadcast-compatible arrays. Invalid divisions return `NaN`.

**Returns**

```python
hardness_ratio, hardness_error
```

---

## `cr`

```python
cr(hard, soft, hard_error, soft_error)
```

Return the colour ratio

$$
CR=\frac{H}{S}
$$

and its propagated one-sigma uncertainty. Invalid divisions return `NaN`.

**Returns**

```python
colour_ratio, colour_error
```

---

## `rebin_snr`

```python
rebin_snr(
    t,
    c,
    sc,
    min_snr=5.0,
    min_bins=1,
    keep_partial=False,
)
```

Accumulate consecutive points and combine them with inverse-variance weights until

```text
abs(weighted_signal) / weighted_uncertainty >= min_snr
```

and at least `min_bins` points have been collected.

- `min_snr` uses the conventional signal/noise definition. `min_snr=5` means S/N at least 5.
- `min_bins` is the minimum number of input points per output group.
- `keep_partial=True` retains the final incomplete group even when it does not reach the target S/N.

**Returns**

```python
rebinned_time, rebinned_values, rebinned_errors
```

---

## `rebin_bins`

```python
rebin_bins(t, c, sc, nbin, keep_partial=True)
```

Group a fixed number of consecutive input points and combine them with inverse-variance weights.

`nbin` is a number of **input points**, not a duration in seconds. With `keep_partial=True`, the final shorter group is retained.

**Returns**

```python
rebinned_time, rebinned_values, rebinned_errors
```

---

## `rebin_resolution`

```python
rebin_resolution(
    t,
    c,
    sc,
    resolution,
    bin_width=None,
    start_time=None,
    keep_partial=True,
)
```

Rebin a rate, flux, or another quantity per unit time onto fixed-width temporal intervals while preserving the integrated quantity through overlap weighting.

- `resolution`: output-bin width in the same units as `t`.
- `bin_width`: scalar or array containing the width of each input bin. If omitted, widths are estimated from time centres.
- `start_time`: left edge of the first output interval.
- `keep_partial`: retain a final output interval shorter than `resolution`.

**Returns**

```python
rebinned_time, rebinned_values, rebinned_errors, coverage
```

`coverage` is the effective temporal coverage contributing to each output bin.

!!! warning "Gapped light curves"
    For gapped data, provide the real `bin_width` or exposure array. Widths inferred only from bin centres cannot reliably distinguish a gap from a broad input bin.

---

## `fold_pulse`

```python
fold_pulse(
    t,
    c,
    sc,
    period,
    snr=None,
    rebin=None,
    min_bins=1,
    keep_partial=False,
)
```

Fold a light curve on a positive period and sort it by pulse phase.

Behaviour:

- `snr` supplied: rebin the folded profile with `rebin_snr`;
- `rebin` supplied: group the folded profile with `rebin_bins`;
- neither supplied: return the unbinned, phase-sorted profile;
- both supplied: raise `ValueError`.

**Returns**

```python
pulse_phase, folded_profile, folded_errors
```

---

## `period_sliding_window`

```python
period_sliding_window(
    t,
    c,
    sc,
    window_sec,
    step_sec,
    max_period=None,
    min_period=None,
    false_alarm_threshold=0.1,
    rel_high_for_error=0.9,
    folded_pulses=False,
    snr_pulse=5.0,
    nbin_pulse=None,
    samples_per_peak=1000,
)
```

Search for uncertainty-weighted Lomb--Scargle periods in successive time windows.

### Important controls

- `window_sec`: duration of each window in the same units as `t`.
- `step_sec`: temporal displacement between consecutive windows.
- `min_period`, `max_period`: searched period interval.
- `false_alarm_threshold`: retain local peaks with FAP below this value.
- `rel_high_for_error`: relative height passed to the peak-width estimator.
- `samples_per_peak`: frequency-grid sampling density. It does not set the number of independent trials or the significance.
- `folded_pulses`: construct a folded profile for every accepted candidate.
- `snr_pulse`, `nbin_pulse`: optional folded-profile rebinning controls.

### Result columns

The result table contains:

```text
min_time
max_time
Frequency
Period
Power
Freq_Error
Period_Error
Power_Error
False_alarm
snr
```

`Freq_Error` is half the measured peak width, `Period_Error` is obtained by propagation through `P=1/f`, `Power_Error` is the robust periodogram-power scatter, and `snr` is the robust peak-power S/N.

**Returns**

```python
result, pulses
```

- `result`: `pandas.DataFrame` containing all retained peaks from all windows;
- `pulses`: dictionary of folded profiles when requested, otherwise an empty dictionary.

If there are too few valid data or no accepted periods, `result` may be `None`; check it before applying DataFrame operations.

---

# Advanced model kernels

The fitting modules also expose pure model functions. These accept explicit physical parameters and are useful for simulations, custom objective functions, numerical gradients, or external optimizers. They do not manage parameter forms or fitting bounds.

## `conic_orbit`

```python
conic_orbit(
    x_data, iphase, semimajor, orbitalperiod, eccentricity,
    periapsis, inclination, Rstar, Mstar1, Mstar2,
    wind_vel, feature, units, method_, extended_binsize,
)
```

Return the Doppler-shifted feature for a simple eccentric orbit.

## `disc_in_orbit`

```python
disc_in_orbit(
    x_data, iphase, semimajor, orbitalperiod, eccentricity,
    periapsis, inclination, Rstar, Mstar1, Mstar2,
    iphase2, semimajor2, orbitalperiod2, eccentricity2,
    periapsis2, inclination2, Mass3, feature, wind_vel,
    units, method_, extended_binsize,
)
```

Return the feature shifted by the summed line-of-sight velocities of the outer and inner orbits.

## `spiral`

```python
spiral(
    x_data, iphase_spiral, semimajor_spiral, b, omega,
    inclination_spiral, feature, units, method_, extended_binsize,
)
```

Return the Doppler-shifted feature for an isolated logarithmic spiral.

## `spiral_orbit`

```python
spiral_orbit(
    x_data, iphase_orbit, semimajor_orbit, orbitalperiod,
    eccentricity, periapsis, inclination_orbit, Rstar,
    Mstar1, Mstar2, iphase_spiral, semimajor_spiral,
    b, omega, inclination_spiral, feature, units,
    method_, extended_binsize,
)
```

Return the feature shifted by the combined binary-orbit and spiral velocities.

## `nh_orbit`

```python
nh_orbit(
    x_data, iphase, semimajor, orbitalperiod, eccentricity,
    periapsis, inclination, Rstar, Mass_loss_rate,
    wind_infinite_velocity, beta, method_, extended_binsize,
)
```

Return the equivalent hydrogen column in units of 10^22 cm^-2.

All corrected kernels support `method_="discrete"` and `method_="extended"` and validate the main dimensions, units, eccentricities, and positive physical scales.

---

## Quick output summary

| Function family | Final statistic returned | Meaning of `Std` |
|---|---|---|
| `fit_*_ps` | weighted `chi_squared` | run-to-run scatter among independent PSO best solutions |
| `fit_*_ls` | `r_squared` | local covariance-based error |
| `fit_nh_ps` | weighted `chi_squared` | run-to-run scatter among independent PSO best solutions |
| `period_sliding_window` | result table and pulse dictionary | peak-width timing errors; unrelated to fitting-table `Std` |

---

## Interpretation warnings

- A parameter at a bound is not securely measured.
- A small PSO `Std` demonstrates repeatability, not a formal confidence interval.
- A finite LS covariance error can still be misleading for a degenerate or boundary-limited fit.
- The absolute value of chi-square requires meaningful observational uncertainties.
- `R^2` is descriptive and does not replace residual inspection or a likelihood-based goodness-of-fit analysis.
- Lomb--Scargle FAP does not automatically account for red noise or every trial introduced by many overlapping windows.
- The spiral and disc models are kinematic descriptions, not hydrodynamic proofs.
- For reproducible work, store `parameter_list`, `bound_list`, the package version, and the random-seed strategy used outside the package.
