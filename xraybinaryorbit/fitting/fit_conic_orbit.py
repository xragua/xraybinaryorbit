"""Fit Doppler modulation produced by an eccentric binary orbit."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from pyswarm import pso
from scipy.optimize import curve_fit

from ..helpers.data_helpers import _define_x_y_sy, _manage_bounds
from ..helpers.math_helpers import _chi_squared_weighted, _orbital_time_to_phase

C_LIGHT = 299_792_458.0  # m s^-1
R_SUN_M = 696_340_000.0


def _canonical_units(units: str) -> str:
    """Return a canonical unit name while accepting the historical typo."""
    units = {"amstrong": "angstrom"}.get(str(units), str(units))
    if units not in {"keV", "s", "angstrom"}:
        raise ValueError("units must be 'keV', 's', or 'angstrom'.")
    return units


def _apply_doppler(feature, velocity, units):
    """Apply the package's non-relativistic Doppler convention."""
    units = _canonical_units(units)
    feature = float(feature)
    velocity = np.asarray(velocity, dtype=float)

    if not np.isfinite(feature) or feature <= 0:
        raise ValueError("feature must be a positive finite value.")

    factor = 1.0 + velocity / C_LIGHT
    if np.any(factor <= 0) or not np.all(np.isfinite(factor)):
        raise ValueError("The fitted velocity produces an invalid Doppler factor.")

    return feature / factor if units == "keV" else feature * factor


def _validate_orbit(
    semimajor,
    orbitalperiod,
    eccentricity,
    Rstar,
    emitter_mass,
    companion_mass,
):
    if semimajor <= 0 or orbitalperiod <= 0 or Rstar <= 0:
        raise ValueError("semimajor, orbitalperiod, and Rstar must be positive.")
    # A zero-mass emitter is useful for the test-particle limit in fit_disc.
    if emitter_mass < 0 or companion_mass <= 0 or emitter_mass + companion_mass <= 0:
        raise ValueError(
            "The emitter mass must be non-negative and the companion mass positive."
        )
    if not 0 <= eccentricity < 1:
        raise ValueError("eccentricity must satisfy 0 <= eccentricity < 1.")


def _anchored_times(times, reference_time):
    """Prepend a common time origin for helpers that internally subtract min(t)."""
    times = np.atleast_1d(np.asarray(times, dtype=float))
    if times.size == 0 or not np.all(np.isfinite(times)):
        raise ValueError("times must contain finite values.")

    if reference_time is None:
        return times, False

    reference_time = float(reference_time)
    if not np.isfinite(reference_time):
        raise ValueError("reference_time must be finite.")
    if reference_time > float(np.min(times)):
        raise ValueError("reference_time cannot be later than the evaluated times.")

    return np.concatenate(([reference_time], times)), True


def _orbital_phase_and_speed(
    times,
    iphase,
    semimajor,
    orbitalperiod,
    eccentricity,
    periapsis,
    Rstar,
    emitter_mass,
    companion_mass,
    reference_time=None,
):
    evaluation_times, prepended = _anchored_times(times, reference_time)
    phase, _, angular_speed = _orbital_time_to_phase(
        evaluation_times,
        iphase,
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        Rstar,
        emitter_mass,
        companion_mass,
        precision=0.01,
    )
    phase = np.asarray(phase, dtype=float)
    angular_speed = np.asarray(angular_speed, dtype=float)
    if prepended:
        phase = phase[1:]
        angular_speed = angular_speed[1:]
    return phase, angular_speed


def _orbital_los_velocity(
    times,
    iphase,
    semimajor,
    orbitalperiod,
    eccentricity,
    periapsis,
    inclination,
    Rstar,
    Mstar1,
    Mstar2,
    wind_vel,
    reference_time=None,
):
    """Line-of-sight velocity of Mstar1 around the system barycentre.

    ``semimajor`` is the relative semimajor axis in units of ``Rstar``.
    ``Mstar1`` carries the feature and ``Mstar2`` is its companion. The full
    derivative of the eccentric radius is included. ``Mstar1=0`` is accepted
    as the test-particle limit required by the inner-disc model.
    """
    _validate_orbit(
        semimajor, orbitalperiod, eccentricity, Rstar, Mstar1, Mstar2
    )

    phase, angular_speed = _orbital_phase_and_speed(
        times,
        iphase,
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        Rstar,
        Mstar1,
        Mstar2,
        reference_time=reference_time,
    )

    theta = 2.0 * np.pi * phase
    true_anomaly = theta - np.deg2rad(periapsis)
    sin_i = np.sin(np.deg2rad(inclination))

    # Barycentric orbit of the feature-carrying object Mstar1.
    a_emitter = semimajor * Mstar2 / (Mstar1 + Mstar2)
    p_emitter = a_emitter * (1.0 - eccentricity**2)
    denominator = 1.0 + eccentricity * np.cos(true_anomaly)
    if np.any(denominator <= 0):
        raise ValueError("Invalid eccentric-orbit denominator.")

    radius = p_emitter / denominator
    dr_dnu = (
        p_emitter
        * eccentricity
        * np.sin(true_anomaly)
        / denominator**2
    )
    radial_speed = dr_dnu * angular_speed  # Rstar s^-1

    # z = R cos(theta) sin(i); differentiate both R and theta.
    orbital_velocity = (
        radial_speed * np.cos(theta)
        - radius * angular_speed * np.sin(theta)
    ) * (Rstar * R_SUN_M) * sin_i

    # Historical optional wind/outflow term retained for API compatibility.
    wind_velocity = float(wind_vel) * 1000.0 * np.cos(theta) * sin_i
    return orbital_velocity + wind_velocity


def _phase_at_times(times, orbit_parameters, reference_time=None):
    phase, _ = _orbital_phase_and_speed(
        times,
        orbit_parameters[0],
        orbit_parameters[1],
        orbit_parameters[2],
        orbit_parameters[3],
        orbit_parameters[4],
        orbit_parameters[6],
        orbit_parameters[7],
        orbit_parameters[8],
        reference_time=reference_time,
    )
    return phase


def _samples_for_bin(phase_width, extended_binsize):
    if extended_binsize <= 0:
        raise ValueError("extended_binsize must be positive.")
    ratio = max(float(phase_width) / float(extended_binsize), 1.0)
    return int(np.clip(np.ceil(8.0 * ratio) + 1, 17, 129))


def _build_sample_blocks(lo, hi, phase_width, extended_binsize):
    blocks = []
    sizes = []
    for start, stop, width in zip(lo, hi, phase_width):
        midpoint = 0.5 * (start + stop)
        if stop > start and width >= extended_binsize:
            block = np.linspace(
                start,
                stop,
                _samples_for_bin(width, extended_binsize),
            )
        else:
            block = np.array([midpoint], dtype=float)
        blocks.append(block)
        sizes.append(block.size)
    return np.concatenate(blocks), np.asarray(sizes, dtype=int)


def _mean_blocks(values, sizes):
    split = np.cumsum(sizes)[:-1]
    return np.asarray(
        [np.mean(block) for block in np.split(np.asarray(values, dtype=float), split)],
        dtype=float,
    )


def conic_orbit(
    x_data,
    iphase,
    semimajor,
    orbitalperiod,
    eccentricity,
    periapsis,
    inclination,
    Rstar,
    Mstar1,
    Mstar2,
    wind_vel,
    feature,
    units,
    method_,
    extended_binsize,
):
    """Evaluate the Doppler-shifted feature for points or time bins.

    The extended workflow follows the historical implementation: first average
    the physical line-of-sight velocity inside each bin, then apply the Doppler
    transformation once to that mean velocity. All phase calculations share
    the same global time origin.
    """
    units = _canonical_units(units)
    method_ = str(method_).lower()
    if method_ not in {"discrete", "extended"}:
        raise ValueError("method_ must be 'discrete' or 'extended'.")

    x = np.asarray(x_data, dtype=float)
    if x.size == 0 or not np.all(np.isfinite(x)):
        raise ValueError("x_data must contain finite times.")

    orbit_parameters = (
        iphase,
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        inclination,
        Rstar,
        Mstar1,
        Mstar2,
        wind_vel,
    )

    if method_ == "discrete":
        if x.ndim != 1:
            raise ValueError("Discrete x_data must be one-dimensional.")
        reference_time = float(np.min(x))
        velocity = _orbital_los_velocity(
            x, *orbit_parameters, reference_time=reference_time
        )
        return _apply_doppler(feature, velocity, units)

    if x.ndim != 2 or x.shape[1] != 2:
        raise ValueError("Extended x_data must have shape (n_bins, 2).")

    lo = np.minimum(x[:, 0], x[:, 1])
    hi = np.maximum(x[:, 0], x[:, 1])
    reference_time = float(np.min(lo))

    edge_times = np.column_stack((lo, hi)).reshape(-1)
    edge_phase = _phase_at_times(
        edge_times,
        orbit_parameters,
        reference_time=reference_time,
    ).reshape(-1, 2)
    phase_width = np.abs(edge_phase[:, 1] - edge_phase[:, 0])

    sample_times, sizes = _build_sample_blocks(
        lo, hi, phase_width, extended_binsize
    )
    sample_velocity = _orbital_los_velocity(
        sample_times,
        *orbit_parameters,
        reference_time=reference_time,
    )
    mean_velocity = _mean_blocks(sample_velocity, sizes)
    return _apply_doppler(feature, mean_velocity, units)


def _prepare_fit_data(x_data, y_data, y_err, method_):
    y = np.asarray(y_data, dtype=float).reshape(-1)
    if y.size == 0 or not np.all(np.isfinite(y)):
        raise ValueError("y_data must be a non-empty finite array.")

    x, sigma = _define_x_y_sy(x_data, y, y_err)
    x = np.asarray(x, dtype=float)
    sigma = np.asarray(sigma, dtype=float).reshape(-1)
    if sigma.size != y.size or np.any(~np.isfinite(sigma)) or np.any(sigma <= 0):
        raise ValueError("y_err must provide one positive finite uncertainty per point.")

    method_ = str(method_).lower()
    if method_ not in {"discrete", "extended"}:
        raise ValueError("method_ must be 'discrete' or 'extended'.")

    if x.ndim == 1:
        if x.size != y.size:
            raise ValueError("x_data and y_data have incompatible lengths.")
        if method_ == "extended":
            warnings.warn(
                "One-dimensional x_data cannot represent bins; using discrete mode.",
                RuntimeWarning,
                stacklevel=2,
            )
        method_ = "discrete"
    elif x.ndim == 2 and x.shape == (y.size, 2):
        if method_ == "discrete":
            x = np.mean(x, axis=1)
    else:
        raise ValueError("x_data must contain points or an (n, 2) bin array.")

    return x, y, sigma, method_


def _parameter_names():
    return [
        "iphase",
        "semimajor",
        "orbitalperiod",
        "eccentricity",
        "periapsis",
        "inclination",
        "Rstar",
        "Mstar1",
        "Mstar2",
        "wind_vel",
        "feature",
    ]


def _validate_bounds(names, lower, upper):
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    if lower.shape != upper.shape or lower.size != len(names):
        raise ValueError("The number of bounds does not match the model parameters.")
    if not np.all(np.isfinite(lower)) or not np.all(np.isfinite(upper)):
        raise ValueError("All parameter bounds must be finite.")
    if np.any(lower >= upper):
        bad = [name for name, lo, hi in zip(names, lower, upper) if lo >= hi]
        raise ValueError(f"Lower bounds must be smaller than upper bounds: {bad}.")

    for name in ("semimajor", "orbitalperiod", "Rstar", "Mstar1", "Mstar2", "feature"):
        if lower[names.index(name)] <= 0:
            raise ValueError(f"The lower bound of {name} must be strictly positive.")

    e_index = names.index("eccentricity")
    if lower[e_index] < 0 or upper[e_index] >= 1:
        raise ValueError("eccentricity bounds must satisfy 0 <= e < 1.")
    return lower, upper


def _results_frame(names, values, errors):
    frame = pd.DataFrame(
        [values, errors], index=["Value", "Std"], columns=names
    )
    frame.columns.name = "Name of the parameter"
    return frame


def _fitted_phase(x_data, params, method_):
    if method_ == "extended":
        times = np.mean(x_data, axis=1)
        reference_time = float(np.min(x_data))
    else:
        times = np.asarray(x_data, dtype=float)
        reference_time = float(np.min(times))

    orbit_parameters = tuple(params[:10])
    return _phase_at_times(
        times, orbit_parameters, reference_time=reference_time
    )


def _safe_score(model, y, sigma, params):
    try:
        prediction = np.asarray(model(params), dtype=float)
        if prediction.shape != y.shape or not np.all(np.isfinite(prediction)):
            return np.inf
        return float(_chi_squared_weighted(y, sigma, prediction))
    except (ValueError, FloatingPointError, OverflowError, ZeroDivisionError):
        return np.inf


def fit_orbit_ps(
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
):
    """Fit the orbital model with particle-swarm optimisation."""
    names = _parameter_names()
    x, y, sigma, method_ = _prepare_fit_data(x_data, y_data, y_err, method_)
    lower, upper = _manage_bounds(
        names,
        "bounds_orbit",
        load_directly=load_directly,
        bound_list=bound_list,
    )
    lower, upper = _validate_bounds(names, lower, upper)
    if num_iterations < 1:
        raise ValueError("num_iterations must be at least 1.")

    def model_from_params(params):
        return conic_orbit(
            x, *params, units, method_, extended_binsize
        )

    def objective(params):
        return _safe_score(model_from_params, y, sigma, params)

    fitted = []
    chi_values = []
    for _ in range(num_iterations):
        params, _ = pso(
            objective,
            lb=lower,
            ub=upper,
            maxiter=maxiter,
            swarmsize=swarmsize,
        )
        params = np.asarray(params, dtype=float)
        score = objective(params)
        if not np.isfinite(score):
            raise RuntimeError(
                "PSO did not find a finite orbital model; verify the bounds."
            )
        fitted.append(params)
        chi_values.append(score)

    fitted = np.asarray(fitted, dtype=float)
    best_params = fitted[int(np.argmin(chi_values))]
    scatter = np.std(fitted, axis=0)
    predicted = np.asarray(model_from_params(best_params), dtype=float)
    chi_squared = float(_chi_squared_weighted(y, sigma, predicted))

    return (
        _results_frame(names, best_params, scatter),
        _fitted_phase(x, best_params, method_),
        predicted,
        chi_squared,
    )


def fit_orbit_ls(
    x_data,
    y_data,
    y_err=0,
    units="keV",
    method_="extended",
    extended_binsize=0.01,
    load_directly=False,
    bound_list=None,
):
    """Fit the orbital model with bounded non-linear least squares."""
    names = _parameter_names()
    x, y, sigma, method_ = _prepare_fit_data(x_data, y_data, y_err, method_)
    lower, upper = _manage_bounds(
        names,
        "bounds_orbit",
        load_directly=load_directly,
        bound_list=bound_list,
    )
    lower, upper = _validate_bounds(names, lower, upper)

    def model(model_x, *params):
        return conic_orbit(
            model_x, *params, units, method_, extended_binsize
        )

    has_errors = np.any(np.asarray(y_err, dtype=float) != 0)
    try:
        fit_params, covariance = curve_fit(
            model,
            x,
            y,
            p0=0.5 * (lower + upper),
            bounds=(lower, upper),
            sigma=sigma if has_errors else None,
            absolute_sigma=bool(has_errors),
            maxfev=200_000,
        )
    except (RuntimeError, ValueError) as exc:
        raise RuntimeError(
            "Curve fitting did not converge; verify the data and parameter bounds."
        ) from exc

    errors = np.sqrt(np.clip(np.diag(covariance), 0.0, np.inf))
    predicted = np.asarray(model(x, *fit_params), dtype=float)
    residuals = y - predicted
    tss = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = np.nan if tss == 0 else 1.0 - float(np.sum(residuals**2)) / tss

    return (
        _results_frame(names, fit_params, errors),
        _fitted_phase(x, fit_params, method_),
        predicted,
        r_squared,
    )
