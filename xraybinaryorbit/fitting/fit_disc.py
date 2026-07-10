"""Fit Doppler modulation from an inner orbit embedded in a binary orbit."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from pyswarm import pso
from scipy.optimize import curve_fit

from ..helpers.data_helpers import _define_x_y_sy, _manage_bounds
from ..helpers.math_helpers import _chi_squared_weighted
from .fit_conic_orbit import (
    _apply_doppler,
    _canonical_units,
    _orbital_los_velocity,
    _phase_at_times,
)


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


def _mean_orbit_velocity_in_bins(
    bins,
    iphase,
    semimajor,
    orbitalperiod,
    eccentricity,
    periapsis,
    inclination,
    Rstar,
    emitter_mass,
    companion_mass,
    wind_vel,
    extended_binsize,
    reference_time,
):
    bins = np.asarray(bins, dtype=float)
    lo = np.minimum(bins[:, 0], bins[:, 1])
    hi = np.maximum(bins[:, 0], bins[:, 1])

    orbit_parameters = (
        iphase,
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        inclination,
        Rstar,
        emitter_mass,
        companion_mass,
        wind_vel,
    )
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
    velocity = _orbital_los_velocity(
        sample_times,
        *orbit_parameters,
        reference_time=reference_time,
    )
    return _mean_blocks(velocity, sizes)


def disc_in_orbit(
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
    iphase2,
    semimajor2,
    orbitalperiod2,
    eccentricity2,
    periapsis2,
    inclination2,
    Mass3,
    feature,
    wind_vel,
    units,
    method_,
    extended_binsize,
):
    """Evaluate Doppler modulation from two hierarchical orbits.

    ``Mstar1`` is the compact object at the centre of the inner orbit,
    ``Mstar2`` is the donor/outer companion, and ``Mass3`` is the body or disc
    element carrying the feature. ``Mass3=0`` is accepted as the test-particle
    limit. ``semimajor`` and ``semimajor2`` are relative semimajor axes in
    units of ``Rstar``.

    Extended bins preserve the historical workflow: average each physical
    velocity in the bin, add the two mean velocities, then apply one Doppler
    transformation. All phase calls use the same global time origin.
    """
    units = _canonical_units(units)
    method_ = str(method_).lower()
    if method_ not in {"discrete", "extended"}:
        raise ValueError("method_ must be 'discrete' or 'extended'.")
    if Mstar1 <= 0 or Mstar2 <= 0 or Mass3 < 0:
        raise ValueError(
            "Mstar1 and Mstar2 must be positive and Mass3 non-negative."
        )

    x = np.asarray(x_data, dtype=float)
    if x.size == 0 or not np.all(np.isfinite(x)):
        raise ValueError("x_data must contain finite times.")

    outer_emitter_mass = Mstar1 + Mass3

    if method_ == "discrete":
        if x.ndim != 1:
            raise ValueError("Discrete x_data must be one-dimensional.")
        reference_time = float(np.min(x))
        v_outer = _orbital_los_velocity(
            x,
            iphase,
            semimajor,
            orbitalperiod,
            eccentricity,
            periapsis,
            inclination,
            Rstar,
            outer_emitter_mass,
            Mstar2,
            wind_vel,
            reference_time=reference_time,
        )
        v_inner = _orbital_los_velocity(
            x,
            iphase2,
            semimajor2,
            orbitalperiod2,
            eccentricity2,
            periapsis2,
            inclination2,
            Rstar,
            Mass3,
            Mstar1,
            0.0,
            reference_time=reference_time,
        )
        return _apply_doppler(feature, v_outer + v_inner, units)

    if x.ndim != 2 or x.shape[1] != 2:
        raise ValueError("Extended x_data must have shape (n_bins, 2).")

    reference_time = float(np.min(x))
    v_outer = _mean_orbit_velocity_in_bins(
        x,
        iphase,
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        inclination,
        Rstar,
        outer_emitter_mass,
        Mstar2,
        wind_vel,
        extended_binsize,
        reference_time,
    )
    v_inner = _mean_orbit_velocity_in_bins(
        x,
        iphase2,
        semimajor2,
        orbitalperiod2,
        eccentricity2,
        periapsis2,
        inclination2,
        Rstar,
        Mass3,
        Mstar1,
        0.0,
        extended_binsize,
        reference_time,
    )
    return _apply_doppler(feature, v_outer + v_inner, units)


def _prepare_fit_data(x_data, y_data, y_err, method_):
    y = np.asarray(y_data, dtype=float).reshape(-1)
    if y.size == 0 or not np.all(np.isfinite(y)):
        raise ValueError("y_data must be a non-empty finite array.")

    x, sigma = _define_x_y_sy(x_data, y, y_err)
    x = np.asarray(x, dtype=float)
    sigma = np.asarray(sigma, dtype=float).reshape(-1)
    if sigma.size != y.size or np.any(~np.isfinite(sigma)) or np.any(sigma <= 0):
        raise ValueError("y_err must contain positive finite uncertainties.")

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
        "iphase2",
        "semimajor2",
        "orbitalperiod2",
        "eccentricity2",
        "periapsis2",
        "inclination2",
        "Mass3",
        "feature",
        "wind_vel",
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

    for name in (
        "semimajor",
        "orbitalperiod",
        "Rstar",
        "Mstar1",
        "Mstar2",
        "semimajor2",
        "orbitalperiod2",
        "feature",
    ):
        if lower[names.index(name)] <= 0:
            raise ValueError(f"The lower bound of {name} must be strictly positive.")

    if lower[names.index("Mass3")] < 0:
        raise ValueError("The lower bound of Mass3 must be non-negative.")

    for name in ("eccentricity", "eccentricity2"):
        index = names.index(name)
        if lower[index] < 0 or upper[index] >= 1:
            raise ValueError(f"{name} bounds must satisfy 0 <= e < 1.")
    return lower, upper


def _results_frame(names, values, errors):
    frame = pd.DataFrame(
        [values, errors], index=["Value", "Std"], columns=names
    )
    frame.columns.name = "Name of the parameter"
    return frame


def _fitted_phases(x_data, params, method_):
    if method_ == "extended":
        times = np.mean(x_data, axis=1)
        reference_time = float(np.min(x_data))
    else:
        times = np.asarray(x_data, dtype=float)
        reference_time = float(np.min(times))

    outer_parameters = (
        params[0],
        params[1],
        params[2],
        params[3],
        params[4],
        params[5],
        params[6],
        params[7] + params[15],
        params[8],
        params[17],
    )
    inner_parameters = (
        params[9],
        params[10],
        params[11],
        params[12],
        params[13],
        params[14],
        params[6],
        params[15],
        params[7],
        0.0,
    )
    ph = _phase_at_times(
        times, outer_parameters, reference_time=reference_time
    )
    ph2 = _phase_at_times(
        times, inner_parameters, reference_time=reference_time
    )
    return ph, ph2


def _safe_score(model, y, sigma, params):
    try:
        prediction = np.asarray(model(params), dtype=float)
        if prediction.shape != y.shape or not np.all(np.isfinite(prediction)):
            return np.inf
        return float(_chi_squared_weighted(y, sigma, prediction))
    except (ValueError, FloatingPointError, OverflowError, ZeroDivisionError):
        return np.inf


def fit_disc_ps(
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
    """Fit the hierarchical two-orbit model using particle swarm optimisation."""
    names = _parameter_names()
    x, y, sigma, method_ = _prepare_fit_data(x_data, y_data, y_err, method_)
    lower, upper = _manage_bounds(
        names,
        "bounds_disc",
        load_directly=load_directly,
        bound_list=bound_list,
    )
    lower, upper = _validate_bounds(names, lower, upper)
    if num_iterations < 1:
        raise ValueError("num_iterations must be at least 1.")

    def model_from_params(params):
        return disc_in_orbit(x, *params, units, method_, extended_binsize)

    def objective(params):
        return _safe_score(model_from_params, y, sigma, params)

    runs = []
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
                "PSO did not find a finite disc model; verify the bounds."
            )
        runs.append(params)
        chi_values.append(score)

    runs = np.asarray(runs, dtype=float)
    best = runs[int(np.argmin(chi_values))]
    scatter = np.std(runs, axis=0)
    predicted = np.asarray(model_from_params(best), dtype=float)
    chi_squared = float(_chi_squared_weighted(y, sigma, predicted))
    ph, ph2 = _fitted_phases(x, best, method_)

    return _results_frame(names, best, scatter), ph, ph2, predicted, chi_squared


def fit_disc_ls(
    x_data,
    y_data,
    y_err=0,
    units="keV",
    method_="extended",
    extended_binsize=0.01,
    load_directly=False,
    bound_list=None,
):
    """Fit the hierarchical two-orbit model using bounded least squares."""
    names = _parameter_names()
    x, y, sigma, method_ = _prepare_fit_data(x_data, y_data, y_err, method_)
    lower, upper = _manage_bounds(
        names,
        "bounds_disc",
        load_directly=load_directly,
        bound_list=bound_list,
    )
    lower, upper = _validate_bounds(names, lower, upper)

    def model(model_x, *params):
        return disc_in_orbit(
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
    ph, ph2 = _fitted_phases(x, fit_params, method_)

    return _results_frame(names, fit_params, errors), ph, ph2, predicted, r_squared
