"""Fit Doppler modulation from an expanding logarithmic spiral."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from pyswarm import pso
from scipy.optimize import curve_fit

from ..helpers.data_helpers import _define_x_y_sy, _manage_bounds
from ..helpers.math_helpers import _chi_squared_weighted
from .fit_conic_orbit import _apply_doppler, _canonical_units

R_SUN_M = 696_340_000.0


def _spiral_phase(times, iphase_spiral, omega, reference_time):
    """Unwrapped spiral phase in cycles; omega is in cycles per second."""
    return (
        np.asarray(times, dtype=float) - float(reference_time)
    ) * float(omega) + float(iphase_spiral)


def _spiral_los_velocity(
    times,
    iphase_spiral,
    semimajor_spiral,
    b,
    omega,
    inclination_spiral,
    reference_time,
):
    """Full line-of-sight velocity of a logarithmic spiral.

    The spiral is R = a exp(b theta), theta = 2 pi phase. The scale
    ``semimajor_spiral`` is in solar radii and ``omega`` is in cycles s^-1.
    Both dR/dt and the azimuthal term are included.
    """
    values = (
        semimajor_spiral,
        b,
        omega,
        inclination_spiral,
        reference_time,
    )
    if not np.all(np.isfinite(np.asarray(values, dtype=float))):
        raise ValueError("Spiral parameters must be finite.")
    if semimajor_spiral <= 0:
        raise ValueError("semimajor_spiral must be positive.")

    phase = _spiral_phase(times, iphase_spiral, omega, reference_time)
    theta = 2.0 * np.pi * phase
    theta_dot = 2.0 * np.pi * float(omega)

    exponent = float(b) * theta
    if np.any(exponent > 700):
        raise OverflowError("The spiral radius overflows for these parameters.")
    radius_m = float(semimajor_spiral) * np.exp(exponent) * R_SUN_M
    sin_i = np.sin(np.deg2rad(inclination_spiral))

    # z = R cos(theta) sin(i), with R_dot = b R theta_dot.
    return (
        radius_m
        * theta_dot
        * (float(b) * np.cos(theta) - np.sin(theta))
        * sin_i
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


def spiral(
    x_data,
    iphase_spiral,
    semimajor_spiral,
    b,
    omega,
    inclination_spiral,
    feature,
    units,
    method_,
    extended_binsize,
):
    """Evaluate a Doppler-shifted feature for logarithmic-spiral motion.

    The extended workflow preserves the historical fitting operation: average
    velocity in each bin first and apply the Doppler transformation afterward.
    The corrected velocity includes both radial expansion and azimuthal motion.
    """
    units = _canonical_units(units)
    method_ = str(method_).lower()
    if method_ not in {"discrete", "extended"}:
        raise ValueError("method_ must be 'discrete' or 'extended'.")

    x = np.asarray(x_data, dtype=float)
    if x.size == 0 or not np.all(np.isfinite(x)):
        raise ValueError("x_data must contain finite times.")

    if method_ == "discrete":
        if x.ndim != 1:
            raise ValueError("Discrete x_data must be one-dimensional.")
        reference_time = float(np.min(x))
        velocity = _spiral_los_velocity(
            x,
            iphase_spiral,
            semimajor_spiral,
            b,
            omega,
            inclination_spiral,
            reference_time,
        )
        return _apply_doppler(feature, velocity, units)

    if x.ndim != 2 or x.shape[1] != 2:
        raise ValueError("Extended x_data must have shape (n_bins, 2).")

    lo = np.minimum(x[:, 0], x[:, 1])
    hi = np.maximum(x[:, 0], x[:, 1])
    reference_time = float(np.min(lo))
    edge_phase = _spiral_phase(
        np.column_stack((lo, hi)).reshape(-1),
        iphase_spiral,
        omega,
        reference_time,
    ).reshape(-1, 2)
    phase_width = np.abs(edge_phase[:, 1] - edge_phase[:, 0])

    sample_times, sizes = _build_sample_blocks(
        lo, hi, phase_width, extended_binsize
    )
    sample_velocity = _spiral_los_velocity(
        sample_times,
        iphase_spiral,
        semimajor_spiral,
        b,
        omega,
        inclination_spiral,
        reference_time,
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
        "iphase_spiral",
        "semimajor_spiral",
        "b",
        "omega",
        "inclination_spiral",
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
    for name in ("semimajor_spiral", "feature"):
        if lower[names.index(name)] <= 0:
            raise ValueError(f"The lower bound of {name} must be strictly positive.")
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
    return _spiral_phase(times, params[0], params[3], reference_time)


def _safe_score(model, y, sigma, params):
    try:
        prediction = np.asarray(model(params), dtype=float)
        if prediction.shape != y.shape or not np.all(np.isfinite(prediction)):
            return np.inf
        return float(_chi_squared_weighted(y, sigma, prediction))
    except (ValueError, FloatingPointError, OverflowError, ZeroDivisionError):
        return np.inf


def fit_spiral_ps(
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
    """Fit the spiral model using particle-swarm optimisation."""
    names = _parameter_names()
    x, y, sigma, method_ = _prepare_fit_data(x_data, y_data, y_err, method_)
    lower, upper = _manage_bounds(
        names,
        "bounds_spiral",
        load_directly=load_directly,
        bound_list=bound_list,
    )
    lower, upper = _validate_bounds(names, lower, upper)
    if num_iterations < 1:
        raise ValueError("num_iterations must be at least 1.")

    def model_from_params(params):
        return spiral(x, *params, units, method_, extended_binsize)

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
                "PSO did not find a finite spiral model; verify the bounds."
            )
        runs.append(params)
        chi_values.append(score)

    runs = np.asarray(runs, dtype=float)
    best = runs[int(np.argmin(chi_values))]
    scatter = np.std(runs, axis=0)
    predicted = np.asarray(model_from_params(best), dtype=float)
    chi_squared = float(_chi_squared_weighted(y, sigma, predicted))

    return (
        _results_frame(names, best, scatter),
        _fitted_phase(x, best, method_),
        predicted,
        chi_squared,
    )


def fit_spiral_ls(
    x_data,
    y_data,
    y_err=0,
    units="keV",
    method_="extended",
    extended_binsize=0.01,
    load_directly=False,
    bound_list=None,
):
    """Fit the spiral model using bounded non-linear least squares."""
    names = _parameter_names()
    x, y, sigma, method_ = _prepare_fit_data(x_data, y_data, y_err, method_)
    lower, upper = _manage_bounds(
        names,
        "bounds_spiral",
        load_directly=load_directly,
        bound_list=bound_list,
    )
    lower, upper = _validate_bounds(names, lower, upper)

    def model(model_x, *params):
        return spiral(model_x, *params, units, method_, extended_binsize)

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
