"""Fit orbital hydrogen column density through a spherical stellar wind."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from pyswarm import pso
from scipy.integrate import quad

from ..helpers.data_helpers import _define_x_y_sy, _manage_bounds
from ..helpers.math_helpers import _chi_squared_weighted, _orbital_time_to_phase

M_SUN_G = 1.98847e33
R_SUN_CM = 6.96340e10
M_PROTON_G = 1.67262192369e-24
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0
HYDROGEN_MASS_FRACTION = 0.70
LOS_UPPER_LIMIT_RSTAR = 1000.0


def _validate_parameters(
    semimajor,
    orbitalperiod,
    eccentricity,
    Rstar,
    mass_loss_rate,
    wind_infinite_velocity,
    beta,
):
    if semimajor <= 0 or orbitalperiod <= 0 or Rstar <= 0:
        raise ValueError("semimajor, orbitalperiod, and Rstar must be positive.")
    if mass_loss_rate <= 0 or wind_infinite_velocity <= 0:
        raise ValueError("Mdot and v_inf must be positive.")
    if beta < 0:
        raise ValueError("beta must be non-negative.")
    if not 0 <= eccentricity < 1:
        raise ValueError("eccentricity must satisfy 0 <= eccentricity < 1.")


def _phase_from_times(
    times,
    iphase,
    semimajor,
    orbitalperiod,
    eccentricity,
    periapsis,
    Rstar,
    reference_time=None,
):
    times = np.atleast_1d(np.asarray(times, dtype=float))
    if times.size == 0 or not np.all(np.isfinite(times)):
        raise ValueError("times must contain finite values.")

    prepended = reference_time is not None
    if prepended:
        reference_time = float(reference_time)
        if reference_time > float(np.min(times)):
            raise ValueError("reference_time cannot be later than the evaluated times.")
        evaluation_times = np.concatenate(([reference_time], times))
    else:
        evaluation_times = times

    # Equal dummy masses avoid imposing a barycentric meaning on a relative
    # donor--compact-object separation while remaining compatible with old helpers.
    phase, _, _ = _orbital_time_to_phase(
        evaluation_times,
        iphase,
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        Rstar,
        1.0,
        1.0,
        precision=0.01,
    )
    phase = np.asarray(phase, dtype=float)
    return phase[1:] if prepended else phase


def _relative_separation_rstar(phase, semimajor, eccentricity, periapsis):
    true_anomaly = 2.0 * np.pi * np.asarray(phase, dtype=float) - np.deg2rad(periapsis)
    denominator = 1.0 + eccentricity * np.cos(true_anomaly)
    if np.any(denominator <= 0):
        raise ValueError("Invalid eccentric-orbit denominator.")
    return semimajor * (1.0 - eccentricity**2) / denominator


def _line_of_sight_intersects_star(separation, cos_alpha):
    """Check the ray from the compact object towards the observer."""
    z_closest = max(float(separation) * float(cos_alpha), 0.0)
    minimum_radius = np.sqrt(
        separation**2
        + z_closest**2
        - 2.0 * separation * z_closest * cos_alpha
    )
    return minimum_radius <= 1.0


def _nh_at_phase(
    phase,
    semimajor,
    eccentricity,
    periapsis,
    inclination,
    Rstar,
    mass_loss_rate,
    wind_infinite_velocity,
    beta,
):
    """Equivalent hydrogen column in units of 1e22 cm^-2."""
    separation = float(
        _relative_separation_rstar(phase, semimajor, eccentricity, periapsis)
    )
    cos_alpha = float(
        np.clip(
            np.cos(2.0 * np.pi * phase) * np.sin(np.deg2rad(inclination)),
            -1.0,
            1.0,
        )
    )

    # Historical package convention: an occulted direct source returns NH=0.
    if separation <= 1.0 or _line_of_sight_intersects_star(separation, cos_alpha):
        return 0.0

    rstar_cm = Rstar * R_SUN_CM
    mdot_g_s = mass_loss_rate * M_SUN_G / SECONDS_PER_YEAR
    vinf_cm_s = wind_infinite_velocity * 1.0e5

    def mass_column_integrand(z_rstar):
        radius_rstar = np.sqrt(
            separation**2
            + z_rstar**2
            - 2.0 * separation * z_rstar * cos_alpha
        )
        if radius_rstar <= 1.0:
            return np.inf
        velocity = vinf_cm_s * (1.0 - 1.0 / radius_rstar) ** beta
        if velocity <= 0 or not np.isfinite(velocity):
            return np.inf
        radius_cm = radius_rstar * rstar_cm
        density = mdot_g_s / (4.0 * np.pi * velocity * radius_cm**2)
        return density * rstar_cm  # dz = Rstar du

    mass_column, _ = quad(
        mass_column_integrand,
        0.0,
        LOS_UPPER_LIMIT_RSTAR,
        epsabs=0.0,
        epsrel=3.0e-5,
        limit=250,
    )
    if not np.isfinite(mass_column) or mass_column < 0:
        raise FloatingPointError("Non-finite wind-column integral.")

    return HYDROGEN_MASS_FRACTION * mass_column / M_PROTON_G / 1.0e22


def _nh_for_phases(phases, parameters):
    phases = np.atleast_1d(np.asarray(phases, dtype=float))
    return np.asarray(
        [
            _nh_at_phase(
                phase,
                parameters[1],
                parameters[3],
                parameters[4],
                parameters[5],
                parameters[6],
                parameters[7],
                parameters[8],
                parameters[9],
            )
            for phase in phases
        ],
        dtype=float,
    )


def _samples_for_bin(phase_width, extended_binsize):
    if extended_binsize <= 0:
        raise ValueError("extended_binsize must be positive.")
    ratio = max(float(phase_width) / float(extended_binsize), 1.0)
    return int(np.clip(np.ceil(4.0 * ratio) + 1, 9, 65))


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


def nh_orbit(
    x_data,
    iphase,
    semimajor,
    orbitalperiod,
    eccentricity,
    periapsis,
    inclination,
    Rstar,
    Mass_loss_rate,
    wind_infinite_velocity,
    beta,
    method_,
    extended_binsize,
):
    """Calculate orbital NH in units of 1e22 cm^-2.

    ``semimajor`` is the relative donor--compact-object semimajor axis in units
    of the donor radius. No barycentric mass correction belongs in this wind
    geometry. Extended bins are sampled in one globally anchored phase call and
    the physical NH values are then averaged inside each bin.
    """
    _validate_parameters(
        semimajor,
        orbitalperiod,
        eccentricity,
        Rstar,
        Mass_loss_rate,
        wind_infinite_velocity,
        beta,
    )

    method_ = str(method_).lower()
    if method_ not in {"discrete", "extended"}:
        raise ValueError("method_ must be 'discrete' or 'extended'.")

    x = np.asarray(x_data, dtype=float)
    if x.size == 0 or not np.all(np.isfinite(x)):
        raise ValueError("x_data must contain finite times.")

    parameters = (
        iphase,
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        inclination,
        Rstar,
        Mass_loss_rate,
        wind_infinite_velocity,
        beta,
    )

    if method_ == "discrete":
        if x.ndim != 1:
            raise ValueError("Discrete x_data must be one-dimensional.")
        reference_time = float(np.min(x))
        phases = _phase_from_times(
            x,
            iphase,
            semimajor,
            orbitalperiod,
            eccentricity,
            periapsis,
            Rstar,
            reference_time=reference_time,
        )
        return _nh_for_phases(phases, parameters)

    if x.ndim != 2 or x.shape[1] != 2:
        raise ValueError("Extended x_data must have shape (n_bins, 2).")

    lo = np.minimum(x[:, 0], x[:, 1])
    hi = np.maximum(x[:, 0], x[:, 1])
    reference_time = float(np.min(lo))
    edge_times = np.column_stack((lo, hi)).reshape(-1)
    edge_phases = _phase_from_times(
        edge_times,
        iphase,
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        Rstar,
        reference_time=reference_time,
    ).reshape(-1, 2)
    phase_width = np.abs(edge_phases[:, 1] - edge_phases[:, 0])

    sample_times, sizes = _build_sample_blocks(
        lo, hi, phase_width, extended_binsize
    )
    sample_phases = _phase_from_times(
        sample_times,
        iphase,
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        Rstar,
        reference_time=reference_time,
    )
    sample_nh = _nh_for_phases(sample_phases, parameters)
    return _mean_blocks(sample_nh, sizes)


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
        "Mdot",
        "v_inf",
        "beta",
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

    for name in ("semimajor", "orbitalperiod", "Rstar", "Mdot", "v_inf"):
        if lower[names.index(name)] <= 0:
            raise ValueError(f"The lower bound of {name} must be strictly positive.")
    if lower[names.index("beta")] < 0:
        raise ValueError("The lower bound of beta must be non-negative.")
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
    return _phase_from_times(
        times,
        params[0],
        params[1],
        params[2],
        params[3],
        params[4],
        params[6],
        reference_time=reference_time,
    )


def _safe_score(model, y, sigma, params):
    try:
        prediction = np.asarray(model(params), dtype=float)
        if prediction.shape != y.shape or not np.all(np.isfinite(prediction)):
            return np.inf
        return float(_chi_squared_weighted(y, sigma, prediction))
    except (ValueError, FloatingPointError, OverflowError, ZeroDivisionError):
        return np.inf


def fit_nh_ps(
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
):
    """Fit the spherical-wind NH model using particle-swarm optimisation."""
    names = _parameter_names()
    x, y, sigma, method_ = _prepare_fit_data(x_data, y_data, y_err, method_)
    lower, upper = _manage_bounds(
        names,
        "bounds_nh",
        load_directly=load_directly,
        bound_list=bound_list,
    )
    lower, upper = _validate_bounds(names, lower, upper)
    if num_iterations < 1:
        raise ValueError("num_iterations must be at least 1.")

    def model_from_params(params):
        return nh_orbit(x, *params, method_, extended_binsize)

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
            phig=2,
        )
        params = np.asarray(params, dtype=float)
        score = objective(params)
        if not np.isfinite(score):
            raise RuntimeError(
                "PSO did not find a finite NH model; verify the bounds and eclipse geometry."
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
