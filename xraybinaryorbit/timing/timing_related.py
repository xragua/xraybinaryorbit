"""Timing-analysis utilities for light curves and folded pulse profiles."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
from scipy.optimize import brentq, minimize_scalar
from scipy.signal import find_peaks


def _clean_timing_arrays(t, y, sy, *extra_arrays, sort=True):
    """Convert aligned timing arrays to NumPy arrays and remove invalid points.

    Parameters
    ----------
    t : array-like
        Time values.
    y : array-like
        Count rate, flux, or another timing quantity.
    sy : array-like
        One-sigma uncertainties on ``y``.
    *extra_arrays : array-like
        Additional one-dimensional arrays aligned with ``t``. A point is
        removed when any supplied array contains a non-finite value there.
    sort : bool, optional
        Sort all retained arrays by increasing time. Default is True.

    Returns
    -------
    tuple of numpy.ndarray
        Cleaned and aligned arrays containing finite values and strictly
        positive uncertainties.

    Raises
    ------
    ValueError
        If an input is not one-dimensional or the arrays have different
        lengths.
    """
    arrays = [
        np.asarray(t, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(sy, dtype=float),
    ]
    arrays.extend(np.asarray(arr, dtype=float) for arr in extra_arrays)

    if any(arr.ndim != 1 for arr in arrays):
        raise ValueError("All timing arrays must be one-dimensional.")

    lengths = {len(arr) for arr in arrays}
    if len(lengths) != 1:
        raise ValueError("All timing arrays must have the same length.")

    valid_mask = (
        np.isfinite(arrays[0])
        & np.isfinite(arrays[1])
        & np.isfinite(arrays[2])
        & (arrays[2] > 0)
    )

    for arr in arrays[3:]:
        valid_mask &= np.isfinite(arr)

    cleaned = [arr[valid_mask] for arr in arrays]

    if sort and len(cleaned[0]) > 1:
        order = np.argsort(cleaned[0], kind="stable")
        cleaned = [arr[order] for arr in cleaned]

    return tuple(cleaned)


def preprocess_data(t, x, sy):
    """Clean timing arrays.

    This backward-compatible wrapper now applies the same validation used by
    all rebinning and period-search routines.
    """
    return _clean_timing_arrays(t, x, sy)


def hr(x, y, ex, ey):
    """Calculate the hardness ratio and its propagated uncertainty.

    The hardness ratio is

    ``HR = (x - y) / (x + y)``.

    Invalid entries are returned as NaN without changing the broadcast shape
    of the inputs.

    Parameters
    ----------
    x, y : float or array-like
        Hard- and soft-band count rates or fluxes.
    ex, ey : float or array-like
        One-sigma uncertainties on ``x`` and ``y``.

    Returns
    -------
    hr_value, hr_error : float or numpy.ndarray
        Hardness ratio and propagated one-sigma uncertainty.
    """
    x, y, ex, ey = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(ex, dtype=float),
        np.asarray(ey, dtype=float),
    )

    denominator = x + y
    valid = (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(ex)
        & np.isfinite(ey)
        & (ex >= 0)
        & (ey >= 0)
        & (denominator != 0)
    )

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        hr_value = (x - y) / denominator
        df_dx = 2.0 * y / denominator**2
        df_dy = -2.0 * x / denominator**2
        hr_error = np.sqrt((df_dx * ex) ** 2 + (df_dy * ey) ** 2)

    hr_value = np.where(valid, hr_value, np.nan)
    hr_error = np.where(valid, hr_error, np.nan)

    return hr_value, hr_error


def cr(x, y, ex, ey):
    """Calculate the color ratio and its propagated uncertainty.

    The color ratio is ``CR = x / y``. Invalid entries are returned as NaN
    without changing the broadcast shape of the inputs.

    Parameters
    ----------
    x, y : float or array-like
        Hard- and soft-band count rates or fluxes.
    ex, ey : float or array-like
        One-sigma uncertainties on ``x`` and ``y``.

    Returns
    -------
    cr_value, cr_error : float or numpy.ndarray
        Color ratio and propagated one-sigma uncertainty.
    """
    x, y, ex, ey = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(ex, dtype=float),
        np.asarray(ey, dtype=float),
    )

    valid = (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(ex)
        & np.isfinite(ey)
        & (ex >= 0)
        & (ey >= 0)
        & (y != 0)
    )

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        cr_value = x / y
        df_dx = 1.0 / y
        df_dy = -x / y**2
        cr_error = np.sqrt((df_dx * ex) ** 2 + (df_dy * ey) ** 2)

    cr_value = np.where(valid, cr_value, np.nan)
    cr_error = np.where(valid, cr_error, np.nan)

    return cr_value, cr_error


def _combine_inverse_variance(t, y, sy):
    """Return an inverse-variance-weighted timing bin."""
    weights = 1.0 / sy**2
    weight_sum = np.sum(weights)

    if not np.isfinite(weight_sum) or weight_sum <= 0:
        return None

    t_weighted = np.sum(t * weights) / weight_sum
    y_weighted = np.sum(y * weights) / weight_sum
    sy_weighted = np.sqrt(1.0 / weight_sum)

    return t_weighted, y_weighted, sy_weighted


def rebin_snr(t, c, sc, min_snr=5.0, min_bins=1, keep_partial=False):
    """Rebin a light curve until a minimum signal-to-noise ratio is reached.

    Consecutive points are accumulated and combined with inverse-variance
    weights. A bin is emitted when

    ``abs(weighted_signal) / weighted_uncertainty >= min_snr``

    and at least ``min_bins`` points have been accumulated.

    Parameters
    ----------
    t : array-like
        Time values.
    c : array-like
        Count rate, flux, or another timing quantity.
    sc : array-like
        One-sigma uncertainties on ``c``.
    min_snr : float, optional
        Minimum conventional signal-to-noise ratio. Default is 5.
    min_bins : int, optional
        Minimum number of input points per output bin. Default is 1.
    keep_partial : bool, optional
        Retain the final group even if it does not reach ``min_snr``.
        Default is False.

    Returns
    -------
    t_new, c_new, sc_new : numpy.ndarray
        Rebinned times, weighted values, and propagated uncertainties.
    """
    t, c, sc = _clean_timing_arrays(t, c, sc)

    min_snr = float(min_snr)
    if not np.isfinite(min_snr) or min_snr <= 0:
        raise ValueError("min_snr must be a positive finite number.")

    if min_snr < 1:
        warnings.warn(
            "min_snr now uses signal/noise. Under the old noise/signal "
            "convention, a threshold of 0.2 corresponds to min_snr=5.",
            UserWarning,
            stacklevel=2,
        )

    if not isinstance(min_bins, (int, np.integer)) or min_bins < 1:
        raise ValueError("min_bins must be a positive integer.")

    t_new, c_new, sc_new = [], [], []
    t_bin, c_bin, sc_bin = [], [], []

    for ti, ci, sci in zip(t, c, sc):
        t_bin.append(ti)
        c_bin.append(ci)
        sc_bin.append(sci)

        combined = _combine_inverse_variance(
            np.asarray(t_bin), np.asarray(c_bin), np.asarray(sc_bin)
        )
        if combined is None:
            continue

        t_weighted, c_weighted, sc_weighted = combined
        snr_now = np.abs(c_weighted) / sc_weighted

        if snr_now >= min_snr and len(c_bin) >= min_bins:
            t_new.append(t_weighted)
            c_new.append(c_weighted)
            sc_new.append(sc_weighted)
            t_bin, c_bin, sc_bin = [], [], []

    if keep_partial and c_bin:
        combined = _combine_inverse_variance(
            np.asarray(t_bin), np.asarray(c_bin), np.asarray(sc_bin)
        )
        if combined is not None:
            t_weighted, c_weighted, sc_weighted = combined
            t_new.append(t_weighted)
            c_new.append(c_weighted)
            sc_new.append(sc_weighted)

    return np.asarray(t_new), np.asarray(c_new), np.asarray(sc_new)


def rebin_bins(t, c, sc, nbin, keep_partial=True):
    """Rebin a light curve by grouping a fixed number of consecutive points.

    Parameters
    ----------
    t : array-like
        Time values.
    c : array-like
        Count rate, flux, or another timing quantity.
    sc : array-like
        One-sigma uncertainties on ``c``.
    nbin : int
        Number of consecutive input points per output bin.
    keep_partial : bool, optional
        Retain the final group if it contains fewer than ``nbin`` points.
        Default is True.

    Returns
    -------
    t_new, c_new, sc_new : numpy.ndarray
        Rebinned times, weighted values, and propagated uncertainties.
    """
    t, c, sc = _clean_timing_arrays(t, c, sc)

    if not isinstance(nbin, (int, np.integer)) or nbin < 1:
        raise ValueError("nbin must be a positive integer.")

    t_new, c_new, sc_new = [], [], []

    for start in range(0, len(c), nbin):
        stop = min(start + nbin, len(c))
        t_bin = t[start:stop]
        c_bin = c[start:stop]
        sc_bin = sc[start:stop]

        if len(c_bin) < nbin and not keep_partial:
            continue

        combined = _combine_inverse_variance(t_bin, c_bin, sc_bin)
        if combined is None:
            continue

        t_weighted, c_weighted, sc_weighted = combined
        t_new.append(t_weighted)
        c_new.append(c_weighted)
        sc_new.append(sc_weighted)

    return np.asarray(t_new), np.asarray(c_new), np.asarray(sc_new)


def _estimate_time_bin_widths(t):
    """Estimate temporal bin widths from ordered bin-center times.

    Notes
    -----
    Large gaps cannot be distinguished reliably from broad input bins when
    only bin centers are supplied. For gapped light curves, pass the actual
    ``bin_width`` or exposure array to :func:`rebin_resolution`.
    """
    t = np.asarray(t, dtype=float)

    if len(t) == 0:
        return np.array([], dtype=float)
    if len(t) == 1:
        return np.ones_like(t, dtype=float)

    edges = np.empty(len(t) + 1, dtype=float)
    edges[1:-1] = 0.5 * (t[:-1] + t[1:])
    edges[0] = t[0] - 0.5 * (t[1] - t[0])
    edges[-1] = t[-1] + 0.5 * (t[-1] - t[-2])

    return np.diff(edges)


def rebin_resolution(
    t,
    c,
    sc,
    resolution,
    bin_width=None,
    start_time=None,
    keep_partial=True,
):
    """Rebin a light curve onto fixed-width temporal intervals.

    This routine is intended for count rates, fluxes, or other quantities per
    unit time. It preserves the integrated quantity by accounting for the
    temporal overlap between input and output bins.

    Parameters
    ----------
    t : array-like
        Input time-bin centers.
    c : array-like
        Count-rate, flux, or another time-density quantity.
    sc : array-like
        One-sigma uncertainties on ``c``.
    resolution : float
        Nominal output-bin width, in the same units as ``t``. For example,
        use ``resolution=100`` for 100-second bins when ``t`` is in seconds.
    bin_width : float or array-like, optional
        Width of each input time bin. A scalar is applied to all points. If
        omitted, widths are estimated from the time centers.
    start_time : float, optional
        Left edge of the first output bin. By default, the first input-bin
        edge is used.
    keep_partial : bool, optional
        Retain the final output interval when it is shorter than
        ``resolution``. Default is True.

    Returns
    -------
    t_new, c_new, sc_new, coverage_new : numpy.ndarray
        Exposure-weighted output times, rebinned values, propagated
        uncertainties, and effective temporal coverage in each output bin.

    Notes
    -----
    Empty intervals are omitted. If an output bin is only partly observed,
    its value is averaged over the observed coverage rather than treating the
    missing interval as zero exposure.
    """
    resolution = float(resolution)
    if not np.isfinite(resolution) or resolution <= 0:
        raise ValueError("resolution must be a positive finite number.")

    t = np.asarray(t, dtype=float)
    c = np.asarray(c, dtype=float)
    sc = np.asarray(sc, dtype=float)

    if bin_width is None:
        t, c, sc = _clean_timing_arrays(t, c, sc)
        bin_width = _estimate_time_bin_widths(t)
    else:
        bin_width = np.asarray(bin_width, dtype=float)
        if bin_width.ndim == 0:
            bin_width = np.full(len(t), float(bin_width))
        t, c, sc, bin_width = _clean_timing_arrays(t, c, sc, bin_width)

    if len(t) == 0:
        empty = np.array([], dtype=float)
        return empty, empty, empty, empty

    valid_width = np.isfinite(bin_width) & (bin_width > 0)
    t = t[valid_width]
    c = c[valid_width]
    sc = sc[valid_width]
    bin_width = bin_width[valid_width]

    if len(t) == 0:
        empty = np.array([], dtype=float)
        return empty, empty, empty, empty

    order = np.argsort(t, kind="stable")
    t = t[order]
    c = c[order]
    sc = sc[order]
    bin_width = bin_width[order]

    input_left = t - 0.5 * bin_width
    input_right = t + 0.5 * bin_width

    output_start = np.min(input_left) if start_time is None else float(start_time)
    if not np.isfinite(output_start):
        raise ValueError("start_time must be finite.")

    output_end = np.max(input_right)
    if output_start >= output_end:
        raise ValueError(
            "start_time must be earlier than the end of the light curve."
        )

    duration = output_end - output_start
    n_full_bins = int(np.floor(duration / resolution))
    edges = output_start + np.arange(n_full_bins + 1, dtype=float) * resolution

    tolerance = np.finfo(float).eps * max(
        abs(output_start), abs(output_end), resolution, 1.0
    )
    remainder = output_end - edges[-1]

    if keep_partial and remainder > tolerance:
        edges = np.append(edges, output_end)

    if len(edges) < 2:
        empty = np.array([], dtype=float)
        return empty, empty, empty, empty

    t_new, c_new, sc_new, coverage_new = [], [], [], []

    for left, right in zip(edges[:-1], edges[1:]):
        overlap_left_all = np.maximum(input_left, left)
        overlap_right_all = np.minimum(input_right, right)
        overlap_all = overlap_right_all - overlap_left_all
        mask = overlap_all > 0

        if not np.any(mask):
            continue

        overlap_left = overlap_left_all[mask]
        overlap_right = overlap_right_all[mask]
        dt = overlap_right - overlap_left
        total_coverage = np.sum(dt)

        if not np.isfinite(total_coverage) or total_coverage <= 0:
            continue

        integrated_c = np.sum(c[mask] * dt)
        integrated_sc = np.sqrt(np.sum((sc[mask] * dt) ** 2))

        c_out = integrated_c / total_coverage
        sc_out = integrated_sc / total_coverage

        overlap_centers = 0.5 * (overlap_left + overlap_right)
        t_out = np.sum(overlap_centers * dt) / total_coverage

        t_new.append(t_out)
        c_new.append(c_out)
        sc_new.append(sc_out)
        coverage_new.append(total_coverage)

    return (
        np.asarray(t_new),
        np.asarray(c_new),
        np.asarray(sc_new),
        np.asarray(coverage_new),
    )


def fold_pulse(
    t,
    c,
    sc,
    period,
    snr=None,
    rebin=None,
    min_bins=1,
    keep_partial=False,
):
    """Fold a light curve on a period and optionally rebin the pulse profile.

    Parameters
    ----------
    t, c, sc : array-like
        Times, light-curve values, and one-sigma uncertainties.
    period : float
        Folding period in the same units as ``t``.
    snr : float, optional
        Minimum conventional signal-to-noise ratio for adaptive phase bins.
    rebin : int, optional
        Number of consecutive phase-sorted points per output bin.
    min_bins : int, optional
        Minimum number of points per adaptive S/N bin. Default is 1.
    keep_partial : bool, optional
        Retain the final incomplete bin. Default is False for S/N rebinning
        and is passed to fixed-count rebinning when requested.

    Returns
    -------
    phase, c_folded, sc_folded : numpy.ndarray
        Folded, phase-sorted profile, optionally rebinned.

    Raises
    ------
    ValueError
        If the period is invalid or both ``snr`` and ``rebin`` are supplied.
    """
    t, c, sc = _clean_timing_arrays(t, c, sc)

    period = float(period)
    if not np.isfinite(period) or period <= 0:
        raise ValueError("period must be a positive finite number.")

    if snr is not None and rebin is not None:
        raise ValueError("Provide either snr or rebin, not both.")

    if len(t) == 0:
        empty = np.array([], dtype=float)
        return empty, empty, empty

    phase = np.mod(t - t[0], period) / period
    order = np.argsort(phase, kind="stable")
    phase = phase[order]
    c = c[order]
    sc = sc[order]

    if snr is not None:
        return rebin_snr(
            phase,
            c,
            sc,
            min_snr=snr,
            min_bins=min_bins,
            keep_partial=keep_partial,
        )

    if rebin is not None:
        return rebin_bins(
            phase,
            c,
            sc,
            nbin=rebin,
            keep_partial=keep_partial,
        )

    return phase, c, sc


def _robust_power_noise(power):
    """Estimate the robust scatter of a periodogram power distribution."""
    power = np.asarray(power, dtype=float)
    finite = power[np.isfinite(power)]

    if len(finite) == 0:
        return np.nan, np.nan

    baseline = np.median(finite)
    mad = np.median(np.abs(finite - baseline))
    noise = 1.4826 * mad

    if not np.isfinite(noise) or noise <= 0:
        noise = np.std(finite)

    return baseline, noise




_PERIODOGRAM_COLUMNS = [
    "min_time",
    "max_time",
    "Frequency",
    "Frequency_Lower",
    "Frequency_Upper",
    "Freq_Error",
    "Freq_Error_Minus",
    "Freq_Error_Plus",
    "Period",
    "Period_Lower",
    "Period_Upper",
    "Period_Error",
    "Period_Error_Minus",
    "Period_Error_Plus",
    "Power",
    "Power_Noise",
    "Power_Error",
    "False_alarm",
    "snr",
    "Delta_Chi2",
]

_SLIDING_WINDOW_COLUMNS = [
    "window",
    "window_start",
    "window_end",
    *_PERIODOGRAM_COLUMNS,
]


def _empty_periodogram_frame(sliding_window=False):
    """Return an empty period-search table with a stable output schema."""
    columns = (
        _SLIDING_WINDOW_COLUMNS
        if sliding_window
        else _PERIODOGRAM_COLUMNS
    )
    return pd.DataFrame({column: pd.Series(dtype=float) for column in columns})


def _profile_frequency_interval(
    ls,
    freq,
    power,
    peak_index,
    chi2_reference,
    delta_chi2=1.0,
):
    """Refine a Lomb-Scargle peak and obtain its local frequency interval.

    The interval is defined from a profile in frequency. At every trial
    frequency, the sinusoidal coefficients are re-optimized internally by
    :class:`astropy.timeseries.LombScargle`. With ``normalization='standard'``
    the relation

    ``chi2(f) - chi2_min = chi2_reference * (power_max - power(f))``

    permits the confidence limits to be found from a chosen ``delta_chi2``.
    For one parameter of interest, ``delta_chi2=1`` is the usual local
    approximately 68.3 per cent interval under Gaussian assumptions.

    Parameters
    ----------
    ls : astropy.timeseries.LombScargle
        Lomb-Scargle object using ``normalization='standard'``.
    freq, power : numpy.ndarray
        Ascending frequency grid and associated periodogram powers.
    peak_index : int
        Index of the sampled candidate maximum.
    chi2_reference : float
        Chi-square of the constant reference model.
    delta_chi2 : float, optional
        Profile-likelihood threshold. Default is 1.

    Returns
    -------
    best_frequency, lower_frequency, upper_frequency, best_power : float
        Refined maximum, lower and upper frequency limits, and refined power.
        A missing limit is returned as NaN when the requested interval reaches
        the boundary of the searched frequency range.
    """
    if (
        peak_index <= 0
        or peak_index >= len(freq) - 1
        or not np.isfinite(chi2_reference)
        or chi2_reference <= 0
    ):
        return (
            float(freq[peak_index]),
            np.nan,
            np.nan,
            float(power[peak_index]),
        )

    def scalar_power(frequency):
        return np.asarray(ls.power(frequency), dtype=float).item()

    lower_bound = float(freq[peak_index - 1])
    upper_bound = float(freq[peak_index + 1])
    frequency_step = min(
        float(freq[peak_index] - freq[peak_index - 1]),
        float(freq[peak_index + 1] - freq[peak_index]),
    )

    refinement = minimize_scalar(
        lambda frequency: -scalar_power(frequency),
        bounds=(lower_bound, upper_bound),
        method="bounded",
        options={
            "xatol": max(
                np.finfo(float).eps,
                frequency_step * 1e-6,
            )
        },
    )

    if refinement.success and np.isfinite(refinement.x):
        best_frequency = float(refinement.x)
        best_power = scalar_power(best_frequency)
    else:
        best_frequency = float(freq[peak_index])
        best_power = float(power[peak_index])

    target_power = best_power - delta_chi2 / chi2_reference

    def root_function(frequency):
        return scalar_power(frequency) - target_power

    def find_crossing(indices):
        previous_frequency = best_frequency
        previous_value = best_power - target_power

        for index in indices:
            current_frequency = float(freq[index])
            current_value = float(power[index] - target_power)

            if current_value <= 0:
                left = min(current_frequency, previous_frequency)
                right = max(current_frequency, previous_frequency)

                try:
                    return float(brentq(root_function, left, right))
                except ValueError:
                    # Linear fallback if tiny numerical differences between
                    # scalar and precomputed powers prevent root bracketing.
                    denominator = previous_value - current_value
                    if denominator == 0:
                        return np.nan

                    fraction = -current_value / denominator
                    return float(
                        current_frequency
                        + fraction
                        * (previous_frequency - current_frequency)
                    )

            previous_frequency = current_frequency
            previous_value = current_value

        return np.nan

    left_indices = np.where(freq < best_frequency)[0][::-1]
    right_indices = np.where(freq > best_frequency)[0]

    lower_frequency = find_crossing(left_indices)
    upper_frequency = find_crossing(right_indices)

    return (
        best_frequency,
        lower_frequency,
        upper_frequency,
        best_power,
    )


def _periodogram_peaks(
    t,
    c,
    sc,
    max_period,
    min_period,
    false_alarm_threshold,
    samples_per_peak,
    delta_chi2=1.0,
):
    """Calculate significant Lomb-Scargle peaks for one time window.

    The frequency limits are local profile-chi-square intervals. Period limits
    are transformed exactly through ``Period = 1 / Frequency`` and are
    therefore retained as asymmetric uncertainties.
    """
    t, c, sc = _clean_timing_arrays(t, c, sc)

    if len(t) < 4 or np.ptp(t) <= 0:
        return _empty_periodogram_frame()

    min_freq = 1.0 / max_period if max_period is not None else None
    max_freq = 1.0 / min_period if min_period is not None else None

    ls = LombScargle(
        t,
        c,
        sc,
        normalization="standard",
    )

    autopower_kwargs = {"samples_per_peak": samples_per_peak}
    if min_freq is not None:
        autopower_kwargs["minimum_frequency"] = min_freq
    if max_freq is not None:
        autopower_kwargs["maximum_frequency"] = max_freq

    freq, power = ls.autopower(**autopower_kwargs)
    freq = np.asarray(freq, dtype=float)
    power = np.asarray(power, dtype=float)

    finite = np.isfinite(freq) & np.isfinite(power) & (freq > 0)
    freq = freq[finite]
    power = power[finite]

    if len(freq) < 3:
        return _empty_periodogram_frame()

    order = np.argsort(freq, kind="stable")
    freq = freq[order]
    power = power[order]

    peak_indices = find_peaks(power)[0]
    if len(peak_indices) == 0:
        return _empty_periodogram_frame()

    fap_kwargs = {
        "method": "baluev",
        "samples_per_peak": samples_per_peak,
        "minimum_frequency": float(freq[0]),
        "maximum_frequency": float(freq[-1]),
    }

    coarse_false_alarm = np.asarray(
        ls.false_alarm_probability(
            power[peak_indices],
            **fap_kwargs,
        ),
        dtype=float,
    )

    significant = (
        np.isfinite(coarse_false_alarm)
        & (coarse_false_alarm < false_alarm_threshold)
    )
    peak_indices = peak_indices[significant]

    if len(peak_indices) == 0:
        return _empty_periodogram_frame()

    # This is the chi-square of the floating constant reference model used by
    # the standard Lomb-Scargle normalization with measurement uncertainties.
    weights = 1.0 / sc**2
    weight_sum = np.sum(weights)
    if not np.isfinite(weight_sum) or weight_sum <= 0:
        return _empty_periodogram_frame()

    weighted_mean = np.sum(weights * c) / weight_sum
    chi2_reference = float(np.sum(((c - weighted_mean) / sc) ** 2))

    if not np.isfinite(chi2_reference) or chi2_reference <= 0:
        return _empty_periodogram_frame()

    baseline, power_noise = _robust_power_noise(power)
    rows = []

    for index in peak_indices:
        (
            best_frequency,
            lower_frequency,
            upper_frequency,
            best_power,
        ) = _profile_frequency_interval(
            ls=ls,
            freq=freq,
            power=power,
            peak_index=index,
            chi2_reference=chi2_reference,
            delta_chi2=delta_chi2,
        )

        refined_false_alarm = np.asarray(
            ls.false_alarm_probability(
                best_power,
                **fap_kwargs,
            ),
            dtype=float,
        ).item()

        if (
            not np.isfinite(refined_false_alarm)
            or refined_false_alarm >= false_alarm_threshold
        ):
            continue

        period = 1.0 / best_frequency

        if np.isfinite(lower_frequency) and lower_frequency > 0:
            frequency_error_minus = best_frequency - lower_frequency
            period_upper = 1.0 / lower_frequency
            period_error_plus = period_upper - period
        else:
            lower_frequency = np.nan
            frequency_error_minus = np.nan
            period_upper = np.nan
            period_error_plus = np.nan

        if np.isfinite(upper_frequency) and upper_frequency > 0:
            frequency_error_plus = upper_frequency - best_frequency
            period_lower = 1.0 / upper_frequency
            period_error_minus = period - period_lower
        else:
            upper_frequency = np.nan
            frequency_error_plus = np.nan
            period_lower = np.nan
            period_error_minus = np.nan

        # Backward-compatible symmetric summaries. These are deliberately NaN
        # if either side of the interval is unconstrained. Scientific output
        # should use the explicit Minus and Plus columns.
        if (
            np.isfinite(frequency_error_minus)
            and np.isfinite(frequency_error_plus)
        ):
            frequency_error = 0.5 * (
                frequency_error_minus + frequency_error_plus
            )
        else:
            frequency_error = np.nan

        if (
            np.isfinite(period_error_minus)
            and np.isfinite(period_error_plus)
        ):
            period_error = 0.5 * (
                period_error_minus + period_error_plus
            )
        else:
            period_error = np.nan

        if np.isfinite(power_noise) and power_noise > 0:
            power_snr = (best_power - baseline) / power_noise
        elif best_power > baseline:
            power_snr = np.inf
        else:
            power_snr = 0.0

        rows.append(
            {
                "min_time": float(np.min(t)),
                "max_time": float(np.max(t)),
                "Frequency": best_frequency,
                "Frequency_Lower": lower_frequency,
                "Frequency_Upper": upper_frequency,
                "Freq_Error": frequency_error,
                "Freq_Error_Minus": frequency_error_minus,
                "Freq_Error_Plus": frequency_error_plus,
                "Period": period,
                "Period_Lower": period_lower,
                "Period_Upper": period_upper,
                "Period_Error": period_error,
                "Period_Error_Minus": period_error_minus,
                "Period_Error_Plus": period_error_plus,
                "Power": best_power,
                "Power_Noise": power_noise,
                # Deprecated compatibility alias: this is periodogram
                # background scatter, not a formal error on peak power.
                "Power_Error": power_noise,
                "False_alarm": refined_false_alarm,
                "snr": power_snr,
                "Delta_Chi2": delta_chi2,
            }
        )

    if not rows:
        return _empty_periodogram_frame()

    return (
        pd.DataFrame(rows, columns=_PERIODOGRAM_COLUMNS)
        .sort_values(by="Power", ascending=False)
        .reset_index(drop=True)
    )


def period_sliding_window(
    t,
    c,
    sc,
    window_sec,
    step_sec,
    max_period=None,
    min_period=None,
    false_alarm_threshold=0.1,
    rel_high_for_error=None,
    folded_pulses=False,
    snr_pulse=5.0,
    nbin_pulse=None,
    samples_per_peak=1000,
    delta_chi2=1.0,
):
    """Search for Lomb-Scargle periods in sliding time windows.

    ``window_sec`` and ``step_sec`` are time durations in the same units as
    ``t``. Peak frequencies are refined continuously around each sampled
    maximum. Their local confidence intervals are estimated from a profile in
    frequency and transformed exactly into asymmetric period intervals.

    Parameters
    ----------
    t, c, sc : array-like
        Times, light-curve values, and one-sigma uncertainties.
    window_sec : float
        Duration of each sliding window in the same units as ``t``.
    step_sec : float
        Time displacement between consecutive windows.
    max_period, min_period : float, optional
        Largest and smallest periods searched.
    false_alarm_threshold : float, optional
        Retain peaks whose global Lomb-Scargle false-alarm probability is
        below this value. Default is 0.1.
    rel_high_for_error : float or None, optional
        Deprecated compatibility parameter. Peak widths are no longer used
        as statistical errors. Supplying a value emits a deprecation warning
        and has no effect.
    folded_pulses : bool, optional
        Return folded profiles for all retained periods. Default is False.
    snr_pulse : float or None, optional
        Minimum conventional S/N for adaptive folded-profile bins. Default
        is 5. Set to None when using ``nbin_pulse``.
    nbin_pulse : int or None, optional
        Number of phase-sorted points per folded-profile bin.
    samples_per_peak : int, optional
        Lomb-Scargle frequency-grid sampling. Default is 1000.
    delta_chi2 : float, optional
        Profile threshold used for frequency limits. ``delta_chi2=1`` gives
        the usual local approximately 68.3 per cent interval for one
        parameter of interest under Gaussian assumptions.

    Returns
    -------
    result : pandas.DataFrame
        Significant periods found in all windows. The table always has a
        stable schema, even when it is empty. Use ``Period_Error_Minus`` and
        ``Period_Error_Plus`` (or ``Period_Lower`` and ``Period_Upper``) for
        scientific reporting. ``Period_Error`` is retained only as a
        backward-compatible average of both sides when both are available.
    pulses : dict
        Folded profiles keyed by result-row number when
        ``folded_pulses=True``. Every profile includes the central period,
        asymmetric frequency and period limits, and detection diagnostics.
        The profile itself is folded at the central best-fit period.
    """
    t, c, sc = _clean_timing_arrays(t, c, sc)

    window_sec = float(window_sec)
    step_sec = float(step_sec)
    delta_chi2 = float(delta_chi2)

    if not np.isfinite(window_sec) or window_sec <= 0:
        raise ValueError("window_sec must be a positive finite number.")
    if not np.isfinite(step_sec) or step_sec <= 0:
        raise ValueError("step_sec must be a positive finite number.")
    if not np.isfinite(delta_chi2) or delta_chi2 <= 0:
        raise ValueError("delta_chi2 must be a positive finite number.")

    for name, value in (("min_period", min_period), ("max_period", max_period)):
        if value is not None and (not np.isfinite(value) or value <= 0):
            raise ValueError(f"{name} must be a positive finite number or None.")

    if min_period is not None and max_period is not None:
        if min_period >= max_period:
            raise ValueError("min_period must be smaller than max_period.")

    if not 0 < false_alarm_threshold < 1:
        raise ValueError("false_alarm_threshold must lie between 0 and 1.")

    if rel_high_for_error is not None:
        warnings.warn(
            "rel_high_for_error is deprecated and ignored. Period and "
            "frequency intervals are now controlled by delta_chi2.",
            DeprecationWarning,
            stacklevel=2,
        )

    if not isinstance(samples_per_peak, (int, np.integer)) or samples_per_peak < 1:
        raise ValueError("samples_per_peak must be a positive integer.")

    if snr_pulse is not None and nbin_pulse is not None:
        raise ValueError(
            "Set either snr_pulse or nbin_pulse for folded profiles, not both."
        )

    if len(t) < 4:
        return _empty_periodogram_frame(sliding_window=True), {}

    positive_dt = np.diff(t)
    positive_dt = positive_dt[positive_dt > 0]

    if min_period is not None and len(positive_dt) > 0:
        sampling_limit = 2.0 * np.median(positive_dt)
        if min_period < sampling_limit:
            warnings.warn(
                "The requested minimum period is shorter than twice the "
                f"median sampling interval ({sampling_limit:.4g}). Its "
                "recovery may be unreliable.",
                UserWarning,
                stacklevel=2,
            )

    representative_dt = np.median(positive_dt) if len(positive_dt) > 0 else 0.0
    effective_data_end = t[-1] + representative_dt
    total_duration = effective_data_end - t[0]

    if max_period is not None and max_period > window_sec:
        warnings.warn(
            "max_period is longer than window_sec; periods comparable to or "
            "longer than a window are poorly constrained.",
            UserWarning,
            stacklevel=2,
        )

    if total_duration <= window_sec:
        window_starts = np.array([t[0]], dtype=float)
    else:
        last_start = effective_data_end - window_sec
        tolerance = np.finfo(float).eps * max(
            abs(t[0]), abs(last_start), window_sec, step_sec, 1.0
        )
        window_starts = np.arange(
            t[0], last_start + step_sec, step_sec, dtype=float
        )
        window_starts = window_starts[window_starts <= last_start + tolerance]

    frames = []

    for window_id, window_start in enumerate(window_starts):
        window_end = window_start + window_sec
        mask = (t >= window_start) & (t < window_end)

        t_window = t[mask]
        c_window = c[mask]
        sc_window = sc[mask]

        if len(t_window) < 4 or np.ptp(t_window) <= 0:
            continue

        df_periods = _periodogram_peaks(
            t_window,
            c_window,
            sc_window,
            max_period=max_period,
            min_period=min_period,
            false_alarm_threshold=false_alarm_threshold,
            samples_per_peak=samples_per_peak,
            delta_chi2=delta_chi2,
        )

        if not df_periods.empty:
            df_periods.insert(0, "window", float(window_id))
            df_periods.insert(1, "window_start", float(window_start))
            df_periods.insert(2, "window_end", float(window_end))
            frames.append(df_periods)

    if frames:
        result = pd.concat(frames, ignore_index=True)
        result = result.reindex(columns=_SLIDING_WINDOW_COLUMNS)
    else:
        result = _empty_periodogram_frame(sliding_window=True)

    pulses = {}

    if folded_pulses and not result.empty:
        for row_number, row in result.iterrows():
            mask = (t >= row["min_time"]) & (t <= row["max_time"])
            phase, c_pulse, sc_pulse = fold_pulse(
                t[mask],
                c[mask],
                sc[mask],
                row["Period"],
                snr=snr_pulse,
                rebin=nbin_pulse,
            )

            pulses[row_number] = {
                "ph_pulse": phase,
                "c_pulse": c_pulse,
                "sc_pulse": sc_pulse,
                "frequency": row["Frequency"],
                "frequency_lower": row["Frequency_Lower"],
                "frequency_upper": row["Frequency_Upper"],
                "frequency_error_minus": row["Freq_Error_Minus"],
                "frequency_error_plus": row["Freq_Error_Plus"],
                "period": row["Period"],
                "period_lower": row["Period_Lower"],
                "period_upper": row["Period_Upper"],
                "period_error_minus": row["Period_Error_Minus"],
                "period_error_plus": row["Period_Error_Plus"],
                "power": row["Power"],
                "power_noise": row["Power_Noise"],
                "false_alarm": row["False_alarm"],
                "snr": row["snr"],
                "delta_chi2": row["Delta_Chi2"],
                "window": row["window"],
                "window_start": row["window_start"],
                "window_end": row["window_end"],
            }

    return result, pulses
