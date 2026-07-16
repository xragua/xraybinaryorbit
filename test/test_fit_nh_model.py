"""Physical and regression tests for the spherical-wind NH model."""

from __future__ import annotations

import inspect

import numpy as np
import pytest

import xraybinaryorbit as xbo

NH = inspect.getmodule(xbo.fit_nh_ps)

def test_eclipsed_line_of_sight_is_marked_as_nan():
    value = NH._nh_at_phase(
        phase=0.0,
        semimajor=2.0,
        eccentricity=0.0,
        periapsis=0.0,
        inclination=90.0,
        Rstar=20.0,
        mass_loss_rate=1.0e-6,
        wind_infinite_velocity=1000.0,
        beta=0.8,
    )

    assert np.isnan(value)


def test_nh_scales_linearly_with_mass_loss_and_inverse_wind_speed():
    common = dict(
        phase=0.5,
        semimajor=2.5,
        eccentricity=0.0,
        periapsis=0.0,
        inclination=60.0,
        Rstar=20.0,
        beta=0.8,
    )

    base = NH._nh_at_phase(
        **common,
        mass_loss_rate=1.0e-6,
        wind_infinite_velocity=1000.0,
    )
    double_mdot = NH._nh_at_phase(
        **common,
        mass_loss_rate=2.0e-6,
        wind_infinite_velocity=1000.0,
    )
    double_velocity = NH._nh_at_phase(
        **common,
        mass_loss_rate=1.0e-6,
        wind_infinite_velocity=2000.0,
    )

    assert base > 0.0
    assert double_mdot == pytest.approx(2.0 * base, rel=2.0e-6)
    assert double_velocity == pytest.approx(0.5 * base, rel=2.0e-6)


def test_extended_nh_uses_two_global_phase_calls_and_averages_samples(monkeypatch):
    calls = []

    def fake_phase(
        times,
        iphase,
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        Rstar,
        reference_time=None,
    ):
        times = np.asarray(times, dtype=float)
        calls.append((times.copy(), reference_time))
        return 0.02 * (times - reference_time)

    def fake_nh(phases, parameters):
        phases = np.asarray(phases, dtype=float)
        return 1.0 + 10.0 * phases**2

    monkeypatch.setattr(NH, "_phase_from_times", fake_phase)
    monkeypatch.setattr(NH, "_nh_for_phases", fake_nh)

    bins = np.array([[10.0, 12.0], [17.0, 16.0]])
    result = NH.nh_orbit(
        bins,
        0.0,
        2.0,
        3.0,
        0.1,
        0.0,
        60.0,
        20.0,
        1.0e-6,
        1000.0,
        0.8,
        "extended",
        0.01,
    )

    lo = np.minimum(bins[:, 0], bins[:, 1])
    hi = np.maximum(bins[:, 0], bins[:, 1])
    width = 0.02 * (hi - lo)
    sample_times, sizes = NH._build_sample_blocks(lo, hi, width, 0.01)
    phases = 0.02 * (sample_times - 10.0)
    expected = NH._mean_blocks(1.0 + 10.0 * phases**2, sizes)

    assert np.allclose(result, expected)
    assert len(calls) == 2
    assert calls[0][1] == calls[1][1] == 10.0


def test_nh_parameter_validation():
    with pytest.raises(ValueError, match="eccentricity"):
        NH.nh_orbit(
            np.array([0.0]),
            0.0,
            2.0,
            3.0,
            1.0,
            0.0,
            60.0,
            20.0,
            1.0e-6,
            1000.0,
            0.8,
            "discrete",
            0.01,
        )
