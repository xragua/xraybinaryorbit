"""Regression tests for the combined spiral-plus-orbit model."""

from __future__ import annotations

import inspect

import numpy as np
import pytest

import xraybinaryorbit as xbo

COMBINED = inspect.getmodule(xbo.fit_spiral_in_orbit_ps)
ORBIT = inspect.getmodule(xbo.fit_orbit_ps)


def _parameters():
    return dict(
        iphase_orbit=0.1,
        semimajor_orbit=2.0,
        orbitalperiod=3.0,
        eccentricity=0.2,
        periapsis=30.0,
        inclination_orbit=70.0,
        Rstar=12.0,
        Mstar1=1.4,
        Mstar2=20.0,
        iphase_spiral=0.2,
        semimajor_spiral=0.3,
        b=0.04,
        omega=2.0e-4,
        inclination_spiral=65.0,
        feature=6.4,
    )


def test_combined_discrete_sums_velocities_and_converts_spiral_scale(monkeypatch):
    captured = {}

    def fake_orbit(times, *args, reference_time=None):
        captured["orbit_reference"] = reference_time
        return np.full(np.atleast_1d(times).size, 1.2e5)

    def fake_spiral(
        times,
        iphase,
        semimajor_solar,
        b,
        omega,
        inclination,
        reference_time,
    ):
        captured["spiral_reference"] = reference_time
        captured["spiral_scale"] = semimajor_solar
        return np.full(np.atleast_1d(times).size, -2.0e4)

    monkeypatch.setattr(COMBINED, "_orbital_los_velocity", fake_orbit)
    monkeypatch.setattr(COMBINED, "_spiral_los_velocity", fake_spiral)

    params = _parameters()
    times = np.array([50.0, 60.0])
    result = COMBINED.spiral_orbit(
        times,
        *params.values(),
        "keV",
        "discrete",
        0.01,
    )

    assert captured["orbit_reference"] == 50.0
    assert captured["spiral_reference"] == 50.0
    assert captured["spiral_scale"] == pytest.approx(
        params["semimajor_spiral"] * params["Rstar"]
    )
    expected = ORBIT._apply_doppler(params["feature"], 1.0e5, "keV")
    assert np.allclose(result, expected)


def test_combined_extended_uses_one_global_origin_and_averages_total_velocity(monkeypatch):
    references = []

    def fake_phase(times, orbit_parameters, reference_time=None):
        references.append(reference_time)
        times = np.asarray(times, dtype=float)
        return 0.02 * (times - reference_time)

    def fake_spiral_phase(times, iphase, omega, reference_time):
        references.append(reference_time)
        times = np.asarray(times, dtype=float)
        return 0.03 * (times - reference_time)

    def fake_orbit(times, *args, reference_time=None):
        references.append(reference_time)
        times = np.asarray(times, dtype=float)
        return 1.0e7 + 2.0e6 * (times - reference_time)

    def fake_spiral_velocity(
        times,
        iphase,
        semimajor,
        b,
        omega,
        inclination,
        reference_time,
    ):
        references.append(reference_time)
        times = np.asarray(times, dtype=float)
        return -4.0e6 + 1.0e6 * (times - reference_time) ** 2

    monkeypatch.setattr(COMBINED, "_phase_at_times", fake_phase)
    monkeypatch.setattr(COMBINED, "_spiral_phase", fake_spiral_phase)
    monkeypatch.setattr(COMBINED, "_orbital_los_velocity", fake_orbit)
    monkeypatch.setattr(COMBINED, "_spiral_los_velocity", fake_spiral_velocity)

    params = _parameters()
    bins = np.array([[10.0, 12.0], [17.0, 16.0]])
    result = COMBINED.spiral_orbit(
        bins,
        *params.values(),
        "s",
        "extended",
        0.01,
    )

    lo = np.minimum(bins[:, 0], bins[:, 1])
    hi = np.maximum(bins[:, 0], bins[:, 1])
    widths = np.maximum(0.02 * (hi - lo), 0.03 * (hi - lo))
    sample_times, sizes = COMBINED._build_sample_blocks(lo, hi, widths, 0.01)
    total_velocity = (
        1.0e7
        + 2.0e6 * (sample_times - 10.0)
        - 4.0e6
        + 1.0e6 * (sample_times - 10.0) ** 2
    )
    expected = ORBIT._apply_doppler(
        params["feature"], COMBINED._mean_blocks(total_velocity, sizes), "s"
    )

    assert np.allclose(result, expected)
    assert references and set(references) == {10.0}


def test_combined_requires_positive_stellar_radius():
    params = _parameters()
    params["Rstar"] = 0.0
    with pytest.raises(ValueError, match="Rstar"):
        COMBINED.spiral_orbit(
            np.array([0.0]),
            *params.values(),
            "keV",
            "discrete",
            0.01,
        )
