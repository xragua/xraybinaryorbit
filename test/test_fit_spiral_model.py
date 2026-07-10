"""Physical and regression tests for the logarithmic-spiral model."""

from __future__ import annotations

import inspect

import numpy as np
import pytest

import xraybinaryorbit as xbo

SPIRAL = inspect.getmodule(xbo.fit_spiral_ps)
ORBIT = inspect.getmodule(xbo.fit_orbit_ps)


def test_spiral_velocity_matches_finite_difference_of_projected_position():
    iphase = 0.13
    semimajor = 1.7
    b = 0.045
    omega = 2.5e-4
    inclination = 63.0
    reference_time = 100.0
    time = 731.0

    calculated = SPIRAL._spiral_los_velocity(
        [time],
        iphase,
        semimajor,
        b,
        omega,
        inclination,
        reference_time,
    )[0]

    sin_i = np.sin(np.deg2rad(inclination))

    def projected_position(t):
        phase = (t - reference_time) * omega + iphase
        theta = 2.0 * np.pi * phase
        radius = semimajor * np.exp(b * theta) * SPIRAL.R_SUN_M
        return radius * np.cos(theta) * sin_i

    step = 1.0e-3
    expected = (projected_position(time + step) - projected_position(time - step)) / (
        2.0 * step
    )

    assert calculated == pytest.approx(expected, rel=2.0e-7, abs=1.0e-3)


def test_spiral_zero_inclination_has_no_projected_velocity():
    velocity = SPIRAL._spiral_los_velocity(
        np.array([0.0, 10.0, 20.0]),
        0.2,
        1.0,
        0.1,
        1.0e-3,
        0.0,
        0.0,
    )
    assert np.allclose(velocity, 0.0)


def test_spiral_extended_averages_velocity_before_doppler(monkeypatch):
    captured = {}

    def fake_velocity(
        times,
        iphase_spiral,
        semimajor_spiral,
        b,
        omega,
        inclination_spiral,
        reference_time,
    ):
        times = np.asarray(times, dtype=float)
        captured["times"] = times.copy()
        captured["reference_time"] = reference_time
        return 1.0e7 + 5.0e6 * (times - reference_time) ** 2

    monkeypatch.setattr(SPIRAL, "_spiral_los_velocity", fake_velocity)

    bins = np.array([[5.0, 7.0], [11.0, 10.0]])
    feature = 1.0
    result = SPIRAL.spiral(
        bins,
        0.0,
        1.0,
        0.02,
        0.02,
        80.0,
        feature,
        "keV",
        "extended",
        0.01,
    )

    lo = np.minimum(bins[:, 0], bins[:, 1])
    hi = np.maximum(bins[:, 0], bins[:, 1])
    phase_width = 0.02 * (hi - lo)
    sample_times, sizes = SPIRAL._build_sample_blocks(lo, hi, phase_width, 0.01)
    sample_velocity = 1.0e7 + 5.0e6 * (sample_times - 5.0) ** 2
    mean_velocity = SPIRAL._mean_blocks(sample_velocity, sizes)
    expected = ORBIT._apply_doppler(feature, mean_velocity, "keV")

    assert np.allclose(result, expected)
    assert captured["reference_time"] == 5.0
    assert np.array_equal(captured["times"], sample_times)

    wrong_order = SPIRAL._mean_blocks(
        ORBIT._apply_doppler(feature, sample_velocity, "keV"), sizes
    )
    assert not np.allclose(result, wrong_order, rtol=0.0, atol=1.0e-12)


def test_spiral_rejects_nonpositive_scale_and_invalid_mode():
    with pytest.raises(ValueError, match="positive"):
        SPIRAL._spiral_los_velocity(
            [0.0], 0.0, 0.0, 0.1, 1.0e-3, 60.0, 0.0
        )

    with pytest.raises(ValueError, match="method"):
        SPIRAL.spiral(
            [0.0], 0.0, 1.0, 0.1, 1.0e-3, 60.0, 6.4, "keV", "bad", 0.01
        )
