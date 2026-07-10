"""Physical and regression tests for the conic-orbit fitting model."""

from __future__ import annotations

import inspect

import numpy as np
import pytest

import xraybinaryorbit as xbo

ORBIT = inspect.getmodule(xbo.fit_orbit_ps)


def test_apply_doppler_preserves_zero_velocity_and_sign():
    feature = 6.4
    velocity = np.array([-2.0e5, 0.0, 2.0e5])

    energy = ORBIT._apply_doppler(feature, velocity, "keV")
    period = ORBIT._apply_doppler(feature, velocity, "s")
    wavelength = ORBIT._apply_doppler(feature, velocity, "amstrong")

    assert energy[1] == pytest.approx(feature)
    assert period[1] == pytest.approx(feature)
    assert wavelength[1] == pytest.approx(feature)
    assert energy[0] > feature > energy[2]
    assert period[0] < feature < period[2]
    assert np.allclose(period, wavelength)


def test_orbital_velocity_matches_finite_difference_of_projected_position(monkeypatch):
    phase = 0.173
    theta_dot = 2.7e-5

    def fake_phase_and_speed(times, *args, **kwargs):
        size = np.atleast_1d(times).size
        return np.full(size, phase), np.full(size, theta_dot)

    monkeypatch.setattr(ORBIT, "_orbital_phase_and_speed", fake_phase_and_speed)

    semimajor = 2.3
    eccentricity = 0.37
    periapsis = 41.0
    inclination = 68.0
    rstar = 12.0
    emitter_mass = 1.4
    companion_mass = 19.0

    calculated = ORBIT._orbital_los_velocity(
        [100.0],
        0.0,
        semimajor,
        2.0,
        eccentricity,
        periapsis,
        inclination,
        rstar,
        emitter_mass,
        companion_mass,
        0.0,
        reference_time=0.0,
    )[0]

    a_emitter = semimajor * companion_mass / (emitter_mass + companion_mass)
    p_emitter = a_emitter * (1.0 - eccentricity**2)
    sin_i = np.sin(np.deg2rad(inclination))
    omega_peri = np.deg2rad(periapsis)

    def projected_position(theta):
        radius = p_emitter / (1.0 + eccentricity * np.cos(theta - omega_peri))
        return radius * np.cos(theta) * rstar * ORBIT.R_SUN_M * sin_i

    theta = 2.0 * np.pi * phase
    step = 1.0e-6
    dz_dtheta = (
        projected_position(theta + step) - projected_position(theta - step)
    ) / (2.0 * step)
    expected = dz_dtheta * theta_dot

    assert calculated == pytest.approx(expected, rel=2.0e-7, abs=1.0e-4)


def test_zero_inclination_gives_zero_orbital_and_wind_velocity(monkeypatch):
    def fake_phase_and_speed(times, *args, **kwargs):
        size = np.atleast_1d(times).size
        return np.full(size, 0.31), np.full(size, 1.0e-4)

    monkeypatch.setattr(ORBIT, "_orbital_phase_and_speed", fake_phase_and_speed)

    velocity = ORBIT._orbital_los_velocity(
        [0.0, 1.0],
        0.0,
        2.0,
        3.0,
        0.2,
        20.0,
        0.0,
        10.0,
        1.4,
        20.0,
        900.0,
    )

    assert np.allclose(velocity, 0.0)


def test_extended_orbit_averages_velocity_before_doppler_and_uses_global_origin(monkeypatch):
    calls = {"phase": [], "velocity": []}

    def fake_phase(times, orbit_parameters, reference_time=None):
        times = np.asarray(times, dtype=float)
        calls["phase"].append((times.copy(), reference_time))
        return 0.02 * (times - reference_time)

    def fake_velocity(times, *args, reference_time=None):
        times = np.asarray(times, dtype=float)
        calls["velocity"].append((times.copy(), reference_time))
        return 2.0e7 + 3.0e6 * (times - reference_time) ** 2

    monkeypatch.setattr(ORBIT, "_phase_at_times", fake_phase)
    monkeypatch.setattr(ORBIT, "_orbital_los_velocity", fake_velocity)

    bins = np.array([[10.0, 12.0], [20.0, 19.0]])
    feature = 6.4
    result = ORBIT.conic_orbit(
        bins,
        0.1,
        2.0,
        3.0,
        0.2,
        40.0,
        70.0,
        10.0,
        1.4,
        20.0,
        0.0,
        feature,
        "keV",
        "extended",
        0.01,
    )

    lo = np.minimum(bins[:, 0], bins[:, 1])
    hi = np.maximum(bins[:, 0], bins[:, 1])
    phase_width = 0.02 * (hi - lo)
    sample_times, sizes = ORBIT._build_sample_blocks(lo, hi, phase_width, 0.01)
    sample_velocity = 2.0e7 + 3.0e6 * (sample_times - 10.0) ** 2
    expected_velocity = ORBIT._mean_blocks(sample_velocity, sizes)
    expected = ORBIT._apply_doppler(feature, expected_velocity, "keV")

    assert np.allclose(result, expected)
    assert len(calls["phase"]) == 1
    assert len(calls["velocity"]) == 1
    assert calls["phase"][0][1] == 10.0
    assert calls["velocity"][0][1] == 10.0

    transformed_samples = ORBIT._apply_doppler(feature, sample_velocity, "keV")
    wrong_order = ORBIT._mean_blocks(transformed_samples, sizes)
    assert not np.allclose(result, wrong_order, rtol=0.0, atol=1.0e-12)


def test_orbit_rejects_invalid_units_and_shapes():
    with pytest.raises(ValueError, match="units"):
        ORBIT._canonical_units("invalid")

    with pytest.raises(ValueError, match="one-dimensional"):
        ORBIT.conic_orbit(
            np.array([[0.0, 1.0]]),
            0.0,
            2.0,
            3.0,
            0.0,
            0.0,
            60.0,
            10.0,
            1.4,
            20.0,
            0.0,
            6.4,
            "keV",
            "discrete",
            0.01,
        )
