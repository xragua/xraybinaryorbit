"""Regression tests for the hierarchical disc-in-orbit model."""

from __future__ import annotations

import inspect

import numpy as np
import pytest

import xraybinaryorbit as xbo

DISC = inspect.getmodule(xbo.fit_disc_ps)
ORBIT = inspect.getmodule(xbo.fit_orbit_ps)


def _disc_parameters(mass3=0.0):
    return dict(
        iphase=0.1,
        semimajor=2.0,
        orbitalperiod=4.0,
        eccentricity=0.2,
        periapsis=30.0,
        inclination=70.0,
        Rstar=10.0,
        Mstar1=1.4,
        Mstar2=20.0,
        iphase2=0.3,
        semimajor2=0.2,
        orbitalperiod2=0.1,
        eccentricity2=0.1,
        periapsis2=50.0,
        inclination2=65.0,
        Mass3=mass3,
        feature=6.4,
        wind_vel=300.0,
    )


def test_disc_uses_correct_mass_convention_and_accepts_test_particle(monkeypatch):
    calls = []

    def fake_velocity(
        times,
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
        reference_time=None,
    ):
        calls.append(
            {
                "emitter_mass": emitter_mass,
                "companion_mass": companion_mass,
                "wind_vel": wind_vel,
                "reference_time": reference_time,
            }
        )
        value = 2.0e5 if len(calls) == 1 else -5.0e4
        return np.full(np.atleast_1d(times).size, value)

    monkeypatch.setattr(DISC, "_orbital_los_velocity", fake_velocity)

    params = _disc_parameters(mass3=0.0)
    times = np.array([100.0, 200.0])
    result = DISC.disc_in_orbit(
        times,
        *params.values(),
        "keV",
        "discrete",
        0.01,
    )

    assert calls[0]["emitter_mass"] == pytest.approx(params["Mstar1"])
    assert calls[0]["companion_mass"] == pytest.approx(params["Mstar2"])
    assert calls[0]["wind_vel"] == pytest.approx(params["wind_vel"])

    assert calls[1]["emitter_mass"] == 0.0
    assert calls[1]["companion_mass"] == pytest.approx(params["Mstar1"])
    assert calls[1]["wind_vel"] == 0.0
    assert calls[0]["reference_time"] == calls[1]["reference_time"] == 100.0

    expected = ORBIT._apply_doppler(params["feature"], 1.5e5, "keV")
    assert np.allclose(result, expected)


def test_disc_extended_evaluates_each_orbit_once_with_common_origin(monkeypatch):
    calls = []

    def fake_mean_velocity(
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
        calls.append((emitter_mass, companion_mass, reference_time))
        if len(calls) == 1:
            return np.array([1.0e5, 2.0e5])
        return np.array([-2.0e4, 3.0e4])

    monkeypatch.setattr(DISC, "_mean_orbit_velocity_in_bins", fake_mean_velocity)

    params = _disc_parameters(mass3=0.2)
    bins = np.array([[30.0, 31.0], [50.0, 49.0]])
    result = DISC.disc_in_orbit(
        bins,
        *params.values(),
        "s",
        "extended",
        0.01,
    )

    assert len(calls) == 2
    assert calls[0] == pytest.approx((1.6, 20.0, 30.0))
    assert calls[1] == pytest.approx((0.2, 1.4, 30.0))

    expected_velocity = np.array([8.0e4, 2.3e5])
    expected = ORBIT._apply_doppler(params["feature"], expected_velocity, "s")
    assert np.allclose(result, expected)


def test_disc_rejects_negative_third_mass():
    params = _disc_parameters(mass3=-1.0e-6)
    with pytest.raises(ValueError, match="Mass3"):
        DISC.disc_in_orbit(
            np.array([0.0]),
            *params.values(),
            "keV",
            "discrete",
            0.01,
        )
