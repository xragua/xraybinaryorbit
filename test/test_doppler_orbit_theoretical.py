import numpy as np
import pytest

from xraybinaryorbit import doppler_orbit_theoretical


@pytest.mark.parametrize("units", ["keV", "s"])
def test_doppler_orbit_theoretical_units(units):
    t = np.array([0.0, 100.0, 200.0, 300.0, 400.0])
    t_out, phase, shifted = doppler_orbit_theoretical(t, units=units, load_directly=True)
    assert len(t_out) == len(phase) == len(shifted) == len(t)
    assert isinstance(shifted, np.ndarray)
    assert np.all(np.isfinite(shifted))


def test_doppler_orbit_theoretical_invalid_units():
    with pytest.raises((KeyError, ValueError)):
        doppler_orbit_theoretical(np.array([0.0, 100.0]), units="invalid", load_directly=True)


def test_doppler_orbit_theoretical_plot():
    doppler_orbit_theoretical(
        np.array([0.0, 100.0, 200.0]),
        units="keV", show_plot=True, load_directly=True,
    )
