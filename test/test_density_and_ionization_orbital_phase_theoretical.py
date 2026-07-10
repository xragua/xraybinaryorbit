import numpy as np
import pytest

from xraybinaryorbit import density_and_ionization_orbital_phase_theoretical


@pytest.mark.parametrize("resolution,size", [(0.01, 10), (0.001, 10), (0.01, 20)])
def test_density_and_ionization_profile(resolution, size):
    output = density_and_ionization_orbital_phase_theoretical(
        resolution=resolution, size=size, show_plot=False, load_directly=True,
    )
    assert len(output) in (3, 4)
    z, density, log_xi = output[:3]
    assert len(z) == len(density) == len(log_xi)
    assert isinstance(density, np.ndarray)
    assert isinstance(log_xi, np.ndarray)
    assert len(z) > 0
    assert np.nanmin(density) >= 0
    if len(output) == 4:
        nh_remaining = output[3]
        assert len(nh_remaining) == len(z)
        assert np.nanmin(nh_remaining) >= 0


def test_density_and_ionization_profile_plot():
    density_and_ionization_orbital_phase_theoretical(
        resolution=0.01, size=10, show_plot=True, load_directly=True,
    )
