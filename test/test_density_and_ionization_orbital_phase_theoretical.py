import pytest
from xraybinaryorbit import *

# Test for the default case with plot generation disabled
def test_density_and_ionization_orbital_phase_theoretical_default():

    # Run the function
    z, density, chi = density_and_ionization_orbital_phase_theoretical(resolution=0.01, size=10, show_plot=False, load_directly=True)

    # Assert outputs
    assert len(z) == len(density) == len(chi), "The output arrays (z, density, chi) should have the same length."
    assert isinstance(density, np.ndarray), "Density should be a NumPy array."
    assert isinstance(chi, np.ndarray), "Chi (ionization parameter) should be a NumPy array."
    assert len(density) > 0, "The density array should not be empty."
    assert len(chi) > 0, "The chi array should not be empty."
    assert density.min() >= 0, "The density values should be non-negative."


# Test for higher resolution
def test_density_and_ionization_orbital_phase_theoretical_high_resolution():

    # Run the function with higher resolution
    z, density, chi = density_and_ionization_orbital_phase_theoretical(resolution=0.001, size=10, show_plot=False, load_directly=True)

    # Assert outputs
    assert len(z) == len(density) == len(chi), "The output arrays (z, density, chi) should have the same length."
    assert len(density) >= 1000, "The density array should have at least 1000 points with higher resolution."
    assert len(chi) >= 1000, "The chi array should have at least 1000 points with higher resolution."
    assert density.min() >= 0, "The density values should be non-negative."
    assert np.all(np.isfinite(chi)), "All chi values should be finite."


# Test for plot generation (it should not raise any exceptions)
def test_density_and_ionization_orbital_phase_theoretical_plot():

    # Run the function with plot generation enabled (it should not raise any exceptions)
    try:
        density_and_ionization_orbital_phase_theoretical(resolution=0.01, size=10, show_plot=True, load_directly=True)
    except Exception as e:
        pytest.fail(f"Plot generation failed with error: {e}")


# Test for custom size parameter
def test_density_and_ionization_orbital_phase_theoretical_custom_size():

    # Run the function with custom size parameter
    z, density, chi = density_and_ionization_orbital_phase_theoretical(resolution=0.01, size=20, show_plot=False, load_directly=True)

    # Assert outputs
    assert len(z) == len(density) == len(chi), "The output arrays (z, density, chi) should have the same length."
    assert density.min() >= 0, "The density values should be non-negative."
    assert np.all(np.isfinite(chi)), "All chi values should be finite."
    assert max(z) > 0, "The maximum z value should be greater than 0."
