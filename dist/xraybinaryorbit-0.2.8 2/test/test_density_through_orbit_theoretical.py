import pytest
from xraybinaryorbit import *

# Test for the default case with plot generation disabled
def test_density_through_orbit_theoretical_default():

    # Run the function with default parameters
    time_out, phase_out, density_out = density_through_orbit_theoretical(resolution=0.01, show_plot=False, load_directly=True)

    # Assert outputs
    assert len(time_out) == len(phase_out) == len(density_out), "The output arrays should have the same length."
    assert isinstance(density_out, np.ndarray), "Density should be a NumPy array."
    assert len(density_out) > 0, "The density array should not be empty."
    assert density_out.min() >= 0, "The density values should be non-negative."


# Test for higher resolution
def test_density_through_orbit_theoretical_high_resolution():

    # Run the function with higher resolution
    time_out, phase_out, density_out = density_through_orbit_theoretical(resolution=0.001, show_plot=False, load_directly=True)

    # Assert outputs
    assert len(time_out) == len(phase_out) == len(density_out), "The output arrays should have the same length."
    assert isinstance(density_out, np.ndarray), "Density should be a NumPy array."
    assert len(density_out) > 0, "The density array should not be empty."
    assert density_out.min() >= 0, "The density values should be non-negative."


# Test for plot generation (it should not raise any exceptions)
def test_density_through_orbit_theoretical_plot():

    # Run the function with plot generation enabled (it should not raise any exceptions)
    try:
        density_through_orbit_theoretical(resolution=0.01, show_plot=True, load_directly=True)
    except Exception as e:
        pytest.fail(f"Plot generation failed with error: {e}")
