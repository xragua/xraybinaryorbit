import pytest
from xraybinaryorbit import *

# Test for the default case with plot generation disabled
def test_absorption_column_through_orbit_theoretical_default():

    # Run the function
    time_out, phase_out, nh_out = absorption_column_through_orbit_theoretical(resolution=0.01, show_plot=False, load_directly=True)

    # Assert outputs
    assert len(time_out) == len(phase_out) == len(nh_out), "The output arrays should have the same length."
    assert isinstance(nh_out, np.ndarray), "NH1 should be a NumPy array."
    assert len(nh_out) > 0, "The NH1 array should not be empty."
    assert nh_out.min() >= 0, "The NH1 values should be non-negative."


# Test for higher resolution
def test_absorption_column_through_orbit_theoretical_high_resolution():

    # Run the function with higher resolution
    time_out, phase_out, nh_out = absorption_column_through_orbit_theoretical(resolution=0.001, show_plot=False, load_directly=True)

    # Assert outputs
    assert len(time_out) == len(phase_out) == len(nh_out), "The output arrays should have the same length."
    assert isinstance(nh_out, np.ndarray), "NH1 should be a NumPy array."
    assert len(nh_out) >= 1000, "The NH1 array should have at least 1000 points with higher resolution."
    assert nh_out.min() >= 0, "The NH1 values should be non-negative."


# Test for plot generation (it should not raise any exceptions)
def test_absorption_column_through_orbit_theoretical_plot():

    # Run the function with plot generation enabled (it should not raise any exceptions)
    try:
        absorption_column_through_orbit_theoretical(resolution=0.01, show_plot=True, load_directly=True)
    except Exception as e:
        pytest.fail(f"Plot generation failed with error: {e}")
