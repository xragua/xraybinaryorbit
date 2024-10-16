import pytest
from xraybinaryorbit import *

# Test for the default case with units in keV
def test_doppler_orbit_theoretical_keV():
    # Example time array
    t = np.array([0, 100, 200, 300, 400])

    # Run the function with keV units and load parameters directly from the file
    t_out, x_out, equation_out = doppler_orbit_theoretical(t, units="keV", load_directly=True)

    # Assert outputs
    assert len(t_out) == len(t), "The output time array should match the input length."
    assert len(x_out) == len(t), "The output orbital phase array should match the input length."
    assert len(equation_out) == len(t), "The Doppler variation array should match the input length."
    assert isinstance(equation_out, np.ndarray), "Doppler variation should be a NumPy array."


# Test for "s" unit conversion
def test_doppler_orbit_theoretical_seconds():
    # Example time array
    t = np.array([0, 100, 200, 300])

    # Run the function with seconds as the unit and load parameters directly
    t_out, x_out, equation_out = doppler_orbit_theoretical(t, units="s", load_directly=True)

    # Assert outputs
    assert len(t_out) == len(t), "The output time array should match the input length."
    assert len(x_out) == len(t), "The output orbital phase array should match the input length."
    assert len(equation_out) == len(t), "The Doppler variation array should match the input length."
    assert isinstance(equation_out, np.ndarray), "Doppler variation should be a NumPy array."


# Test for invalid unit handling
def test_doppler_orbit_theoretical_invalid_units():
    # Example time array
    t = np.array([0, 100, 200])

    # Expect a KeyError when an invalid unit is passed
    with pytest.raises(KeyError):
        doppler_orbit_theoretical(t, units="invalid_unit", load_directly=True)


# Test for plot generation
def test_doppler_orbit_theoretical_plot():
    # Example time array
    t = np.array([0, 100, 200, 300])

    # Run the function with plot generation enabled (should not raise any exceptions)
    try:
        doppler_orbit_theoretical(t, units="keV", show_plot=True, load_directly=True)
    except Exception as e:
        pytest.fail(f"Plot generation failed with error: {e}")
import pytest
import numpy as np
import sys

sys.path.append('/Users/graci/Desktop/git/xraybinaryorbit/xraybinaryorbit')
from xraybinaryorbit import *

# Test for the default case with units in keV
def test_doppler_orbit_theoretical_keV():
    # Example time array
    t = np.array([0, 100, 200, 300, 400])

    # Run the function with keV units and load parameters directly from the file
    t_out, x_out, equation_out = doppler_orbit_theoretical(t, units="keV", load_directly=True)

    # Assert outputs
    assert len(t_out) == len(t), "The output time array should match the input length."
    assert len(x_out) == len(t), "The output orbital phase array should match the input length."
    assert len(equation_out) == len(t), "The Doppler variation array should match the input length."
    assert isinstance(equation_out, np.ndarray), "Doppler variation should be a NumPy array."


# Test for "s" unit conversion
def test_doppler_orbit_theoretical_seconds():
    # Example time array
    t = np.array([0, 100, 200, 300])

    # Run the function with seconds as the unit and load parameters directly
    t_out, x_out, equation_out = doppler_orbit_theoretical(t, units="s", load_directly=True)

    # Assert outputs
    assert len(t_out) == len(t), "The output time array should match the input length."
    assert len(x_out) == len(t), "The output orbital phase array should match the input length."
    assert len(equation_out) == len(t), "The Doppler variation array should match the input length."
    assert isinstance(equation_out, np.ndarray), "Doppler variation should be a NumPy array."


# Test for invalid unit handling
def test_doppler_orbit_theoretical_invalid_units():
    # Example time array
    t = np.array([0, 100, 200])

    # Expect a KeyError when an invalid unit is passed
    with pytest.raises(KeyError):
        doppler_orbit_theoretical(t, units="invalid_unit", load_directly=True)


# Test for plot generation
def test_doppler_orbit_theoretical_plot():
    # Example time array
    t = np.array([0, 100, 200, 300])

    # Run the function with plot generation enabled (should not raise any exceptions)
    try:
        doppler_orbit_theoretical(t, units="keV", show_plot=True, load_directly=True)
    except Exception as e:
        pytest.fail(f"Plot generation failed with error: {e}")
