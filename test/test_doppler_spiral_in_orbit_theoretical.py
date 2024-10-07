import pytest
from xraybinaryorbit import *

# Test for doppler_spiral_in_orbit_theoretical with keV units
def test_doppler_spiral_in_orbit_theoretical_keV():
    # Example time array
    t = np.arange(0, 1000, 1)

    # Run the function with keV units and load parameters directly from the file
    t_out, x_out, x2_out, equation_out = doppler_spiral_in_orbit_theoretical(t, units="keV", load_directly=True)

    # Assert outputs
    assert len(t_out) == len(t), "The output time array should match the input length."
    assert len(x_out) == len(t), "The output orbital phase array should match the input length."
    assert len(x2_out) == len(t), "The second orbital phase array should match the input length."
    assert len(equation_out) == len(t), "The Doppler variation array should match the input length."
    assert isinstance(equation_out, np.ndarray), "Doppler variation should be a NumPy array."


# Test for plot generation with keV units
def test_doppler_spiral_in_orbit_theoretical_plot_keV():
    # Example time array
    t = np.arange(0, 1000, 1)

    # Run the function with plot generation enabled (should not raise any exceptions)
    try:
        doppler_spiral_in_orbit_theoretical(t, units="keV", show_plot=True, load_directly=True)
    except Exception as e:
        pytest.fail(f"Plot generation failed with error: {e}")
