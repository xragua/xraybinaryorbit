import pytest
import numpy as np
from xraybinaryorbit import doppler_disc_theoretical

# Mock functions for _manage_parameters and _orbital_time_to_phase
def _mock_manage_parameters(parameter_names, context):
    """Mock function to return dummy fixed values for disc parameters."""
    return [0.1, 10, 100, 0.5, 30, 45, 1.2, 1.5, 0.8, 0.2, 5, 50, 0.2, 15, 35, 2.0, 500, 6.4]  # Replace with reasonable defaults

def _mock_orbital_time_to_phase(t, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01):
    """Mock function to simulate phase calculation for the disc."""
    x = np.linspace(0, 1, len(t))  # Simulate orbital phase values
    W = 1.0  # A simple placeholder for W
    return x, None, W

# Test for doppler_disc_theoretical with keV units
def test_doppler_disc_theoretical_keV(monkeypatch):
    monkeypatch.setattr('xraybinaryorbit._manage_parameters', _mock_manage_parameters)
    monkeypatch.setattr('xraybinaryorbit._orbital_time_to_phase', _mock_orbital_time_to_phase)

    t = np.arange(0, 1000, 1)
    t_out, x_out, x2_out, equation_out = doppler_disc_theoretical(t, units="keV")

    assert len(t_out) == len(t), "The output time array should match the input length."
    assert len(x_out) == len(t), "The output orbital phase array for the first orbit should match the input length."
    assert len(x2_out) == len(t), "The output orbital phase array for the second orbit should match the input length."
    assert len(equation_out) == len(t), "The Doppler variation array should match the input length."
    assert isinstance(equation_out, np.ndarray), "Doppler variation should be a NumPy array."

# Test for plot generation with keV units
def test_doppler_disc_theoretical_plot_keV(monkeypatch):
    monkeypatch.setattr('xraybinaryorbit._manage_parameters', _mock_manage_parameters)
    monkeypatch.setattr('xraybinaryorbit._orbital_time_to_phase', _mock_orbital_time_to_phase)

    t = np.arange(0, 1000, 1)
    try:
        doppler_disc_theoretical(t, units="keV", show_plot=True)
    except Exception as e:
        pytest.fail(f"Plot generation failed with error: {e}")
