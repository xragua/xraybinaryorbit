import pytest
import numpy as np
from xraybinaryorbit import doppler_spiral_theoretical

# Mock functions for _manage_parameters and _orbital_time_to_phase
def _mock_manage_parameters(parameter_names, context):
    """Mock function to return dummy fixed values for parameters."""
    return [0.1, 10, 0.05, 2 * np.pi, 45, 6.4]  # Example values for spiral parameters

def _mock_orbital_time_to_phase(t, iphase, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2, precision=0.01):
    """Mock function to simulate phase calculation."""
    x = np.linspace(0, 1, len(t))  # Simulate orbital phase values
    W = 1.0  # A simple placeholder for W
    return x, None, W

# Test for doppler_spiral_theoretical with keV units
def test_doppler_spiral_theoretical_keV(monkeypatch):
    monkeypatch.setattr('xraybinaryorbit._manage_parameters', _mock_manage_parameters)
    monkeypatch.setattr('xraybinaryorbit._orbital_time_to_phase', _mock_orbital_time_to_phase)

    t = np.arange(0, 1000, 1)
    t_out, x_out, equation_out = doppler_spiral_theoretical(t, units="keV")

    assert len(t_out) == len(t), "The output time array should match the input length."
    assert len(x_out) == len(t), "The output orbital phase array should match the input length."
    assert len(equation_out) == len(t), "The Doppler variation array should match the input length."
    assert isinstance(equation_out, np.ndarray), "Doppler variation should be a NumPy array."

# Test for invalid unit handling
def test_doppler_spiral_theoretical_invalid_units(monkeypatch):
    monkeypatch.setattr('xraybinaryorbit._manage_parameters', _mock_manage_parameters)
    monkeypatch.setattr('xraybinaryorbit._orbital_time_to_phase', _mock_orbital_time_to_phase)

    t = np.arange(0, 1000, 1)
    with pytest.raises(KeyError):
        doppler_spiral_theoretical(t, units="invalid_unit")

# Test for plot generation with keV units
def test_doppler_spiral_theoretical_plot_keV(monkeypatch):
    monkeypatch.setattr('xraybinaryorbit._manage_parameters', _mock_manage_parameters)
    monkeypatch.setattr('xraybinaryorbit._orbital_time_to_phase', _mock_orbital_time_to_phase)

    t = np.arange(0, 1000, 1)
    try:
        doppler_spiral_theoretical(t, units="keV", show_plot=True)
    except Exception as e:
        pytest.fail(f"Plot generation failed with error: {e}")
