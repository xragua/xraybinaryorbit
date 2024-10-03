import pytest
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from xraybinaryorbit import ionization_map_phase

# Mock functions for _manage_parameters
def _mock_manage_parameters(parameter_names, context):
    """Mock function to return dummy fixed values for parameters."""
    return [0.5, 10, 0.1, 30, 1.0, 1.5, 0.8, 1000, 2e-5, 0.8, 1e38, 1e3, 1e5]  # Example values for parameters

# Test for ionization_map_phase with default parameters
def test_ionization_map_phase_default(monkeypatch):
    # Patch the internal functions
    monkeypatch.setattr('xraybinaryorbit._manage_parameters', _mock_manage_parameters)

    # Run the function with default parameters
    chi_result, area_between_bounds = ionization_map_phase()

    # Assert outputs
    assert isinstance(chi_result, pd.DataFrame), "The output chi_result should be a pandas DataFrame."
    assert '0.0' in chi_result.columns, "The DataFrame should contain orbital phase columns."  # Updated
    assert not chi_result.empty, "The DataFrame should not be empty."
    assert isinstance(area_between_bounds, float), "The area between bounds should be a float."
    assert area_between_bounds >= 0, "The area should be a non-negative number."  # Updated

# Test for plot generation
def test_ionization_map_phase_plot(monkeypatch):
    # Patch the internal functions
    monkeypatch.setattr('xraybinaryorbit._manage_parameters', _mock_manage_parameters)

    # Run the function with save_plot=True (it should not raise any exceptions)
    try:
        ionization_map_phase(save_plot=True, name="test_ionization_map")
    except Exception as e:
        pytest.fail(f"Plot generation failed with error: {e}")

# Test for ionization map with custom size and color scale
def test_ionization_map_phase_custom_size_color(monkeypatch):
    # Patch the internal functions
    monkeypatch.setattr('xraybinaryorbit._manage_parameters', _mock_manage_parameters)

    # Run the function with custom size and color scale
    chi_result, area_between_bounds = ionization_map_phase(size_in_Rstar=3.0, min_color=1e-5, max_color=1e2)

    # Assert outputs
    assert isinstance(chi_result, pd.DataFrame), "The output chi_result should be a pandas DataFrame."
    assert not chi_result.empty, "The DataFrame should not be empty."
    assert isinstance(area_between_bounds, float), "The area between bounds should be a float."
    assert area_between_bounds >= 0, "The area should be a non-negative number."  # Updated
