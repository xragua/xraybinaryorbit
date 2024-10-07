import pytest
from xraybinaryorbit import *


# Test for ionization_map_phase with default parameters
def test_ionization_map_phase_default():
    # Run the function with default parameters and loading data directly
    chi_result, area_between_bounds = ionization_map_phase(load_directly=True)

    # Assert outputs
    assert isinstance(chi_result, pd.DataFrame), "The output chi_result should be a pandas DataFrame."
    assert '0.0' in chi_result.columns, "The DataFrame should contain orbital phase columns."
    assert not chi_result.empty, "The DataFrame should not be empty."
    assert isinstance(area_between_bounds, float), "The area between bounds should be a float."
    assert area_between_bounds >= 0, "The area should be a non-negative number."

    # Check if the plot is saved correctly
    plot_filename = "test_ionization_map.png"
    ionization_map_phase(save_plot=True, name="test_ionization_map", load_directly=True)

    # Assert that the plot file has been created
    assert os.path.isfile(plot_filename), "The plot should be saved as a file."

    # Clean up: remove the file after the test
    if os.path.exists(plot_filename):
        os.remove(plot_filename)


# Test for ionization_map_phase with custom size and color scale
def test_ionization_map_phase_custom_size_color():
    # Run the function with custom size and color scale
    chi_result, area_between_bounds = ionization_map_phase(size_in_Rstar=3.0, min_color=1e-5, max_color=1e2, load_directly=True)

    # Assert outputs
    assert isinstance(chi_result, pd.DataFrame), "The output chi_result should be a pandas DataFrame."
    assert not chi_result.empty, "The DataFrame should not be empty."
    assert isinstance(area_between_bounds, float), "The area between bounds should be a float."
    assert area_between_bounds >= 0, "The area should be a non-negative number."

    # Check some specific properties of the custom color scale and size
    assert chi_result.shape[0] > 0, "The DataFrame should have at least one row."
