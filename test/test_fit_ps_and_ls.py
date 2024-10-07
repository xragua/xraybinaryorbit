import pytest
from xraybinaryorbit import *
#from xraybinaryorbit import _manage_parameters, _orbital_phase_to_time, _orbital_time_to_phase, _manage_bounds

# Test data
x = np.array([[7.97000557e+08, 7.97002557e+08], [7.97004557e+08, 7.97006557e+08]])
y = np.array([6.68379, 6.80991])
y_err = np.array([[0.01962, 0.01599], [0.09961, 0.01009]])

def test_fit_orbit_ps():
    # Run the function
    result, phase, predicted_data, r_squared = fit_orbit_ps(
        x, y, y_err=y_err, num_iterations=1, maxiter=100,
        swarmsize=10, units="keV", method_="extended",
        extended_binsize=0.01, load_directly=True
    )

    # Assertions
    assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame."
    assert len(phase) == len(x), "Phase should match input length."
    assert len(predicted_data) == len(x), "Predicted data should match input length."
    assert isinstance(r_squared, float), "R-squared should be a float."

def test_fit_spiral_in_orbit_ps():
    # Run the function
    result, phase, predicted_data, r_squared = fit_spiral_in_orbit_ps(
        x, y, y_err=y_err, num_iterations=1, maxiter=100,
        swarmsize=10, units="keV", method_="extended",
        extended_binsize=0.01, load_directly=True)

    # Assertions
    assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame."
    assert len(phase) == len(x), "Phase should match input length."
    assert len(predicted_data) == len(x), "Predicted data should match input length."
    assert isinstance(r_squared, float), "R-squared should be a float."

def test_fit_spiral_ps():
    # Run the function
    result, phase, predicted_data, r_squared = fit_spiral_ps(
        x, y, y_err=y_err, num_iterations=1, maxiter=100,
        swarmsize=10, units="keV", method_="extended",
        extended_binsize=0.01, load_directly=True
    )

    # Assertions
    assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame."
    assert len(phase) == len(x), "Phase should match input length."
    assert len(predicted_data) == len(x), "Predicted data should match input length."
    assert isinstance(r_squared, float), "R-squared should be a float."

def test_fit_disc_ps():
    # Run the function
    result, phase, dphase, predicted_data, r_squared = fit_disc_ps(
        x, y, y_err=y_err, num_iterations=1, maxiter=100,
        swarmsize=10, units="keV", method_="extended",
        extended_binsize=0.01, load_directly=True
    )
    
    # Assertions
    assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame."
    assert len(phase) == len(x), "Phase should match input length."
    assert len(predicted_data) == len(x), "Predicted data should match input length."
    assert isinstance(r_squared, float), "R-squared should be a float."

# Additional test data for NH fitting
ti = np.array([7.39011885e+08, 7.69450823e+08])
te = np.array([7.39043275e+08, 7.69480053e+08])
nh1 = [1.1215199999999999e+22, 1.27078e+22]
enh1 = np.array([1, 1])

x_data = np.transpose([ti, te])
y_data = nh1
y_err = enh1

def test_fit_nh():
    # Run the function
    result, phase, predicted_data, r_squared = fit_nh_ps(
        x_data, y_data, num_iterations=1, maxiter=2,
        swarmsize=2, method_="extended", extended_binsize=0.1,
        load_directly=True
    )

    # Assertions
    assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame."
    assert len(phase) == len(x_data), "Phase should match input length."
    assert len(predicted_data) == len(x_data), "Predicted data should match input length."
    assert isinstance(r_squared, float), "R-squared should be a float."
