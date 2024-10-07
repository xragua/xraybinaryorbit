import pytest
import numpy as np
import sys

sys.path.append('/Users/graci/Desktop/git/xraybinaryorbit/xraybinaryorbit')
from xraybinaryorbit import *


x = np.array([[7.97000557e+08, 7.97002557e+08],
       [7.97004557e+08, 7.97006557e+08]])

y = np.array([6.68379, 6.80991])

y_err = np.array([[0.01962, 0.01599],[0.09961, 0.01009]])



###############################PS


def test_fit_orbit_ps():
    # Run the function
    result, phase, predicted_data, r_squared = fit_orbit_ps(x, y, y_err=y_err, num_iterations=1, maxiter=100, swarmsize=10, units="keV", method_="extended", extended_binsize=0.01, load_directly=True)

    # Assert the result is a pandas DataFrame
    assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame."
    assert len(phase) == len(x), "Phase should match input length."
    assert len(predicted_data) == len(x), "Predicted data should match input length."
    assert isinstance(r_squared, float), "R-squared should be a float."
    
#def test_fit_spiral_in_orbit_ps():
    # Run the function
 #   result, phase, predicted_data, r_squared = fit_spiral_in_orbit_ps(x, y, y_err=y_err, num_iterations=1, maxiter=100, swarmsize=10, units="keV", method_="extended", extended_binsize=0.01, load_directly=True)

    # Assert the result is a pandas DataFrame
   # assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame."
   # assert len(phase) == len(x), "Phase should match input length."
    #assert len(predicted_data) == len(x), "Predicted data should match input length."
    #assert isinstance(r_squared, float), "R-squared should be a float."
    
#def test_fit_spiral_ps():
    # Run the function
 #   result, phase, predicted_data, r_squared = fit_spiral_ps(x, y, y_err=y_err, num_iterations=1, maxiter=100, swarmsize=10, units="keV", method_="extended", extended_binsize=0.01, load_directly=True)

    # Assert the result is a pandas DataFrame
   # assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame."
  #  assert len(phase) == len(x), "Phase should match input length."
   # assert len(predicted_data) == len(x), "Predicted data should match input length."
   # assert isinstance(r_squared, float), "R-squared should be a float."
    
#def test_fit_disc_ps():
    # Run the function
#    result, phase,dphase, predicted_data, r_squared = fit_disc_ps(x, y, y_err=y_err, num_iterations=1, maxiter=100, swarmsize=10, units="keV", method_="extended", extended_binsize=0.01, load_directly=True)

    # Assert the result is a pandas DataFrame
 #   assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame."
 #   assert len(phase) == len(x), "Phase should match input length."
 #   assert len(predicted_data) == len(x), "Predicted data should match input length."
#    assert isinstance(r_squared, float), "R-squared should be a float."


###############################LS


def test_fit_orbit_ls():
    # Run the function
    result, phase, predicted_data, r_squared = fit_orbit_ls(x, y, y_err=y_err, units="keV", method_="extended", extended_binsize=0.01, load_directly=True)

    # Assert the result is a pandas DataFrame
    assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame."
    assert len(phase) == len(x), "Phase should match input length."
    assert len(predicted_data) == len(x), "Predicted data should match input length."
    assert isinstance(r_squared, float), "R-squared should be a float."
    
#def test_fit_spiral_in_orbit_ls():
    # Run the function
 #   result, phase, predicted_data, r_squared = fit_spiral_in_orbit_ls(x, y, y_err=y_err, units="keV", method_="extended", extended_binsize=0.01, load_directly=True)

    # Assert the result is a pandas DataFrame
 #   assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame."
 #   assert len(phase) == len(x), "Phase should match input length."
 #   assert len(predicted_data) == len(x), "Predicted data should match input length."
 #   assert isinstance(r_squared, float), "R-squared should be a float."
    
#def test_fit_spiral_ls():
    # Run the function
#    result, phase, predicted_data, r_squared = fit_spiral_ls(x, y, y_err=y_err,  units="keV", method_="extended", extended_binsize=0.01, load_directly=True)

    # Assert the result is a pandas DataFrame
 #   assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame."
 #   assert len(phase) == len(x), "Phase should match input length."
 #   assert len(predicted_data) == len(x), "Predicted data should match input length."
  #  assert isinstance(r_squared, float), "R-squared should be a float."
    
def test_fit_disc_ls():
    # Run the function
    result, phase,dphase, predicted_data, r_squared = fit_disc_ls(x, y, y_err=y_err, units="keV", method_="extended", extended_binsize=0.01, load_directly=True)

    # Assert the result is a pandas DataFrame
    assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame."
    assert len(phase) == len(x), "Phase should match input length."
    assert len(predicted_data) == len(x), "Predicted data should match input length."
    assert isinstance(r_squared, float), "R-squared should be a float."


ti= np.array([ 7.39011885e+08, 7.69450823e+08])
te = np.array([7.39043275e+08,7.69480053e+08])
nh1=([1.1215199999999999e+22, 1.27078e+22])
enh1 = np.array([1, 1])

x_data = np.transpose([ti, te])
y_data = nh1
y_err= enh1

def test_fit_nh():
    # Run the function
    result, phase,predicted_data, r_squared = fit_nh_ps(x_data, y_data, num_iterations=1, maxiter=5, swarmsize=2, method_="extended", extended_binsize=0.1,load_directly=True)

    # Assert the result is a pandas DataFrame
    assert isinstance(result, pd.DataFrame), "Result should be a pandas DataFrame."
    assert len(phase) == len(x), "Phase should match input length."
    assert len(predicted_data) == len(x), "Predicted data should match input length."
    assert isinstance(r_squared, float), "R-squared should be a float."
