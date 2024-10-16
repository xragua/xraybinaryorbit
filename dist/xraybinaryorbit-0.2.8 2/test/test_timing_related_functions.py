import pytest
from xraybinaryorbit import *

# Fake light curve generation
t = np.arange(0, 3 * 24 * 60 * 60)  # time in seconds over 3 days
c = 3 * np.sin(2 * np.pi * t / (0.2 * 24 * 60 * 60)) + 10  # sinusoidal flux with a baseline
sc = 0.1 * c  # 10% error on flux


def test_hr():
    # Generate fake data for hard and soft bands
    hard_band = c  # Use the same light curve for the hard band
    soft_band = 0.5 * c  # Make the soft band half the value of the hard band
    hard_error = sc  # Error is 10% of the hard band
    soft_error = 0.5 * sc  # Error is also 10% of the soft band
    
    # Calculate hardness ratio and errors
    hr_value, hr_error = hr(hard_band, soft_band, hard_error, soft_error)
    
    # Assertions
    assert len(hr_value) == len(c), "Hardness ratio should have the same length as the input data"
    assert len(hr_error) == len(c), "Hardness ratio error should have the same length as the input data"
    assert np.all(hr_error >= 0), "Errors should be non-negative"


def test_cr():
    # Generate fake data for hard and soft bands
    hard_band = c
    soft_band = 0.5 * c
    hard_error = sc
    soft_error = 0.5 * sc
    
    # Calculate color ratio and errors
    cr_value, cr_error = cr(hard_band, soft_band, hard_error, soft_error)
    
    # Assertions
    assert len(cr_value) == len(c), "Color ratio should have the same length as the input data"
    assert len(cr_error) == len(c), "Color ratio error should have the same length as the input data"
    assert np.all(cr_error >= 0), "Errors should be non-negative"


def test_rebin_snr():
    # Test rebinning with a signal-to-noise ratio threshold
    snr_threshold = 0.2
    t_rebinned, c_rebinned, sc_rebinned = rebin_snr(t, c, sc, snr_threshold)
    
    # Assertions
    assert len(t_rebinned) < len(t), "Rebinned time array should have fewer points than the original"
    assert len(c_rebinned) == len(t_rebinned), "Rebinned flux should match rebinned time"
    assert len(sc_rebinned) == len(t_rebinned), "Rebinned errors should match rebinned time"
    
def test_rebin_bins():
    # Number of bins for rebinning
    nbin = 10

    # Run the rebin_bins function
    t_rebinned, c_rebinned, sc_rebinned = rebin_bins(t, c, sc, nbin)

    # Assert the output lengths are as expected after rebinning
    assert len(t_rebinned) == len(c_rebinned) == len(sc_rebinned), "The rebinned arrays should have the same length."
    assert len(t_rebinned) == len(t) // nbin, "The rebinned arrays should be shorter based on the bin size."

    # Assert the rebinned time array has reasonable values
    assert all(np.diff(t_rebinned) > 0), "The rebinned time array should be monotonically increasing."

    # Assert that the rebinned errors are positive
    assert all(np.array(sc_rebinned) > 0), "The rebinned errors should be positive."

    # Assert that the count rate is rebinned correctly
    assert np.allclose(np.mean(c_rebinned), np.mean(c), atol=0.1), "The mean of the rebinned count rate should be close to the original."

    
def test_fold_pulse():
    # Period of the fake signal (0.2 days converted to seconds)
    period = 0.2 * 24 * 60 * 60
    
    # Fold the pulse
    t_folded, c_folded, sc_folded = fold_pulse(t, c, sc, period, snr=0.2)
    
    # Assertions
    assert len(t_folded) < len(t), "Folded time array should have fewer points than the original"
    assert len(c_folded) == len(t_folded), "Folded flux should match folded time"
    assert len(sc_folded) == len(t_folded), "Folded errors should match folded time"    
    
def test_period_sliding_window():
    # Sliding window parameters
    window_sec = 100  # Window size in seconds
    step_sec = 50  # Step size in seconds
    min_period = 0.1 * 24 * 60 * 60  # Minimum period in seconds
    max_period = 0.3 * 24 * 60 * 60  # Maximum period in seconds

    # Run the sliding window period analysis
    result, pulses = period_sliding_window(t, c, sc, window_sec, step_sec, max_period=max_period, min_period=min_period, folded_pulses=True)

    # Assertions for period analysis
    assert result is not None, "The result should not be None."
    assert not result.empty, "The result DataFrame should not be empty."
    assert 'Period' in result.columns, "The result should contain 'Period' column."
    assert 'Power' in result.columns, "The result should contain 'Power' column."

    # Check for the presence of folded pulses
    assert pulses is not None, "Pulses dictionary should not be None."
    assert len(pulses) > 0, "There should be folded pulses calculated."
