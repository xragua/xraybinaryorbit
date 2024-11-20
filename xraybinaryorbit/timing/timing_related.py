import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.interpolate import interp1d
import math
import matplotlib.pyplot as plt
from tkinter import messagebox
import tkinter as tk
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks, peak_widths

from ..helpers.data_helpers import _manage_parameters,_define_x_y_sy,_copy_fields, _load_values_to_interface, _manage_parameters,_load_bounds_to_interface, _manage_bounds

from ..helpers.math_helpers import _gaussian,_time_pairs,_interpolate_pchip,_chi_squared_weighted,_chi_squared,_orbital_phase_to_time,_orbital_time_to_phase,scale

c = 299792458

msun = (1.98847*10**30)*1000 #gr
rsun_m = 696340*1000 #
rsun_cm = 696340*1000*100 #cm

kev_ams = 1.23984193

na = 6.02214076*10**23/1.00797
mu = 0.5
mp = 1.67E-24

#.................................................. Hardness ratio
def hr(x, y, ex, ey):
    """
    Calculates the hardness ratio and its associated errors.

    The hardness ratio (HR) is calculated using the formula:
    
    HR = (h - l) / (h + l)

    where:
    - h is the count rate or flux in the hard band.
    - l is the count rate or flux in the soft band.

    Parameters
    ----------
    x : float or array-like
        Count rate or flux in the hard band.
    y : float or array-like
        Count rate or flux in the soft band.
    ex : float or array-like
        Errors associated with the hard band.
    ey : float or array-like
        Errors associated with the soft band.

    Returns
    -------
    hr : float or array-like
        The calculated hardness ratio.
    hr_error : float or array-like
        The error in the hardness ratio.
    """


    f_value = (x - y) / (x + y)
    df_dx = 2 * y / (x + y)**2
    df_dy = -2 * x / (x + y)**2
    
    sigma_f = np.sqrt((df_dx * ex)**2 + (df_dy * ey)**2)
    
    return f_value, sigma_f

#.................................................. Color ratio
def cr(x, y, ex, ey):
    """
    Calculates the color ratio and its associated errors.

    The color ratio (CR) is calculated using the formula:
    
    CR = h / l

    where:
    - h is the count rate or flux in the hard band.
    - l is the count rate or flux in the soft band.

    Parameters
    ----------
    x : float or array-like
        Count rate or flux in the hard band.
    y : float or array-like
        Count rate or flux in the soft band.
    ex : float or array-like
        Errors associated with the hard band.
    ey : float or array-like
        Errors associated with the soft band.

    Returns
    -------
    cr : float or array-like
        The calculated color ratio.
    cr_error : float or array-like
        The error in the color ratio.
    """


    f_value = x / y
    df_dx = 1 / y
    df_dy = -x / (y**2)
    
    sigma_f = np.sqrt((df_dx * ex)**2 + (df_dy * ey)**2)
    
    return f_value, sigma_f


#.................................................. Rebin signal to noise
def rebin_snr(t, x, sy, snr_threshold):
    """
    Calculates a rebinned signal-to-noise ratio (SNR) lightcurve.

    Parameters
    ----------
    t : array-like
        Time array.
    x : array-like
        Count rate or flux array.
    sy : array-like
        Error array associated with the count rate or flux.
    snr_threshold : float
        Minimum required signal-to-noise ratio (typically between 0.2 and 0.05).

    Returns
    -------
    t_rebinned : array-like
        Rebinned time array.
    x_rebinned : array-like
        Rebinned lightcurve.
    sy_rebinned : array-like
        Rebinned errors.
    """

    
    w=[]
    c_bin=[]
    t_bin=[]
    sc_bin=[]

    c_new=[]
    t_new=[]
    sc_new=[]
    
    mask = np.where(sy > 0)[0]
    
    t=t[mask]
    x=x[mask]
    sy=sy[mask]
    
    for i in range(len(x)-1):

        w.append(pow(1/sy[i],2))
        t_bin.append(t[i])
        c_bin.append(x[i])
        sc_bin.append(sy[i])
        
        sc_weight = pow(1/(sum(np.array(w))),0.5)
        c_weight = sum(np.array(c_bin)*np.array(w))/sum(np.array(w))
        
        snr_now = sc_weight/c_weight

        if (snr_now <= snr_threshold):

            w = np.array(w)
            c_bin = np.array(c_bin)
            sc_bin = np.array(sc_bin)
            t_bin = np.array(t_bin)

            sc_new.append(pow(1/(sum(w)),0.5))
            c_new.append(sum(c_bin*w)/sum(w))
            t_new.append(sum(t_bin*w)/sum(w))

            w=[]
            c_bin=[]
            t_bin=[]
            sc_bin=[]
     
        
    return t_new,c_new,sc_new

#.................................................. Rebin by bins
def rebin_bins(t, c, sc, nbin):
    """
    Calculates a rebinned lightcurve using a specified number of time bins.

    Parameters
    ----------
    t : array-like
        Time array.
    x : array-like
        Count rate or flux array.
    sy : array-like
        Error array associated with the count rate or flux.
    nbin : int
        Number of time bins for rebinning (requires a larger time bin compared to the original one).

    Returns
    -------
    t_rebinned : array-like
        Rebinned time array.
    x_rebinned : array-like
        Rebinned lightcurve.
    sy_rebinned : array-like
        Rebinned errors.
    """


    c_new=[]
    t_new=[]
    sc_new=[]
    
    t, c, sc = preprocess_data(t, c, sc)
    
    for i in range(len(c)-nbin-1):

        w = (pow(1/sc[i*nbin:(i+1)*nbin],2))
        
        if (sum(w) >0):
        
            t_bin = t[i*nbin:(i+1)*nbin]
            c_bin = c[i*nbin:(i+1)*nbin]
            sc_bin = sc[i*nbin:(i+1)*nbin]

            #...............................

            sc_new.append(pow(1/(sum(w)),0.5))
            c_new.append(sum(c_bin*w)/sum(w))
            t_new.append(sum(t_bin*w)/sum(w))

    return t_new,c_new,sc_new
#.................................................. Periods sliding window
def fold_pulse(t, c, sc, period, snr=None, rebin=None):
    """
    Folds a lightcurve data array based on a given period and optionally rebins it using either signal-to-noise
    ratio (SNR) or uniform time binning.

    Parameters
    ----------
    t : array-like
        Time array of the lightcurve.
    c : array-like
        Count rate or flux array corresponding to the time array.
    sc : array-like
        Errors (standard deviation) associated with the count rate or flux.
    period : float
        Period to fold the lightcurve, in the same units as the time array `t`.

    Optional Parameters
    -------------------
    snr : float, optional
        Minimum signal-to-noise ratio threshold. If provided, applies signal-to-noise rebinning using `rebin_snr`.
    rebin : int, optional
        Number of bins to rebin the folded lightcurve. If provided, applies uniform time binning using `rebin_bins`.

    Returns
    -------
    If `snr` is specified:
        t_rebinned : array-like
            Rebinned time array after signal-to-noise ratio rebinning.
        c_rebinned : array-like
            Rebinned folded lightcurve.
        sc_rebinned : array-like
            Rebinned errors after signal-to-noise ratio rebinning.

    If `rebin` is specified:
        t_rebinned : array-like
            Rebinned time array after uniform time binning.
        c_rebinned : array-like
            Rebinned folded lightcurve.
        sc_rebinned : array-like
            Rebinned errors after uniform time binning.

    Notes
    -----
    Either `snr` or `rebin` must be specified to proceed with the function.
    """

    
    phase = (t - min(t)) / period - np.floor((t - min(t)) / period)

    idx_sort = np.argsort(phase)
    
    phase_sorted = phase[idx_sort]
    c_sorted = c[idx_sort]
    sc_sorted = sc[idx_sort]
    
    if snr:
        pulse = rebin_snr(phase_sorted, c_sorted, sc_sorted, snr)
        return pulse
        
    if rebin:
        pulse = rebin_bins(phase_sorted, c_sorted, sc_sorted, rebin)
        return pulse
        
    else:
        print("Please provide a Minimum signal-to-noise ratio threshold (snr) or a bin time (rebin) to proceed.")



def preprocess_data(t, x, sy):
    """ Preprocess the data by removing entries where sc <= 0, NaN, or inf. """
    mask = (sy > 0) & np.isfinite(sy)
    t = t[mask]
    x = x[mask]
    sy = sy[mask]
    
    return t, x, sy
#.................................................. Periods sliding window


def period_sliding_window(t, c, sc, window_sec, step_sec, max_period=None, min_period=None, false_alarm_threshold=0.1,
                          rel_high_for_error=0.9, folded_pulses=False, snr_pulse=0.2, nbin_pulse=None):
    """
    Performs period analysis using a sliding window approach on a lightcurve dataset.

    Parameters
    ----------
    t : array-like
        Time array of the lightcurve.
    c : array-like
        Count rate or flux array corresponding to the time array.
    sc : array-like
        Errors (standard deviation) associated with the count rate or flux.
    window_sec : float
        Size of the sliding window in seconds for period analysis.
    step_sec : float
        Step size in seconds between consecutive windows.
    max_period : float, optional
        Maximum period to consider in the periodogram analysis. Default is None.
    min_period : float, optional
        Minimum period to consider in the periodogram analysis. Default is None.
    false_alarm_threshold : float, optional
        Threshold value for false alarm probability in Lomb-Scargle periodogram. Default is 0.1.
    rel_high_for_error : float, optional
        Relative height for error estimation in the `peak_widths` function. Default is 0.9.
    folded_pulses : bool, optional
        If True, folds the lightcurve using `fold_pulse` for each identified period. Default is False.
    snr_pulse : float, optional
        Minimum signal-to-noise ratio threshold for folding using `fold_pulse`. Default is 0.2.
    nbin_pulse : int, optional
        Number of bins for uniform time binning in folding using `fold_pulse`. Default is None.

    Returns
    -------
    result : DataFrame
        DataFrame containing the results of the period analysis, including:
        - periods,
        - frequencies,
        - powers,
        - errors in period,
        - errors in power,
        - false alarm probabilities,
        - time range.
    pulses : dict, optional
        Dictionary containing folded pulse data for each identified period if `folded_pulses` is True.

    Notes
    -----
    - For each window, the function returns a list of all peaks found in the periodogram.
    - The function performs Lomb-Scargle periodogram analysis within each sliding window of the specified size.
    - It filters periods based on false alarm probability and sorts them by power.
    - If `folded_pulses` is True, it folds the lightcurve for each identified period using `fold_pulse`
      and stores the results.
    """


    def preprocess_data(t, c, sc):
        """ Preprocess the data by removing entries where sc <= 0, NaN, or inf. """
        mask = (sc > 0) & np.isfinite(sc)
        return t[mask], c[mask], sc[mask]

    def lb_period_freq(t, c, sc, step, window, max_period, min_period, false_alarm_threshold):
    
        min_freq = 1 / max_period if max_period else None
        max_freq = 1 / min_period if min_period else None

        # Lomb-Scargle periodogram
        freq, power = LombScargle(t, c, sc).autopower(maximum_frequency=max_freq, minimum_frequency=min_freq, samples_per_peak=1000)
        ls = LombScargle(t, c, sc)
        
        pos = find_peaks(power)[0]  # Peak positions
        df_list = []

        if len(pos) > 0:
            # False alarm probability for the peaks
            fa_prob = ls.false_alarm_probability(power[pos])
            index_fa = fa_prob < false_alarm_threshold
            peaks_indices = pos[index_fa]

            if len(peaks_indices) > 0:
            
                # Calculate widths and errors for the peaks
                results_widths = peak_widths(power, peaks_indices, rel_height=rel_high_for_error)
                
                # Loop through each peak and store data
                for i, index in enumerate(peaks_indices):
                    
                    # Ensure that the frequency slice contains enough elements
                    start_idx = max(int(results_widths[2][i]) - 3, 0)
                    end_idx = min(int(results_widths[2][i]) + 3, len(freq) - 1)
                    
                    if end_idx > start_idx:  # Ensure there are at least two elements
                        freq_error = max(np.diff(freq[start_idx:end_idx])) * 3
                    else:
                        freq_error = 0  # Assign 0 if there's not enough data to calculate frequency error
                
                    power_error = results_widths[1][i]  # Width corresponds to the peak
                    snr=power[index]/np.median(power)

                    df_list.append({
                        'min_time': min(t),
                        'max_time': max(t),
                        'Frequency': freq[index],
                        'Period': 1 / freq[index],
                        'Power': power[index],
                        'Freq_Error': freq_error,  # Properly indexed freq_error
                        'Period_Error': (freq_error / freq[index]**2)**(1/2),  # Error in period
                        'Power_Error': power_error,  # Corresponding power error
                        'False_alarm': fa_prob[index_fa][i],  # False alarm probability
                        'snr': snr
                    })

        # Return DataFrame with results
        df = pd.DataFrame(df_list).reset_index(drop=True)

        if len(df) > 1:
            df = df.sort_values(by='Power', ascending=False).reset_index(drop=True)

        return df

    # Preprocess the data
    t = np.asarray(t)
    c = np.asarray(c)
    sc = np.asarray(sc)
    
    t, c, sc = preprocess_data(t, c, sc)

    # Collect results using sliding window
    periods = {}
    
    
    min_period_ = 2 * np.mean(np.diff(t))
    
    if (min_period_ >= min_period):

        print(f"\033[91mWarning: The minimum period given is too low to be accurately captured with the current sampling rate. According to Nyquist-Shannon Sampling Theorem, the minimum period allowed to be identified is {min_period_:.4f} seconds.\033[0m")

    
    for i in range(0, len(t) - window_sec - 1, step_sec):
        t_window = t[i:i + window_sec]
        c_window = c[i:i + window_sec]
        sc_window = sc[i:i + window_sec] + 1e-9  # Avoid division by zero

        if len(t_window) == 0:
            continue

        df_periods = lb_period_freq(t_window, c_window, sc_window, step_sec, window_sec, max_period, min_period, false_alarm_threshold).reset_index(drop=True)

        if len(df_periods) > 0:
            periods[i] = df_periods

    # Combine results
    if len(periods) >= 1:
        result = pd.concat(periods).reset_index(drop=True)
    else:
        result = None

    # Fold pulses if requested
    pulses = {}
    if folded_pulses and result is not None:
        for i in range(len(result)):
            idx = (t >= result.min_time[i]) & (t <= result.max_time[i])

            t_ = np.array(t[idx])
            c_ = np.array(c[idx])
            sc_ = np.array(sc[idx])

            ph_pulse, c_pulse, sc_pulse = fold_pulse(t_, c_, sc_, result.Period[i], snr=snr_pulse, rebin=nbin_pulse)
            pulse_data = {'ph_pulse': ph_pulse, 'c_pulse': c_pulse, 'sc_pulse': sc_pulse}

            pulses[i] = pulse_data

    return result, pulses
