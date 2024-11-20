import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import tkinter as tk
from tkinter import messagebox


from ..helpers.data_helpers import _manage_parameters,_define_x_y_sy,_copy_fields, _load_values_to_interface, _manage_parameters,_load_bounds_to_interface, _manage_bounds

from ..helpers.math_helpers import _gaussian,_time_pairs,_interpolate_pchip,_chi_squared_weighted,_chi_squared,_orbital_phase_to_time,_orbital_time_to_phase

c = 299792458

msun = (1.98847*10**30)*1000 #gr
rsun_m = 696340*1000 #
rsun_cm = 696340*1000*100 #cm

kev_ams = 1.23984193

na = 6.02214076*10**23/1.00797
mu = 0.5
mp = 1.67E-24


# DENSITY IN THE ORBIT #############################################################################
def density_through_orbit_theoretical(resolution=0.01, show_plot=False, load_directly=False, parameter_list=None):
    """
    Visualizes the density (gr/cm^2) encountered by a compact object along its orbit, assuming a spherically
    distributed stellar wind based on the CAK (Castor-Abbott-Klein) model.

    Parameters
    ----------
    resolution : float, optional
        Resolution for the phase array. Default is 0.01.
    show_plot : bool, optional
        If True, displays and saves a plot of the density through the orbit. The plot is saved as
        "density_through_the_orbit.png." Default is False.

    Notes
    -----
    A form will appear to input the necessary orbital parameters. These parameters are saved in a .txt file
    in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter parameters,
    allowing modification of only those that require adjustment.
    By setting load_directly=True the data will be authomatically loaded into the function without rising a form. By providing parameter_list=(list of parameters) the parameters will be trated as an imput and saved within the current directory for subsequent runs.

    Returns
    -------
    time : array-like
        The time array corresponding to the orbital movement.
    phase : array-like
        The orbital phase array.
    density : array-like
        The density encountered by the compact object through the orbit, measured in gr/cm^2.

    """

    parameter_names = ["semimajor","orbitalperiod" ,"eccentricity", "periapsis", "Rstar","Mstar1","Mstar2","wind_infinite_velocity","Mass_loss_rate","beta" ]
    
    fixed_values = _manage_parameters(parameter_names, "density_through_orbit",load_directly=load_directly,parameter_list=parameter_list )
    semimajor,orbitalperiod , eccentricity,periapsis, Rstar, Mstar1, Mstar2, wind_infinite_velocity, Mass_loss_rate, beta = fixed_values

    th = np.arange(0,1,resolution)

    _,time,_ = _orbital_phase_to_time(th,0, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,precision=0.01)

    #........................................................
    abar = semimajor * max(Mstar1,Mstar2)/(Mstar1+Mstar2)
    
    M_dot_grams = Mass_loss_rate * msun/(365*24*60*60) #Mdot gr/s
    vinf_cm_s = wind_infinite_velocity*100000 #cm/s
    Rstar_cm = Rstar*rsun_cm #R* in cm
    
    Rorb = (abar*(1-eccentricity**2)/(1+eccentricity*np.cos((th-periapsis/360)*2*np.pi)))* Rstar * rsun_cm #In cm

    v = vinf_cm_s*(1-Rstar_cm/Rorb)**beta
    
    rho = (M_dot_grams/(4 * np.pi * v[Rorb>Rstar_cm] * Rorb[Rorb>Rstar_cm]**2))
    
    #...................................................
    if show_plot:
    
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

        ax1.errorbar(time[Rorb>Rstar_cm], rho,  label='Data', color='b')
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("Density through the orbit gr/cm$^2$")
        ax1.legend()

        ax2.errorbar(th[Rorb>Rstar_cm], rho,  label='Data', color='b')
        ax2.set_xlabel("Orbital phase")
        ax2.set_ylabel("Density through the orbit gr/cm$^2$")
        ax2.legend()

        ax3 = plt.subplot(1, 3, 3, projection='polar')
        ax3.plot(th[Rorb>0] * 2 * np.pi, Rorb[Rorb>0],"b")
        ax3.set_theta_zero_location("N")
        ax3.set_theta_direction(-1)
        ax3.set_rlabel_position(90)

        plt.tight_layout()
        
        plt.savefig("density_through_the_orbit.png")
    #.....................................................

    return  time[Rorb > Rstar_cm], th[Rorb > Rstar_cm],  rho

# ABSOPTION COLUMN #############################################################################
def absorption_column_through_orbit_theoretical(resolution=0.01, show_plot=True, load_directly=False, parameter_list=None):
    """
    Visualizes the column density (NH1, x 10^22 cm^-2) encountered by radiation emitted at each orbital
    phase as it travels towards an observer. Assumes a spherically distributed, neutral (unionized)
    stellar wind based on the CAK (Castor-Abbott-Klein) model.

    Parameters
    ----------
    resolution : float, optional
        Resolution for the phase array. Default is 0.01.
    show_plot : bool, optional
        If True, displays and saves a plot of the absorption column density through the orbit.
        Default is False.

    Notes
    -----
    A form will appear to input the necessary orbital parameters. These parameters are saved in a .txt file
    in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter
    parameters, allowing modification of only those that require adjustment.
    By setting load_directly=True the data will be authomatically loaded into the function without rising a form. By providing parameter_list=(list of parameters) the parameters will be trated as an imput and saved within the current directory for subsequent runs.
    If the distance to the star is smaller than the stellar radius, the result will be 0.

    Returns
    -------
    time : array-like
        Time array corresponding to the orbital movement.
    phase : array-like
        Orbital phase array.
    NH1 : array-like
        Absorption column density (NH1, x 10^22 cm^-2) through the orbit.
    
    """


    parameter_names = ["semimajor","orbitalperiod" ,"eccentricity", "periapsis" ,"inclination", "Rstar","Mstar1","Mstar2","wind_infinite_velocity","Mass_loss_rate","beta" ]
    
    fixed_values = _manage_parameters(parameter_names, "absorption_column_through_orbit",load_directly=load_directly,parameter_list=parameter_list )
    semimajor, orbitalperiod,eccentricity, periapsis, inclination, Rstar, Mstar1, Mstar2, wind_infinite_velocity, Mass_loss_rate, beta = fixed_values

    th = np.arange(0,1,resolution)
    
    _,time,_ = _orbital_phase_to_time(th,0, semimajor, orbitalperiod, eccentricity, periapsis, Rstar, Mstar1, Mstar2,precision=0.01)

    #........................................................
    abar = semimajor * max(Mstar1,Mstar2)/(Mstar1+Mstar2)
    
    M_dot_grams = Mass_loss_rate * msun/(365*24*60*60) #Mdot gr/s
    vinf_cm_s = wind_infinite_velocity*100000 #cm/s
    Rstar_cm = Rstar*rsun_cm

    Rorb_plot = (abar*(1-eccentricity**2)/(1+eccentricity*np.cos((th-periapsis/360)*2*np.pi)))*Rstar*rsun_cm #In cm
    
    nh=[]

    for i in range(len(th)):
        
        Rorb = (abar*(1-eccentricity**2)/(1+eccentricity*np.cos((th[i]-periapsis/360)*2*np.pi)))*Rstar*rsun_cm #In cm

        def integrand(z):
        
            alpha = np.arccos(np.cos(th[i]*2*np.pi)*np.sin(inclination*2*np.pi/360))
            x = np.sqrt(Rorb**2+z**2-2*Rorb*z*np.cos(alpha))
            v = (vinf_cm_s*(1-Rstar_cm/x)**beta)
            rho = (M_dot_grams/(4 * np.pi * v * x**2))
            return rho

        ne, _ = quad(integrand, 1, Rstar_cm * 1000)

        nh.append(ne*na/1e22)
        
    nh = np.nan_to_num(nh)
    #...................................................
    if show_plot:
    
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10, 5))
        
        ax1.errorbar(time, nh,  label='Data', color='b')
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel("NH1 (x 10$^{22}$ cm$^{-2}$)")
        ax1.legend()
        
        ax2.errorbar(th, nh,  label='Data', color='b')
        ax2.set_xlabel("Orbital phase")
        ax2.set_ylabel("NH1 (x 10$^{22}$ cm$^{-2}$)")
        ax2.legend()

        ax3 = plt.subplot(1, 3, 3, projection='polar')
        ax3.plot(th * 2 * np.pi, Rorb_plot,"b")
        ax3.set_theta_zero_location("N")
        ax3.set_theta_direction(-1)
        ax3.set_rlabel_position(90)

        plt.tight_layout()
        
        plt.savefig("NH_through_the_orbit.png")
    #...................................................
    
    return time, th, nh
         
         

# DENSITY AND LOGCHI #############################################################################

def density_and_ionization_orbital_phase_theoretical(resolution=0.01, size=10, show_plot=True, load_directly=False, parameter_list=None):
    """
    Calculates and visualizes the density and ionization parameter (log(ξ)) encountered by radiation emitted
    at each orbital phase as it travels towards an observer. Assumes a spherically distributed, neutral stellar
    wind based on the CAK (Castor, Abbott, Klein) model. The density profile and ionization parameter are calculated
    along the path from a neutron star (NS) through the stellar wind of its companion.

    Parameters
    ----------
    resolution : float, optional
        The resolution for the phase array (orbital phase). Default is 0.01.
    size : float, optional
        The scaling factor for the size of the path from the NS through the stellar wind. Determines how far into
        the wind the path is extended. Default is 10.
    show_plot : bool, optional
        If True, displays and saves plots of the density, ionization parameter, and orbital path. Default is True.
    load_directly : bool, optional
        If True, attempts to load previously saved orbital parameters from a file. If False, prompts the user to
        input the parameters. Default is False.
By providing parameter_list=(list of parameters) the parameters will be trated as an imput and saved within the current directory for subsequent runs.
        

    Returns
    -------
    z : array-like
        Array of distances along the path from the neutron star to the observer (in cm).
    density : array-like
        Density profile of the stellar wind (in cm$^{-3}$) along the path.
    chi : array-like
        Ionization parameter (log(ξ)) calculated at each point along the path.

    """

    parameter_names = [ "orb_phase", "luminosity","semimajor", "eccentricity", "periapsis", "inclination","Rstar", "Mstar1", "Mstar2", "wind_infinite_velocity", "Mass_loss_rate", "beta"]
    
    fixed_values = _manage_parameters(parameter_names, "den_chi_orbphase",load_directly=load_directly,parameter_list=parameter_list )
    orb_phase, luminosity, semimajor, eccentricity, periapsis, inclination ,Rstar, Mstar1, Mstar2, wind_infinite_velocity, Mass_loss_rate, beta = fixed_values

    th = np.arange(0, 1, resolution)
    abar = semimajor * max(Mstar1, Mstar2) / (Mstar1 + Mstar2)#BARICENTER CORRECTION
    
    luminosity_ = luminosity * 1e+32
    M_dot_grams = Mass_loss_rate * msun / (365 * 24 * 60 * 60)  # Mass loss rate in grams/s
    vinf_cm_s = wind_infinite_velocity * 100000  # Wind velocity in cm/s
    Rstar_cm = Rstar * rsun_cm  # Convert stellar radius to cm

    # Calculate orbital radius
    Rorb_plot = (abar * (1 - eccentricity**2) / (1 + eccentricity * np.cos((th - periapsis / 360) * 2 * np.pi))) * Rstar * rsun_cm  # In cm
    Rorb = (abar * (1 - eccentricity**2) / (1 + eccentricity * np.cos((orb_phase - periapsis / 360) * 2 * np.pi))) * Rstar * rsun_cm  # In cm
    
    # Create path from the NS to some distance towards the observer (z=distrance travelled by the emitted radiation fromm the NS towards the observer)
    z = np.arange(0, Rorb * size, Rorb * size / 10000)  # Range for distance in the wind

    # Calculate the angle `alpha`
    alpha = np.arccos(np.cos(orb_phase * 2 * np.pi) * np.sin(inclination * 2 * np.pi / 360))
    
    # Calculate x (distance from z points to donnor)
    cosalpha = np.round(np.cos(alpha),10)
    x = np.sqrt(Rorb**2 + z**2 - 2 * Rorb * z * cosalpha)
    
    # Velocity and density in the wind depending on distance to the donnor (depending on distance travelled from NS)
    v = (vinf_cm_s * (1 - Rstar_cm / x)**beta)
    density = (M_dot_grams / (4 * np.pi * v * x**2  * 1.67E-24* 0.5))  # Renamed from rho to density (mu, mp)
    
    # Calculate the chi parameter for each z
    chi = np.log(luminosity_ / (density * 4 * np.pi * abs(z)**2))
    #...................................................
    if show_plot:
    
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
        
        ax1.errorbar(z/Rstar_cm, density,  label='Data', color='b')
        ax1.set_xlabel("Path from NS (R*)")
        ax1.set_ylabel("Density cm$^{-3}$)")
        ax1.legend()
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        
        ax2.errorbar(z/Rstar_cm, chi,  label='Data', color='b')
        ax2.set_xlabel("Path from NS (R*)")
        ax2.set_ylabel("Ionization Parameter (log(ξ))")
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.legend()

        ax3 = plt.subplot(1, 3, 3, projection='polar')
        ax3.plot(th * 2 * np.pi, Rorb_plot,"b")
        ax3.set_theta_zero_location("N")
        ax3.set_theta_direction(-1)
        ax3.set_rlabel_position(90)

        plt.tight_layout()
        
        plt.savefig("density_and_chi_orbital_phase.png")
    #...................................................
    
    return z, density, chi
         
         
         
# Ionization parameter map ###############################################################


def ionization_map_phase(size_in_Rstar=0, min_color=None, max_color=None, save_plot=False, name="ionization_map", load_directly=False, parameter_list=None):
    """
    Generates a logarithmic ionization parameter map based on the stellar wind density, luminosity, and
    orbital parameters. The uncolored area in the map represents the X-ray shadow.

    Parameters
    ----------
    size_in_Rstar : float, optional
        Extent of the map from the stellar center in stellar radii. Default is 2 times the semimajor axis.
    min_color : float, optional
        Minimum value for the color scale of the ionization parameter. Default is None.
    max_color : float, optional
        Maximum value for the color scale of the ionization parameter. Default is None.
    save_plot : bool, optional
        If True, saves the generated plot. Default is False.
    name : str, optional
        Name of the file to save the plot. Default is "ionization_map".



    Notes
    -----
    A form will appear to input the necessary orbital parameters. These parameters are saved in a .txt file
    in the current directory and automatically loaded in subsequent runs. This avoids the need to re-enter parameters,
    allowing modification of only those that require adjustment.
    By setting load_directly=True the data will be authomatically loaded into the function without rising a form. By providing parameter_list=(list of parameters) the parameters will be trated as an imput and saved within the current directory for subsequent runs.
        


    Returns
    -------
    chi_result : pd.DataFrame
        DataFrame containing the ionization parameter map.
    area : float
        The calculated area between bounds in cm^2.
    """


    parameter_names = [
        "phase", "semimajor", "eccentricity", "periapsis", "Rstar", "Mstar1", "Mstar2",
        "wind_infinite_velocity", "Mass_loss_rate", "beta", "luminosity", "bound1", "bound2"
    ]
    
    # Load fixed values
    fixed_values = _manage_parameters(parameter_names, "ionization_map_phase",load_directly=load_directly,parameter_list=parameter_list )
    phase, semimajor, eccentricity, periapsis, Rstar, Mstar1, Mstar2, wind_infinite_velocity, Mass_loss_rate, beta, luminosity, bound1, bound2 = fixed_values
    
    # Calculate various parameters
    th = np.arange(0, 1, 0.001)
    abar = semimajor * max(Mstar1, Mstar2) / (Mstar1 + Mstar2)
    M_dot_grams = Mass_loss_rate * msun / (365 * 24 * 60 * 60)
    vinf_cm_s = wind_infinite_velocity * 100000
    Rstar_cm = Rstar * rsun_cm
    
    # Orbital radius calculations
    Rorb = (abar * (1 - eccentricity**2) / (1 + eccentricity * np.cos((phase - periapsis / 360) * 2 * np.pi))) * Rstar_cm
    Rorb_plot = (abar * (1 - eccentricity**2) / (1 + eccentricity * np.cos((th - periapsis / 360) * 2 * np.pi))) * Rstar_cm
    
    if size_in_Rstar==0:
        size_in_Rstar = 2*max(Rorb_plot/Rstar_cm)
    
    
    x = np.arange(1, size_in_Rstar, 0.01) * Rstar_cm
    v = vinf_cm_s * (1 - Rstar_cm / x)**beta
    ro = M_dot_grams / (4 * np.pi * v * x**2)
    
    # Calculate electron number density
    na = 6.02214076e23 / 1.00797
    ne = ro * na
    
    #..............................................Angles for shadow

    phase_ns = phase * 2 * np.pi
    phase_ns_degrees = np.degrees(phase_ns)

    alpha = np.arcsin(1/(Rorb/ Rstar_cm))
    alpha_degrees = np.degrees(alpha)

    alpha2 = np.arcsin(1/size_in_Rstar)
    alpha_degrees2 = np.degrees(alpha2)

    gamma_ = 180-(alpha_degrees2 + alpha_degrees)
    rho_ = gamma_+phase_ns_degrees

    rho = rho_*2*np.pi/360
    gamma = gamma_*2*np.pi/360

    phase_touch = (180-90-alpha_degrees)/360

    #..............................................
    
    # Create a DataFrame to store chi results
    chi_result = pd.DataFrame(index=np.round(x / Rstar_cm, 3))
    cmap = plt.get_cmap('rainbow')

    # Compute chi values
    for i in range(len(th)):
        distance = np.sqrt(Rorb**2 + x**2 - 2 * Rorb * x * np.cos(2 * np.pi * (th[i] - phase)))
        chi = luminosity / (ne * abs(distance)**2)
        chi_result[str(round(th[i], 3))] = chi

    # Determine color scale limits
    if not max_color:
        max_color = np.percentile(np.concatenate(chi_result.values), 90)
        print("max color coefficient is", round(max_color,2))  # Corrected print statement
        
    if not min_color:
        min_color = np.percentile(np.concatenate(chi_result.values), 10)
        print("min color coefficient is", round(min_color,2))

    # Initialize plot
    fig, axs = plt.subplots(1, 1, figsize=(20, 10), subplot_kw={'projection': 'polar'})
    norm = Normalize(vmin=np.log(min_color), vmax=np.log(max_color))
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    
    # Calculate area between bounds and plot
    area_between_bounds = 0
    
    bound_limit_1=[]
    bound_limit_2=[]
    bound_limit_3=[]
    bound_limit_4=[]
    phase_limit=[]
    phase_limit2=[]
    #....................................................................................
    th_ = np.arange(phase_touch+phase ,phase + gamma_/360, 0.001)
    for i in range(len(th_)):

        x_ = (Rorb/Rstar_cm) * np.sin(alpha_degrees * 2 * np.pi / 360) / np.sin((180 - (th_[i] - phase) * 360 - alpha_degrees) * 2 * np.pi / 360)
            
        if (x_ >= 1):
                
            x = np.arange(x_, size_in_Rstar, 0.01)*Rstar_cm
            v = vinf_cm_s * (1 - Rstar_cm / x)**beta
            ro = M_dot_grams / (4 * np.pi * v * x**2)
            ne = ro * na
        
            distance = np.sqrt(Rorb**2 + x**2 - 2 * Rorb * x * np.cos(2 * np.pi * (th_[i] - phase)))
            chi = luminosity / (ne * abs(distance)**2)
            x_chi_select = x[(chi >= bound1) & (chi <= bound2)] / Rstar_cm
            colors = cmap(norm(np.log(chi)))
        

            axs.scatter(np.tile(th_[i] * 2 * np.pi, len(x)), x / Rstar_cm, c=colors, cmap='rainbow',alpha=0.3)

            if len(x_chi_select) > 1:
                
                bound_limit_1.append(min(x_chi_select))
                bound_limit_2.append(max(x_chi_select))
                phase_limit.append(th_[i])
                
                if max(np.diff(x_chi_select))>0.02:
                
                    idx3 = np.argmax(np.diff(x_chi_select))
                    idx4 = np.argmax(np.diff(x_chi_select))+1
                    
                    bound_limit_3.append(x_chi_select[idx3])
                    bound_limit_4.append(x_chi_select[idx4])
                    phase_limit2.append(th_[i])
                    
    #....................................................................................
    th_ = np.arange(phase - gamma_/360 ,-phase_touch + phase, 0.001)
    for i in range(len(th_)):
            
        x_ = (Rorb/Rstar_cm)*np.sin(alpha_degrees*2*np.pi/360)/(np.sin((180-(phase-th_[i])*360-alpha_degrees)*2*np.pi/360))
            
        if (x_ >= 1):
                
            x = np.arange(x_, size_in_Rstar, 0.01)*Rstar_cm
            v = vinf_cm_s * (1 - Rstar_cm / x)**beta
            ro = M_dot_grams / (4 * np.pi * v * x**2)
            ne = ro * na
        
            distance = np.sqrt(Rorb**2 + x**2 - 2 * Rorb * x * np.cos(2 * np.pi * (th_[i] - phase)))
            chi = luminosity / (ne * abs(distance)**2)
            x_chi_select = x[(chi >= bound1) & (chi <= bound2)] / Rstar_cm
            colors = cmap(norm(np.log(chi)))
            distance = luminosity/(ne*chi)**0.5

            axs.scatter(np.tile(th_[i] * 2 * np.pi, len(x)), x / Rstar_cm, c=colors, cmap='rainbow',alpha=0.1)
            
            
            if len(x_chi_select) > 1:
            
                bound_limit_1.append(min(x_chi_select))
                bound_limit_2.append(max(x_chi_select))
                phase_limit.append(th_[i])
                
                if max(np.diff(x_chi_select))>0.02:
                
                    idx3 = np.argmax(np.diff(x_chi_select))
                    idx4 = np.argmax(np.diff(x_chi_select))+1
                    
                    bound_limit_3.append(x_chi_select[idx3])
                    bound_limit_4.append(x_chi_select[idx4])
                    phase_limit2.append(th_[i])

                    
    #....................................................................................
    th_ = np.arange(-phase_touch + phase ,phase_touch+phase, 0.001)
    for i in range(len(th_)):

        x = np.arange(1, size_in_Rstar,0.01)*Rstar_cm
        v = vinf_cm_s * (1 - Rstar_cm / x)**beta
        ro = M_dot_grams / (4 * np.pi * v * x**2)
        ne = ro * na
        
        distance = np.sqrt(Rorb**2 + x**2 - 2 * Rorb * x * np.cos(2 * np.pi * (th_[i] - phase)))
        chi = luminosity / (ne * abs(distance)**2)
        x_chi_select = x[(chi >= bound1) & (chi <= bound2)] / Rstar_cm
        colors = cmap(norm(np.log(chi)))
        

        axs.scatter(np.tile(th_[i] * 2 * np.pi, len(x)), x / Rstar_cm, c=colors, cmap='rainbow',alpha=0.1)

        if len(x_chi_select) > 1:
        
            bound_limit_1.append(min(x_chi_select))
            bound_limit_2.append(max(x_chi_select))
            phase_limit.append(th_[i])

            
            if max(np.diff(x_chi_select))>0.02:
                
                idx3 = np.argmax(np.diff(x_chi_select))
                idx4 = np.argmax(np.diff(x_chi_select))+1
                    
                bound_limit_3.append(x_chi_select[idx3])
                bound_limit_4.append(x_chi_select[idx4])
                phase_limit2.append(th_[i])

    
    # Plot additional elements
    axs.plot(np.linspace(0, 2 * np.pi, 100), np.ones(100), color='black')
    axs.plot(th * 2 * np.pi, Rorb_plot / Rstar_cm, color='black', linestyle='--',alpha=0.1)
    axs.plot(phase * 2 * np.pi, Rorb / Rstar_cm, color='black', marker='.')
    
    #Plot  and area between bounds
    ph_lim=np.array(phase_limit)
    ph_lim2=np.array(phase_limit2)
    x_bound1 = np.array(bound_limit_1)
    x_bound2 = np.array(bound_limit_2)
    x_bound3 = np.array(bound_limit_3)
    x_bound4 = np.array(bound_limit_4)
    
    sorted_indices = np.argsort(ph_lim)

    ph_lim = ph_lim[sorted_indices]*2*np.pi
    x_bound1 = x_bound1[sorted_indices]
    x_bound2 = x_bound2[sorted_indices]
    
    sorted_indices2 = np.argsort(ph_lim2)
    
    ph_lim2 = ph_lim2[sorted_indices2]*2*np.pi
    x_bound3 = x_bound3[sorted_indices2]
    x_bound4 = x_bound4[sorted_indices2]

    axs.plot(ph_lim ,x_bound1 ,"ko", markersize=1)
    axs.plot(ph_lim ,x_bound2,"ko", markersize=1)
    
    axs.plot(ph_lim2 ,x_bound3 ,"ko", markersize=1)
    axs.plot(ph_lim2 ,x_bound4,"ko", markersize=1)
    
    dph_lim = np.diff(ph_lim)
    dph_lim2 = np.diff(ph_lim2)

    area1 = 0.5 * np.sum((x_bound1[:-1] * x_bound1[1:]) * np.sin(dph_lim))
    area2 = 0.5 * np.sum((x_bound2[:-1] * x_bound2[1:]) * np.sin(dph_lim))
    
    area3 = 0.5 * np.sum((x_bound3[:-1] * x_bound3[1:]) * np.sin(dph_lim2))
    area4 = 0.5 * np.sum((x_bound4[:-1] * x_bound4[1:]) * np.sin(dph_lim2))
    
    area_sec_1 = np.abs(area2-area1)*Rstar_cm**2
    area_sec_2 = np.abs(area3-area4)*Rstar_cm**2
    
    area_between_bounds = np.abs(area_sec_1-area_sec_2)
       #.......................
    axs.set_theta_direction(-1)
    axs.set_theta_offset(np.pi / 2)
    
    cbar = plt.colorbar(sm, ax=axs, orientation='vertical')
    cbar.set_label(r'Log $\chi$')

    
    if save_plot:
        plt.savefig(f"{name}.png")
    
    return chi_result, area_between_bounds
