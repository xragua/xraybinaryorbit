import numpy as np
import pandas as pd
from scipy.integrate import quad, cumulative_trapezoid
from pyswarm import pso
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import warnings
import tkinter as tk
from tkinter import messagebox
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks, peak_widths, peak_prominences, find_peaks_cwt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import inspect
import math
import matplotlib.gridspec as gridspec
from scipy.interpolate import CubicSpline
from ..helpers.data_helpers import (
    _manage_parameters,
    _define_x_y_sy,
    _copy_fields,
    _load_values_to_interface,
    _load_bounds_to_interface,
    _manage_bounds,
)
from ..helpers.math_helpers import (
    _gaussian,
    _time_pairs,
    _interpolate_pchip,
    _chi_squared_weighted,
    _chi_squared,
    _orbital_phase_to_time,
    _orbital_time_to_phase,
)

c = 299792458

msun = (1.98847 * 10**30) * 1000  # g
rsun_m = 696340 * 1000
rsun_cm = 696340 * 1000 * 100  # cm

kev_ams = 1.23984193

na = 6.02214076 * 10**23 / 1.00797
mu = 0.5
mp = 1.67e-24  # g

# Default hydrogen mass fraction used to convert mass density into n_H and N_H.
# A value near 0.70 is appropriate for a roughly solar-composition OB-star wind.
# Change this constant for systems with a measured non-standard H/He composition.
HYDROGEN_MASS_FRACTION = 0.70


def _relative_orbital_separation(phase, semimajor, eccentricity, periapsis, Rstar_cm):
    """Relative separation between the two stars, measured from the donor centre."""
    return (
        semimajor
        * (1 - eccentricity**2)
        / (
            1
            + eccentricity
            * np.cos((np.asarray(phase) - periapsis / 360) * 2 * np.pi)
        )
        * Rstar_cm
    )


def _wind_mass_density(radius, Rstar_cm, vinf_cm_s, M_dot_grams, beta):
    """CAK-like wind mass density in g cm^-3 for radii outside the donor."""
    radius = np.asarray(radius, dtype=float)
    rho = np.full(radius.shape, np.nan, dtype=float)
    valid = radius > Rstar_cm

    if np.any(valid):
        velocity = vinf_cm_s * (1 - Rstar_cm / radius[valid]) ** beta
        rho[valid] = M_dot_grams / (4 * np.pi * velocity * radius[valid] ** 2)

    return rho


def _validate_hydrogen_mass_fraction(hydrogen_mass_fraction):
    """Validate the hydrogen mass fraction used to convert rho into n_H."""
    if not 0 < hydrogen_mass_fraction <= 1:
        raise ValueError("hydrogen_mass_fraction must satisfy 0 < X_H <= 1.")


def _segment_intersects_star(Rorb, cos_alpha, z_start, z_end, Rstar_cm):
    """Return True when a line-of-sight segment crosses the stellar photosphere."""
    z_closest = np.clip(Rorb * cos_alpha, z_start, z_end)
    minimum_radius = np.sqrt(
        Rorb**2
        + z_closest**2
        - 2 * Rorb * z_closest * cos_alpha
    )
    return minimum_radius <= Rstar_cm


# DENSITY IN THE ORBIT ###########################################################################
def density_through_orbit_theoretical(
    resolution=0.01,
    show_plot=False,
    load_directly=False,
    parameter_list=None,
):
    """
    Calculates the wind mass density encountered by a compact object along
    its orbit, assuming a spherically symmetric CAK-like stellar wind.

    ``semimajor`` is the relative orbital semimajor axis in units of the
    donor radius. No barycentric correction is applied.

    Returns
    -------
    time : array-like
        Time corresponding to the orbital positions, in seconds.
    phase : array-like
        Orbital phase.
    density : array-like
        Wind mass density in g cm^-3.
    """

    parameter_names = [
        "semimajor",
        "orbitalperiod",
        "eccentricity",
        "periapsis",
        "Rstar",
        "wind_infinite_velocity",
        "Mass_loss_rate",
        "beta",
    ]

    fixed_values = _manage_parameters(
        parameter_names,
        "density_through_orbit",
        load_directly=load_directly,
        parameter_list=parameter_list,
    )

    (
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        Rstar,
        wind_infinite_velocity,
        Mass_loss_rate,
        beta,
    ) = fixed_values

    th = np.arange(0, 1, resolution)

    # The unused spatial arguments remain in the helper signature only for
    # backward compatibility.
    _, time, _ = _orbital_phase_to_time(
        th,
        0,
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        Rstar,
        1.0,
        1.0,
        precision=0.01,
    )

    M_dot_grams = Mass_loss_rate * msun / (365 * 24 * 60 * 60)
    vinf_cm_s = wind_infinite_velocity * 100000
    Rstar_cm = Rstar * rsun_cm

    Rorb = _relative_orbital_separation(
        th,
        semimajor,
        eccentricity,
        periapsis,
        Rstar_cm,
    )

    rho_all = _wind_mass_density(
        Rorb,
        Rstar_cm,
        vinf_cm_s,
        M_dot_grams,
        beta,
    )
    valid = np.isfinite(rho_all)
    rho = rho_all[valid]

    if show_plot:
        fig = plt.figure(figsize=(9, 3))
        gs = gridspec.GridSpec(1, 3, width_ratios=[1.2, 1.2, 1.2])

        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(time[valid], rho)
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel(r"Wind mass density (g cm$^{-3}$)")

        ax2 = fig.add_subplot(gs[0, 1])
        ax2.plot(th[valid], rho)
        ax2.set_xlabel("Orbital phase")
        ax2.set_ylabel(r"Wind mass density (g cm$^{-3}$)")

        ax3 = fig.add_subplot(gs[0, 2], projection="polar")
        ax3.plot(th * 2 * np.pi, Rorb / Rstar_cm)
        ax3.set_theta_zero_location("N")
        ax3.set_theta_direction(-1)
        ax3.set_frame_on(False)

        plt.tight_layout()
        plt.savefig("density_through_the_orbit.png", bbox_inches="tight")

    return time[valid], th[valid], rho


# ABSORPTION COLUMN #############################################################################
def absorption_column_through_orbit_theoretical(
    resolution=0.01,
    show_plot=True,
    load_directly=False,
    parameter_list=None,
):
    """
    Calculates the equivalent hydrogen column density along the line of sight.

    ``HYDROGEN_MASS_FRACTION`` is the hydrogen mass fraction X_H used in

        n_H = X_H rho / m_p.

    A value close to 0.70 is a suitable generic default for a nearly solar
    H/He mixture, but it may be changed for a helium-enriched donor.

    Returns
    -------
    time : array-like
        Time array in seconds.
    phase : array-like
        Orbital phase array.
    NH : array-like
        Equivalent hydrogen column density in units of 10^22 cm^-2.

    Notes
    -----
    If the line of sight crosses the stellar photosphere, the direct source is
    occulted and NH is not directly observable. For a continuous orbital model,
    eclipse samples are assigned the maximum finite NH calculated outside the
    eclipse over the sampled orbit.
    """

    parameter_names = [
        "semimajor",
        "orbitalperiod",
        "eccentricity",
        "periapsis",
        "inclination",
        "Rstar",
        "wind_infinite_velocity",
        "Mass_loss_rate",
        "beta",
    ]

    fixed_values = _manage_parameters(
        parameter_names,
        "absorption_column_through_orbit",
        load_directly=load_directly,
        parameter_list=parameter_list,
    )

    (
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        inclination,
        Rstar,
        wind_infinite_velocity,
        Mass_loss_rate,
        beta,
    ) = fixed_values

    hydrogen_mass_fraction = HYDROGEN_MASS_FRACTION
    _validate_hydrogen_mass_fraction(hydrogen_mass_fraction)

    th = np.arange(-0.1, 1.1, resolution)

    _, time, _ = _orbital_phase_to_time(
        th,
        0,
        semimajor,
        orbitalperiod,
        eccentricity,
        periapsis,
        Rstar,
        1.0,
        1.0,
        precision=0.01,
    )

    M_dot_grams = Mass_loss_rate * msun / (365 * 24 * 60 * 60)
    vinf_cm_s = wind_infinite_velocity * 100000
    Rstar_cm = Rstar * rsun_cm

    Rorb_plot = _relative_orbital_separation(
        th,
        semimajor,
        eccentricity,
        periapsis,
        Rstar_cm,
    )

    z_upper = Rstar_cm * 1000
    nh = []
    eclipse_mask = []

    for phase, Rorb in zip(th, Rorb_plot):
        cos_alpha = np.clip(
            np.cos(phase * 2 * np.pi)
            * np.sin(np.deg2rad(inclination)),
            -1.0,
            1.0,
        )

        occulted = Rorb <= Rstar_cm or _segment_intersects_star(
            Rorb,
            cos_alpha,
            0.0,
            z_upper,
            Rstar_cm,
        )
        eclipse_mask.append(occulted)

        if occulted:
            # The direct source is not observable in eclipse. Store NaN for now
            # and replace it below with the maximum visible orbital NH.
            nh.append(np.nan)
            continue

        def integrand(z):
            radius = np.sqrt(
                Rorb**2
                + z**2
                - 2 * Rorb * z * cos_alpha
            )
            velocity = vinf_cm_s * (1 - Rstar_cm / radius) ** beta
            return M_dot_grams / (4 * np.pi * velocity * radius**2)

        mass_column, _ = quad(
            integrand,
            0.0,
            z_upper,
            limit=200,
        )

        nh.append(
            hydrogen_mass_fraction
            * mass_column
            / mp
            / 1e22
        )

    nh = np.asarray(nh, dtype=float)
    eclipse_mask = np.asarray(eclipse_mask, dtype=bool)

    if np.any(eclipse_mask):
        visible_mask = (~eclipse_mask) & np.isfinite(nh)
        if not np.any(visible_mask):
            raise ValueError(
                "The compact object is occulted at every sampled orbital phase; "
                "a finite out-of-eclipse NH maximum cannot be defined."
            )
        nh[eclipse_mask] = np.nanmax(nh[visible_mask])

    if show_plot:
        fig = plt.figure(figsize=(12, 3))
        gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1.2])

        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(time, nh)
        ax1.set_xlabel("Time (s)")
        ax1.set_ylabel(r"$N_{\rm H}$ ($10^{22}$ cm$^{-2}$)")

        ax2 = fig.add_subplot(gs[0, 1])
        ax2.plot(th, nh)
        ax2.set_xlabel("Orbital phase")
        ax2.set_ylabel(r"$N_{\rm H}$ ($10^{22}$ cm$^{-2}$)")

        ax3 = fig.add_subplot(gs[0, 2], projection="polar")
        ax3.plot(th * 2 * np.pi, Rorb_plot / Rstar_cm)
        ax3.set_theta_zero_location("N")
        ax3.set_theta_direction(-1)
        ax3.set_frame_on(False)

        plt.tight_layout()
        plt.savefig("NH_through_the_orbit.png")

    return time, th, nh


# DENSITY AND LOG XI ###########################################################################
def density_and_ionization_orbital_phase_theoretical(
    resolution=0.01,
    size=10,
    show_plot=True,
    load_directly=False,
    parameter_list=None,
):
    """
    Calculates n_H, log10(xi), and the remaining hydrogen column along a
    ray from the compact object towards the observer.

    The ionization parameter is defined as

        xi = L / (n_H d^2),

    where d is the distance from the compact object. The input ``luminosity``
    retains the original convention of units of 10^32 erg s^-1.

    Returns
    -------
    z : array-like
        Distance from the compact object in units of the donor radius.
    density : array-like
        Hydrogen number density n_H in cm^-3.
    chi : array-like
        log10(xi), with xi in erg cm s^-1.
    NH : array-like
        Remaining equivalent hydrogen column in units of 10^22 cm^-2.
    """

    parameter_names = [
        "orb_phase",
        "luminosity",
        "semimajor",
        "eccentricity",
        "periapsis",
        "inclination",
        "Rstar",
        "wind_infinite_velocity",
        "Mass_loss_rate",
        "beta",
    ]

    fixed_values = _manage_parameters(
        parameter_names,
        "den_chi_orbphase",
        load_directly=load_directly,
        parameter_list=parameter_list,
    )

    (
        orb_phase,
        luminosity,
        semimajor,
        eccentricity,
        periapsis,
        inclination,
        Rstar,
        wind_infinite_velocity,
        Mass_loss_rate,
        beta,
    ) = fixed_values

    hydrogen_mass_fraction = HYDROGEN_MASS_FRACTION
    _validate_hydrogen_mass_fraction(hydrogen_mass_fraction)

    if size <= 0:
        raise ValueError("size must be greater than zero.")

    th = np.arange(0, 1, resolution)

    luminosity_ = luminosity * 1e32
    M_dot_grams = Mass_loss_rate * msun / (365 * 24 * 60 * 60)
    vinf_cm_s = wind_infinite_velocity * 100000
    Rstar_cm = Rstar * rsun_cm

    Rorb_plot = _relative_orbital_separation(
        th,
        semimajor,
        eccentricity,
        periapsis,
        Rstar_cm,
    )
    Rorb = float(
        _relative_orbital_separation(
            orb_phase,
            semimajor,
            eccentricity,
            periapsis,
            Rstar_cm,
        )
    )

    z = np.linspace(0.0, Rorb * size, 10000, endpoint=False)

    cos_alpha = np.clip(
        np.cos(orb_phase * 2 * np.pi)
        * np.sin(np.deg2rad(inclination)),
        -1.0,
        1.0,
    )

    radius_from_donor = np.sqrt(
        Rorb**2
        + z**2
        - 2 * Rorb * z * cos_alpha
    )

    rho = _wind_mass_density(
        radius_from_donor,
        Rstar_cm,
        vinf_cm_s,
        M_dot_grams,
        beta,
    )

    # n_H, not total particle density.
    density = hydrogen_mass_fraction * rho / mp

    chi = np.full(z.shape, np.nan, dtype=float)
    valid_chi = (
        np.isfinite(density)
        & (density > 0)
        & (z > 0)
    )
    chi[valid_chi] = np.log10(
        luminosity_
        / (density[valid_chi] * z[valid_chi] ** 2)
    )

    # Integrate the mass column from every z point to the end of the ray.
    rho_for_integral = np.nan_to_num(rho, nan=0.0, posinf=0.0, neginf=0.0)
    mass_column_from_z = -cumulative_trapezoid(
        rho_for_integral[::-1],
        z[::-1],
        initial=0.0,
    )[::-1]

    nh_from_z = (
        hydrogen_mass_fraction
        * mass_column_from_z
        / mp
        / 1e22
    )

    # A point whose remaining line of sight intersects the star is occulted.
    z_closest = np.clip(Rorb * cos_alpha, z, z[-1])
    minimum_radius = np.sqrt(
        Rorb**2
        + z_closest**2
        - 2 * Rorb * z_closest * cos_alpha
    )
    occulted = minimum_radius <= Rstar_cm
    nh_from_z[occulted] = 0.0

    if show_plot:
        fig = plt.figure(figsize=(20, 5))
        gs = fig.add_gridspec(1, 4)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[0, 2])
        ax4 = fig.add_subplot(gs[0, 3], projection="polar")

        z_rstar = z / Rstar_cm
        positive_z = z_rstar > 0

        ax1.plot(z_rstar[positive_z], density[positive_z])
        ax1.set_xlabel(r"Path from compact object ($R_\star$)")
        ax1.set_ylabel(r"$n_{\rm H}$ (cm$^{-3}$)")
        ax1.set_xscale("log")
        ax1.set_yscale("log")

        ax2.plot(z_rstar[positive_z], chi[positive_z])
        ax2.set_xlabel(r"Path from compact object ($R_\star$)")
        ax2.set_ylabel(r"$\log_{10}\xi$")
        ax2.set_xscale("log")

        ax3.plot(z_rstar[positive_z], nh_from_z[positive_z])
        ax3.set_xlabel(r"Path from compact object ($R_\star$)")
        ax3.set_ylabel(r"$N_{\rm H}$ ($10^{22}$ cm$^{-2}$)")
        ax3.set_xscale("log")

        ax4.plot(th * 2 * np.pi, Rorb_plot / Rstar_cm)
        ax4.set_theta_zero_location("N")
        ax4.set_theta_direction(-1)
        ax4.spines["polar"].set_visible(False)
        ax4.set_xticks([])

        plt.tight_layout()
        plt.savefig("density_and_chi_orbital_phase.png")

    return z / Rstar_cm, density, chi, np.asarray(nh_from_z)


# IONIZATION PARAMETER MAP ####################################################################
def ionization_map_phase(
    size_in_Rstar=0,
    min_color=None,
    max_color=None,
    save_plot=False,
    name="ionization_map",
    load_directly=False,
    parameter_list=None,
):
    """
    Generates an ionization-parameter map for a CAK-like stellar wind.

    ``semimajor`` is the relative semimajor axis in donor-radius units.
    ``luminosity`` is supplied in erg s^-1, preserving the original
    convention of this function. ``HYDROGEN_MASS_FRACTION`` is used to
    calculate n_H for
    xi = L / (n_H d^2). The second returned DataFrame remains the electron
    density for compatibility with the original output.

    Returns
    -------
    chi_result : pandas.DataFrame
        Ionization parameter xi in erg cm s^-1.
    ne_result : pandas.DataFrame
        Approximate electron density in cm^-3, assuming fully ionized H/He
        and neglecting the small metal contribution.
    area_between_bounds : float
        Area satisfying the selected xi bounds, in cm^2.
    """

    parameter_names = [
        "phase",
        "semimajor",
        "eccentricity",
        "periapsis",
        "Rstar",
        "wind_infinite_velocity",
        "Mass_loss_rate",
        "beta",
        "luminosity",
        "bound1",
        "bound2",
    ]

    fixed_values = _manage_parameters(
        parameter_names,
        "ionization_map_phase",
        load_directly=load_directly,
        parameter_list=parameter_list,
    )

    (
        phase,
        semimajor,
        eccentricity,
        periapsis,
        Rstar,
        wind_infinite_velocity,
        Mass_loss_rate,
        beta,
        luminosity,
        bound1,
        bound2,
    ) = fixed_values

    hydrogen_mass_fraction = HYDROGEN_MASS_FRACTION
    _validate_hydrogen_mass_fraction(hydrogen_mass_fraction)

    th = np.arange(0, 1, 0.001)
    M_dot_grams = Mass_loss_rate * msun / (365 * 24 * 60 * 60)
    vinf_cm_s = wind_infinite_velocity * 100000
    Rstar_cm = Rstar * rsun_cm

    Rorb = float(
        _relative_orbital_separation(
            phase,
            semimajor,
            eccentricity,
            periapsis,
            Rstar_cm,
        )
    )
    Rorb_plot = _relative_orbital_separation(
        th,
        semimajor,
        eccentricity,
        periapsis,
        Rstar_cm,
    )

    if Rorb <= Rstar_cm:
        raise ValueError("The compact object lies inside the donor radius.")

    if size_in_Rstar == 0:
        size_in_Rstar = 2 * np.max(Rorb_plot / Rstar_cm)

    if size_in_Rstar <= 1.001:
        raise ValueError("size_in_Rstar must be greater than 1.")

    radial_grid = np.arange(1.001, size_in_Rstar, 0.01) * Rstar_cm
    rho_grid = _wind_mass_density(
        radial_grid,
        Rstar_cm,
        vinf_cm_s,
        M_dot_grams,
        beta,
    )

    n_h_grid = hydrogen_mass_fraction * rho_grid / mp

    # For a fully ionized H/He mixture with Z neglected:
    # n_e = (X_H + Y/2) rho/m_p = (1 + X_H) rho/(2 m_p).
    n_e_grid = (
        (1 + hydrogen_mass_fraction)
        * rho_grid
        / (2 * mp)
    )

    phase_ns_degrees = np.degrees(phase * 2 * np.pi)
    alpha = np.arcsin(Rstar_cm / Rorb)
    alpha_degrees = np.degrees(alpha)

    alpha2 = np.arcsin(1 / size_in_Rstar)
    alpha_degrees2 = np.degrees(alpha2)

    gamma_degrees = 180 - (alpha_degrees2 + alpha_degrees)
    phase_touch = (90 - alpha_degrees) / 360

    chi_result = pd.DataFrame(
        index=np.round(radial_grid / Rstar_cm, 3)
    )
    ne_result = pd.DataFrame(
        index=np.round(radial_grid / Rstar_cm, 3)
    )
    cmap = plt.get_cmap("rainbow")

    for theta in th:
        distance = np.sqrt(
            Rorb**2
            + radial_grid**2
            - 2
            * Rorb
            * radial_grid
            * np.cos(2 * np.pi * (theta - phase))
        )
        chi = luminosity / (n_h_grid * distance**2)
        column_name = str(round(theta, 3))
        chi_result[column_name] = chi
        ne_result[column_name] = n_e_grid

    finite_chi = chi_result.to_numpy().ravel()
    finite_chi = finite_chi[np.isfinite(finite_chi) & (finite_chi > 0)]

    if finite_chi.size == 0:
        raise ValueError("No finite positive ionization-parameter values were produced.")

    if max_color is None:
        max_color = np.percentile(finite_chi, 90)
        print("max color coefficient is", round(max_color, 2))

    if min_color is None:
        min_color = np.percentile(finite_chi, 10)
        print("min color coefficient is", round(min_color, 2))

    if min_color <= 0 or max_color <= min_color:
        raise ValueError("Color limits must satisfy 0 < min_color < max_color.")

    fig, axs = plt.subplots(
        1,
        1,
        figsize=(20, 10),
        subplot_kw={"projection": "polar"},
    )
    norm = Normalize(
        vmin=np.log10(min_color),
        vmax=np.log10(max_color),
    )
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    bound_limit_1 = []
    bound_limit_2 = []
    bound_limit_3 = []
    bound_limit_4 = []
    phase_limit = []
    phase_limit2 = []

    def process_ray(theta, minimum_radius_rstar, alpha_value):
        if not np.isfinite(minimum_radius_rstar):
            return
        if minimum_radius_rstar < 1 or minimum_radius_rstar >= size_in_Rstar:
            return

        x = np.arange(minimum_radius_rstar, size_in_Rstar, 0.01) * Rstar_cm
        if x.size < 2:
            return

        rho = _wind_mass_density(
            x,
            Rstar_cm,
            vinf_cm_s,
            M_dot_grams,
            beta,
        )
        n_h = hydrogen_mass_fraction * rho / mp

        distance = np.sqrt(
            Rorb**2
            + x**2
            - 2
            * Rorb
            * x
            * np.cos(2 * np.pi * (theta - phase))
        )
        chi = luminosity / (n_h * distance**2)

        valid = np.isfinite(chi) & (chi > 0)
        if not np.any(valid):
            return

        colors = cmap(norm(np.log10(chi[valid])))
        axs.scatter(
            np.full(np.count_nonzero(valid), theta * 2 * np.pi),
            x[valid] / Rstar_cm,
            c=colors,
            alpha=alpha_value,
            s=5,
        )

        x_chi_select = x[
            valid & (chi >= bound1) & (chi <= bound2)
        ] / Rstar_cm

        if len(x_chi_select) > 1:
            bound_limit_1.append(np.min(x_chi_select))
            bound_limit_2.append(np.max(x_chi_select))
            phase_limit.append(theta)

            gaps = np.diff(x_chi_select)
            if gaps.size and np.max(gaps) > 0.02:
                gap_index = np.argmax(gaps)
                bound_limit_3.append(x_chi_select[gap_index])
                bound_limit_4.append(x_chi_select[gap_index + 1])
                phase_limit2.append(theta)

    # Rays outside one side of the geometrical shadow.
    theta_values = np.arange(
        phase_touch + phase,
        phase + gamma_degrees / 360,
        0.001,
    )
    for theta in theta_values:
        denominator = np.sin(
            np.deg2rad(
                180
                - (theta - phase) * 360
                - alpha_degrees
            )
        )
        if np.isclose(denominator, 0):
            continue
        x_start = (
            (Rorb / Rstar_cm)
            * np.sin(np.deg2rad(alpha_degrees))
            / denominator
        )
        process_ray(theta, x_start, 0.3)

    # Rays outside the other side of the geometrical shadow.
    theta_values = np.arange(
        phase - gamma_degrees / 360,
        -phase_touch + phase,
        0.001,
    )
    for theta in theta_values:
        denominator = np.sin(
            np.deg2rad(
                180
                - (phase - theta) * 360
                - alpha_degrees
            )
        )
        if np.isclose(denominator, 0):
            continue
        x_start = (
            (Rorb / Rstar_cm)
            * np.sin(np.deg2rad(alpha_degrees))
            / denominator
        )
        process_ray(theta, x_start, 0.1)

    # Rays in front of the donor, starting at the photosphere.
    theta_values = np.arange(
        -phase_touch + phase,
        phase_touch + phase,
        0.001,
    )
    for theta in theta_values:
        process_ray(theta, 1.001, 0.1)

    axs.plot(np.linspace(0, 2 * np.pi, 100), np.ones(100), color="black")
    axs.plot(
        th * 2 * np.pi,
        Rorb_plot / Rstar_cm,
        color="black",
        linestyle="--",
        alpha=0.1,
    )
    axs.plot(
        phase * 2 * np.pi,
        Rorb / Rstar_cm,
        color="black",
        marker=".",
    )

    ph_lim = np.asarray(phase_limit)
    ph_lim2 = np.asarray(phase_limit2)
    x_bound1 = np.asarray(bound_limit_1)
    x_bound2 = np.asarray(bound_limit_2)
    x_bound3 = np.asarray(bound_limit_3)
    x_bound4 = np.asarray(bound_limit_4)

    if ph_lim.size:
        sorted_indices = np.argsort(ph_lim)
        ph_lim = ph_lim[sorted_indices] * 2 * np.pi
        x_bound1 = x_bound1[sorted_indices]
        x_bound2 = x_bound2[sorted_indices]
        axs.plot(ph_lim, x_bound1, "ko", markersize=1)
        axs.plot(ph_lim, x_bound2, "ko", markersize=1)

    if ph_lim2.size:
        sorted_indices2 = np.argsort(ph_lim2)
        ph_lim2 = ph_lim2[sorted_indices2] * 2 * np.pi
        x_bound3 = x_bound3[sorted_indices2]
        x_bound4 = x_bound4[sorted_indices2]
        axs.plot(ph_lim2, x_bound3, "ko", markersize=1)
        axs.plot(ph_lim2, x_bound4, "ko", markersize=1)

    area1 = 0.0
    area2 = 0.0
    area3 = 0.0
    area4 = 0.0

    if ph_lim.size > 1:
        dph_lim = np.diff(ph_lim)
        area1 = 0.5 * np.sum(
            x_bound1[:-1] * x_bound1[1:] * np.sin(dph_lim)
        )
        area2 = 0.5 * np.sum(
            x_bound2[:-1] * x_bound2[1:] * np.sin(dph_lim)
        )

    if ph_lim2.size > 1:
        dph_lim2 = np.diff(ph_lim2)
        area3 = 0.5 * np.sum(
            x_bound3[:-1] * x_bound3[1:] * np.sin(dph_lim2)
        )
        area4 = 0.5 * np.sum(
            x_bound4[:-1] * x_bound4[1:] * np.sin(dph_lim2)
        )

    area_sec_1 = abs(area2 - area1) * Rstar_cm**2
    area_sec_2 = abs(area3 - area4) * Rstar_cm**2
    area_between_bounds = abs(area_sec_1 - area_sec_2)

    axs.set_theta_direction(-1)
    axs.set_theta_offset(np.pi / 2)

    cbar = plt.colorbar(sm, ax=axs, orientation="vertical")
    cbar.set_label(r"$\log_{10}\xi$")

    if save_plot:
        plt.savefig(f"{name}.png")

    return chi_result, ne_result, area_between_bounds


###################################### ORBITAL PHASE TO TIME ##############################
# Orbital phase to time approximation (constant areolar velocity)
# Orbital time to phase (constant areolar velocity and interpolation)
##########################################################################################
