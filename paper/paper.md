---
title: 'Orbital Dynamics in X-RAY stellar binary systems'
tags:
  - Python
  - x-ray astronomy
  - dynamics
  - doppler effect
  - Binary stars
  - Neutron stars
  - Black holes
authors:
  - name: Graciela Sanjurjo-Ferrín
    orcid: 0009-0001-0378-7879
    equal-contrib: true
    affiliation: 1
  - name: Jessica Planelles Villalva 
    equal-contrib: true 
    affiliation: 1
  - name: Jose Miguel Torrejón
    equal-contrib: true
    affiliation: 1
  - name: Jose Joaquín Rodes-Roca
    equal-contrib: true
    affiliation: 1

affiliations:
 - name: Instituto Universitario de Física Aplicada a las Ciencias y las Tecnologías, Universidad de Alicante, 03690 Alicante, Spain
   index: 1

date: 21 March 2024
bibliography: paper.bib

---


# Summary

X-ray astronomy is a young discipline, spanning no more than a few decades. The Earth's atmosphere is opaque to this type of radiation, so observations in this part of the spectrum were delayed until the space era began, when rocket launchers carrying X-ray telescopes revealed the universe from a brand-new point of view.

X-ray binaries are tremendously interesting systems. In these extreme environments, a compact object—either a neutron star or black hole—draws in matter from a companion star, producing intense X-ray emissions. These systems offer a unique window into extreme physics, from the effects of strong gravity and relativistic jets to the presence of intense magnetic fields.

Orbital modulations are observed in nearly all X-ray binary systems. These variations arise from the orbital motions of the system, driven by the relative velocities of the two stars and their changing configurations with respect to each other and the observer.

To aid in the study of these modulations, we introduce **xraybinaryorbit**, a user-friendly Python package designed to simplify the analysis of orbital modulations in X-ray binaries.
There are three groups of funtions within this package, theoretical, fitting and timing. The software provides a user-friendly interface for managing the various parameters that influence orbital modulations. Upon first use, the user inputs parameters through a form, which are then saved in a file within the running directory, this file is automatically loaded in future runs, eliminating the need to re-enter all parameters. The software offers two fitting approaches: least squares (LS) and particle swarm optimization (PSO) [@10.1162/EVCO_r_00180,@pyswarms]. The LS methods is faster but may fail to converge in certain cases, whereas the PS0 is more robust but computationally intensive.  Some results obtained with this packages are shown in Figure 1. 

There are many software that can be applied to X-ray binary system analysis, including **Stingray** [@2019ApJ88139H] which focus on timing analysis of astronomical data, **Jaxspec** which specializes in spectral fitting using Bayesian inference. **Lightkurve**,  facilitates the analysis of time-series data from space missions like Kepler and TESS. 
A well know comprehensive package for general astronomy computations and data handling is **Astropy** [@astropy:2022]. 
In contrast **xraybinaryorbit** is highly specialized, providing detailed insights into orbital modulations within X-ray binary systems. 

# Science Behind
### Conservation of Angular Momentum
In systems with non-zero eccentricity, the orbital phase varies non-linearly with observational time due to increased speed at periastron, governed by the conservation of angular momentum. Kepler’s laws further explain that objects sweep out equal areas in equal times [@2006imabookC].

$$ r^2 \cdot \omega = Area $$

We account for this in all functions and provide tools to convert phase to time and vice versa.

### CAK Model
The CAK model [@1975ApJ195157C] describes radiation-driven winds in massive stars, with wind velocity, density, and ionization state varying with distance to the companion star.

$$ \rho = \frac{\dot{M}}{4 \pi v R^2} $$

Here, $\rho$  is wind density, $\dot{M}$ is the mass loss rate, v is the wind velocity and R is the distance to the companion. We assume a spherical, smooth and un-ionized wind distribution.

### Accretion Luminosity and Ionization
Accretion powers many close binary systems [@Frank_King_Raine_2002]. The accretion luminosity is given by:

$$ L_{ac} = \frac{G M \dot{M}}{R} $$

where L$_{ac}$ is the luminosity, G is the gravitational constant, M is the compact object mass, and R is the compact object radius. The ionization parameter $\xi$ is defined as:

$$ \xi = \frac{L_{\rm X}}{n(r_{\rm X}) r_{\rm X}^{2}} $$

### Doppler Effect
The Doppler effect refers to the change in frequency or wavelength of a wave in relation to an observer moving relative to the source. 

The Doppler velocity provoked by orbital motion can be calculated as:

$$ v_{D} = (-r\omega \sin\phi \sin i) $$

The observed Doppler length wave can be written as:

$$ \lambda_{D} = \lambda_{\text{rest}}\left(1+\frac{v_{D}}{c}\right) $$

where r is the orbital radius, $\omega$ is angular velocity, $\i$ is inclination, and $\lambda_{\rm D}$ and $\lambda_{\rm rest}$ are the Doppler-shifted and rest wavelengths, respectively.


![Some results obtained with the functions contained in this package.](joss.jpg){#sylt width="100%"}


# Statement of Need

Within our observations, we might detect slight oscillating shifts in the center of emission lines and/or small variations in a neutron star's spin period. This phenomenon is likely caused by the Doppler effect. Our code helps translate this data into the orbital parameters responsible for such Doppler shifts.

However, this seemingly simple process becomes complex when accounting for the several factors involved: inclination, eccentricity, periapsis, and distance to the barycenter (which depend on the masses of the stars involved). If the eccentricity is greater than zero, the velocity around the orbit is not constant, making the analysis even more intricate.

When considering stellar wind (the matter accreted by compact objects), various combinations of these factors result in orbital modulations. Eccentricity causes changes in wind density throughout the orbit, affecting the amount of accreted matter. Both eccentricity and inclination influence the absorption column faced by the emitted radiation and the ionization of the stellar wind, which also varies with orbital phase.

Although these orbital modulations are conceptually simple, they are challenging to analyze. At the same time, they can provide tremendous insights into the orbital mechanics and wind properties of our systems, may help to complete the puzzle when we already have some pieces of information or to give us a list of possible scenarios compatible with our data. This is where our tools provide valuable support.

Historically, the primary limitation in this type of analysis has been the lack of resolution for detailed phase-resolved observations. However, upcoming high-resolution missions like XRISM [@2022arXiv220205399X] and New Athena [2016SPIE.9905E..2FB] promise to significantly improve the quality of these analyses. In addition to better resolution, advances in computational power have been crucial. Many of these tools have already been successfully applied to studies using XMM-Newton and Chandra data, enabling analyses that were previously impossible. [@2022MNRAS512304S;@2021MNRAS.501.5892S;@newcenx3]


# Acknowledgements

This research has been funded by the ASFAE/2022/02 project from the Generalitat Valenciana. 


# References
