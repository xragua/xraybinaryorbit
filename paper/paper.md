---
title: 'xraybinaryorbit: A Python Package for Analyzing Orbital Modulations in X-ray Binaries'
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
    orcid: 0009-0007-9988-0202
    equal-contrib: true 
    affiliation: 1
  - name: Jose Miguel Torrejón
    orcid: 0000-0002-5967-5163
    equal-contrib: true
    affiliation: 1
  - name: Jose Joaquín Rodes-Roca
    orcid: 0000-0003-4363-8138
    equal-contrib: true
    affiliation: 1

affiliations:
 - name: Instituto Universitario de Física Aplicada a las Ciencias y las Tecnologías, Universidad de Alicante, 03690 Alicante, Spain
   index: 1

date: 11 November 2024
bibliography: paper.bib

---


# Summary

X-ray astronomy is a young discipline, spanning no more than a few decades. The Earth's atmosphere is opaque to this type of radiation, so observations in this part of the spectrum were delayed until the space era began, when rocket launchers carrying X-ray telescopes revealed the universe from a brand-new point of view [@1962PhRvL.9.439G].

X-ray binaries are tremendously interesting systems. In these extreme environments, a compact object (either a neutron star, a black hole or a white dwarf) draws in matter from a companion star, producing X-ray radiation. These systems offer a unique window into extreme physics, from the effects of strong gravity and relativistic jets to the presence of intense magnetic fields [@Carroll_Ostlie_2017;@1983bhwd.book.S;@lewin1997xray].

Orbital modulations are observed in nearly all X-ray binary systems. These variations arise from the orbital motions within the system, driven by the relative velocities of the two stars and/or their changing configurations with respect to each other and the observer. 

For instance, we might detect slight oscillating shifts in the center of emission lines and/or small variations in a neutron star's spin period. This phenomenon is likely caused by the Doppler effect. Our code translates this data into orbital parameters that explain these Doppler shifts. Although the process seems simple, it becomes complex due to factors like inclination, eccentricity, periapsis, and the masses of the stars involved. If the eccentricity is greater than zero, the velocity around the orbit is not constant, making the analysis even more intricate.

When considering stellar wind (the matter accreted by compact objects), various combinations of these factors result in orbital modulations. Eccentricity causes changes in wind density throughout the orbit, affecting the amount of accreted matter. Both eccentricity and inclination influence the absorption column faced by the emitted radiation and the ionization of the stellar wind, which also varies with orbital phase.

Although these orbital modulations are conceptually simple, they are challenging to analyze. At the same time, they can provide tremendous insights into the orbital mechanics and wind properties of our systems, may help to complete the puzzle when we already have some pieces of information or to give us a list of possible scenarios compatible with our data. This is where our tools provide valuable support.

Historically, the primary limitation in this type of analysis has been the lack of resolution for detailed phase-resolved observations. However, upcoming high-resolution missions like XRISM [@2022arXiv220205399X] and New Athena [@2016SPIE.9905E.2FB] promise to significantly improve the quality of these analyses. In addition to better resolution, advances in computational power have been crucial. Many of these tools have already been successfully applied to studies using XMM-Newton and Chandra data, enabling analyses that were previously impossible. [@2022MNRAS.512.304S;@2021MNRAS.501.5892S;@newcenx3]


# Statement of Need

Several software packages are available for the analysis of X-ray astronomy. Some notable Python-based options include Stingray [@2019ApJ88139H], which is dedicated to the timing analysis of astronomical data, and Jaxspec [@2024A&A.690A.317D], which specializes in spectral fitting using Bayesian inference. Lightkurve simplifies the analysis of time-series data from space missions like Kepler and TESS. A well-known and widely used comprehensive package for general astronomy computations and data handling is Astropy [@astropy:2022].

xraybinaryorbit is a dedicated software package specifically designed for studying orbital modulations in X-ray binary systems.
The tool features an intuitive interface for managing the parameters that affect orbital modulations. During the initial use, users input parameters via a form, and these are stored in a configuration file within the working directory. This file is automatically loaded in subsequent sessions, removing the need for repeated data entry (direct loads are also permited).

The software offers two fitting methods: least squares (LS) and particle swarm optimization (PSO) [@10.1162/EVCO_r_00180, @pyswarms]. The LS method is faster but may fail to converge in some cases, whereas PSO is more robust but computationally intensive. Some results obtained with this package are shown in Figures 1, 2 and 3.



# Background theory

This package primarily relies on the following theoretical frameworks:


### CAK Model
The CAK model [@1975ApJ195157C] describes radiation-driven winds in massive stars, with wind velocity, density, and ionization state varying with distance to the companion star.

$$ \rho = \frac{\dot{M}}{4 \pi v R^2} $$

Here, $\rho$  is wind density, $\dot{M}$ is the mass loss rate, v is the wind velocity and R is the distance to the companion. within this package we assume a spherical, smooth and un-ionized wind distribution.

### Accretion Luminosity and Ionization
Accretion powers many close binary systems [@Frank_King_Raine_2002]. The accretion luminosity is given by:

$$ L_{ac} = \frac{G M \dot{M}}{R} $$

where L$_{ac}$ is the luminosity, G is the gravitational constant, M is the compact object mass, and R is the compact object radius. The ionization parameter $\xi$ provides information about the ionization state of the stellar wind, which, in turn, affects its absorption power and the emission lines we might observe. This parameter can be approximated as:

$$ \xi = \frac{L_{\rm X}}{n(r_{\rm X}) r_{\rm X}^{2}} $$

### Doppler Effect

The Doppler effect refers to the change in frequency or wavelength of a wave as observed by someone who is moving relative to the wave source. The Doppler velocity induced by orbital motion and the observed Doppler-shifted wavelength can be written as:

$$ v_{D} = (-r\omega \sin\phi \sin i);  \lambda_{D} = \lambda_{\text{rest}}\left(1+\frac{v_{D}}{c}\right)$$


where r is the orbital radius, $\omega$ is angular velocity, $\i$ is inclination, and $\lambda_{\rm D}$ and $\lambda_{\rm rest}$ are the Doppler-shifted and rest wavelengths, respectively.

# Examples of Package Performance

![Reconstructed trajectory of Fe XXV emitting plasma in Cen X-3, derived from the observed Doppler shifts of the emission line center.](trajectory_and_shifts.png){#sylt width="100%"}

![Ionization map of the stellar wind in the orbital plane of 4U 0114+65, showing variations in different observed luminosities.](new_iopn_map.png){#sylt width="100%"}

![Theoretical local absorption column for a hypothetical X-Ray binary system, vs both time and orbital phase, calculated through the orbit.](NH_through_the_orbit.png){#sylt width="100%"}


# Acknowledgements

This research has been funded by the ASFAE/2022/02 project from the Generalitat Valenciana. 


# References
