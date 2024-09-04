---
title: 'Orbital Dynamics in X-RAY stellar binary systems'
tags:
  - Python
  - x-ray astronomy
  - dynamics
  - doppler effectt
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

X-ray astronomy is a young discipline, spanning no more than a few decades. The Earth's atmosphere is opaque to this type of radiation, so observations in this part of the spectrum had to wait until the beginning of the space era, with rocket launchers carrying X-ray telescopes to reveal the universe from a brand-new point of view.

X-ray binary systems consist of two stars in close orbit around each other, where one of the stars is typically a compact object such as a neutron star or a black hole. The compact object accretes matter from its companion star, which can be a main sequence star, a giant star, or even another compact object. The X-ray radiation in these systems is generated through the accretion of matter from the companion's powerful stellar wind, typical of these early-type stars. Close binaries may become compact-object mergers and eventually sources of gravitational waves and/or short gamma-ray bursts. They will also provide insight into the behavior of matter under extreme gravitational and magnetic fields. Understanding these processes is fundamental to modern astrophysics and has driven numerous theoretical and observational studies.

The xraybinaryoorbit package helps to unveil X-Ray binary orbital dynamics based on key theories including the conservation of angular momentum in orbital mechanics, the CAK model for radiation-driven stellar winds, accretion luminosity, ionization parameters, and the Doppler effect. 


# Science behind

The functions contained in this package rely in the following theories:

## Conservation of angular momentum in orbital mechanics:

If the eccentricity of our system is different than 0, the orbital phase will not vary linearly with the observational time, as the speed will increase at periastron primarily due to the conservation of angular momentum, which dictates that as the compact object moves closer to the central star, it must travel faster to maintain the total angular momentum of the system. This relationship is further influenced by Kepler’s laws of planetary motion, which describe how objects sweep out equal areas in equal times (see [@2006ima..book.....C] as an example).

$$ r^2 \cdot \omega = h $$

We will take this fact into consideration in all our functions and provide dedicated functions to transform phase into time and vice versa.

## CAK model:

The CAK model, proposed by Castor, Abbott, and Klein in 1975 [@1975ApJ195157C], is a theoretical framework used to describe radiation-driven winds in massive stars. These stars have strong stellar winds driven by the interaction between their radiation and the surrounding material.

The CAK model provides a quantitative description of how the wind velocity, density, and ionization state vary with distance from the companion.

$$ \rho = \frac{\dot{M}}{4 \pi v R^2} $$

where $\rho$ is the density of the wind at a given distance R $\dot{M}$ is the mass accretion rate in units of mass per unit of time, v is the orbital speed at distances greater than the stellar radius and R is the distance to the star.
In this package, we assume that the wind is spherically distributed and unionized.

## Accretion Luminosity and Ionization Parameter:

Accretion is the process by which gravitational potential energy is extracted from material accreting onto a gravitating body  (see [@Frank_King_Raine_2002]). This phenomenon serves as the primary power source in various types of close binary systems and is also believed to fuel active galactic nuclei and quasars. When considering a flux of matter with an accretion rate $\dot{M}$, the resulting luminosity (assuming all mechanical energy is radiated) is defined as the accretion luminosity:

$$ L_{ac} = \frac{GM \dot{M}}{R} $$

where $L_{ac}$ is the accretion luminosity, G is the gravitational constant, M is the mass of the gravitating body, $\dot{M}$ is the accretion rate, and R is the characteristic radius associated with the accretion process.

The ionization parameter $\xi$ is defined as:

$$ \xi = \frac{L_{\rm X}}{n(r_{\rm X}) r_{\rm X}^{2}} $$

where $L_{\rm X}$ is the X-ray luminosity, $n(r_{\rm X})$ is the local particle density at a distance  $ r_{\rm X}$ from the X-ray source (such as a neutron star), $r_{\rm X}$ is the distance from the X-ray source.

This parameter quantifies the ionization state of the surrounding medium due to X-ray radiation from the neutron star. We provide a function which calculates the ionization map if the binary system plane taking into account these calculations within the CAK frame.

## Doppler Effect:

The Doppler effect, named after the Austrian physicist Christian Doppler who first proposed it in 1842, is the change in frequency or wavelength of a wave in relation to an observer moving relative to the source of the wave.

In astronomy, the Doppler effect is used to analyze the motion of celestial objects by observing shifts in their emitted light. By measuring Doppler shifts in the spectra of stars and galaxies, astronomers determine radial velocities, study galactic rotation, identify exoplanets, and explore the expansion of the universe through cosmological redshift. The Doppler effect plays a pivotal role in deciphering cosmic motions and unraveling the mysteries of the cosmos.

In the context of X-ray binaries, the Doppler effect is evident in the pulsations of a neutron star (NS) orbiting its companion, allowing precise determination of orbital parameters like radius, mass, inclination, and eccentricity. Additionally, the Doppler effect influences emission line energies when the emitting plasma is in motion.

The general equation for the Doppler velocity in terms of the orbital phase is:

$$ v_{D} = (-r\omega \sin\phi \sin i) $$

$$ \lambda_{D} = \lambda_{\text{rest}}\left(1+\frac{v_{D}}{c}\right) $$

where r is the orbital radius, a is the semimajor axis, b is the distance to the barycenter (the semimajor axis corrected by the reduced mass of the stellar system), e is the eccentricity, \phi is the orbital phase,  W is the angle to the periapsis, $\omega$ is the angular velocity, i is the inclination, and $\lambda_{\rm D}$ and $\lambda_{\rm rest}$ are the center of the emission line, Doppler shifted and at rest, respectively, in wavelength units.


Within the Fitting functions, we use a particle swarm approach ([@pyswarms], [@10.1162/EVCO_r_00180]) as a classical least squares algorithm does not always converge.

![Some results obtained with the functions contained in this package.](joss.jpg){#sylt width="100%"}

# Statement of Need

The study of orbital modulations in X-ray binaries is crucial for understanding their physical properties and dynamics. Upcoming telescopes, such as Athena's X-IFU [@2016SPIE9905E2FB] and XRISM [@2022arXiv220205399X], with their significantly higher resolution, are expected to greatly enhance these analyses, providing deeper insights into the intricate dynamics of X-ray systems.

However, it's not just the next-generation instruments that will advance our understanding. New analytical techniques and improved computing capacities will allow us to revisit and reanalyze existing data from current telescopes from a different point of view. This approach has already been successfully applied to data from Chandra and XMM, as demonstrated in recent studies [@2022MNRAS512304S; @2021MNRAS.501.5892S], and is also the focus of ongoing research currently under review.



# Acknowledgements

This research has been funded by the ASFAE/2022/02 project from the Generalitat Valenciana. 

# References

