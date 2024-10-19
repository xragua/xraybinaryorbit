# xraybinaryorbit


**X-ray binaries are truly fascinating!** In these extreme environments, a compact object—either a neutron star or black hole—draws in matter from a companion star, producing intense X-ray emissions. These systems offer a unique window into extreme physics, from the effects of strong gravity and relativistic jets to the presence of intense magnetic fields.

Orbital modulations are observed in nearly all X-ray binary systems. These variations arise from the orbital motions of the system, driven by the relative velocities of the two stars and their changing configurations with respect to each other and the observer.

To aid in the study of these modulations, we introduce **xraybinaryorbit** —a user-friendly Python package designed to simplify the analysis of orbital modulations in X-ray binaries. Whether you are studying the orbital parameters or working to refine your data analysis, this tool is built to help you extract the most from your observations.

If you have any questions or need assistance, please feel free to reach out: graciela.sanjurjo@ua.es.



## Getting started

### Installation
---------------------------------------------------------


You can install the package directly from PyPI using pip: **pip install xraybinaryorbit**.

Or download the code from [here](https://github.com/xragua/xraybinaryorbit/releases/tag/0.2.9).

Some examples of their usage are presented [here](https://github.com/xragua/xraybinaryorbit/tree/main/examples).


### Which orbital modulations are we talking about?
---------------------------------------------------------

We can observe how the center of an emission line slightly changes during a phase resolved analysis, or how the NS spin period slightly varies following a trend. These phenomena can be caused by Doppler effect. Our code will help you turn these observations into the orbital parameters that cause those Doppler shifts.

But this simple idea can get tricky when you consider all the factors involved. Inclination, eccentricity, periapsis, distance to the barycenter (which depends on the masses of the stars involved), and orbital phases really matter in this analysis. Plus, if there’s eccentricity different than 0, the velocity around the orbit isn’t constant. So, it’s not as straightforward as it might seem, right?.

If we think about stellar wind (the matter accreted by our compact objects) there are many combinations that can result into orbital modulations. With the eccentricity the density around the orbit changes, and thus, the accreted matter. With eccentricity and inclination the absorption column faced by the emitted radiation varies depending on the orbital phase, so does the ionization of the wind.

These orbital modulations are easy to grasp but not so straightforward to analyze—yet that’s where our tools come in.

The primary challenge in this type of analysis has long been the lack of sufficient resolution for detailed phase-resolved observations. However, upcoming high-resolution missions, such as **XRISM** and **New Athena**, promise to take these analyses to the next level. But it’s not just about improved resolution—advances in computational power are equally crucial. Many of these tools have already been successfully applied in studies using **XMM-Newton** and **Chandra** data, enabling analyses that would have been impossible just a few years ago.


So, dive in out **NOTEBOOKS** where we show some interesting examples, theoretial and real data, and start exploring the fascinating world of X-ray binaries with our tools!



