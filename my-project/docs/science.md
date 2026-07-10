# Scientific background and numerical methods

`xraybinaryorbit` provides simplified, interpretable models for orbital modulation in X-ray binaries. This page summarizes the physical assumptions and numerical methods used by the package. Function signatures and return values are documented separately in [Functions](functions.md).

---

## Conventions and units

The current development version uses the following conventions.

| Quantity | Symbol | Package convention |
|---|---:|---|
| Time | $t$ | seconds |
| Orbital phase | $\phi$ | cycles; one complete orbit spans 1 |
| Initial phase | `iphase` | cycles |
| Orbital period | $P_{\rm orb}$ | days in the orbital models |
| Eccentricity | $e$ | $0\leq e<1$ |
| Argument of periapsis | $\omega_{\rm p}$ | degrees |
| Inclination | $i$ | degrees |
| Stellar radius | $R_*$ | solar radii |
| Stellar masses | $M_1,M_2,M_3$ | solar masses |
| Wind velocity | $v_\infty$ | km s$^{-1}$ |
| Mass-loss rate | $\dot M$ | $M_\odot\,{\rm yr}^{-1}$ |
| Energy | $E$ | keV |
| Wavelength | $\lambda$ | angstrom |
| Hydrogen column | $N_{\rm H}$ | $10^{22}\,{\rm cm}^{-2}$ |

### Relative and barycentric semi-major axes

In the binary models, `semimajor` is the **relative semi-major axis**

$$
a_{\rm rel}=a_1+a_2,
$$

expressed in units of the donor radius $R_*$. If `Mstar1` is the object carrying the observed feature and `Mstar2` is its companion, the barycentric semi-major axis of `Mstar1` is

$$
a_1=a_{\rm rel}\frac{M_2}{M_1+M_2}.
$$

This distinction is important:

- Doppler motion of an emitter uses its **barycentric** orbit.
- Wind density and absorption use the **relative separation** from the donor centre.

---

## Orbital motion and the phase--time relation

### Eccentric separation

Let

$$
\theta=2\pi\phi,
\qquad
\nu=\theta-\omega_{\rm p},
$$

where $\nu$ is the true-anomaly-like angle used by the model. The relative separation is

$$
r_{\rm rel}(\nu)=
\frac{a_{\rm rel}(1-e^2)}{1+e\cos\nu}.
$$

When `semimajor` is given in donor-radius units,

$$
r_{\rm rel,cm}=r_{\rm rel}R_*.
$$

### Kepler's second law

For an eccentric orbit, geometrical phase does not increase linearly with time. The package uses Kepler's second law,

$$
dA=\frac{1}{2}r^2d\theta,
\qquad
\frac{dA}{dt}=\mathrm{constant},
$$

to convert between orbital phase and time. The swept area is integrated on a numerical phase grid and normalized so that one complete revolution equals the supplied orbital period.

The conversion also returns the instantaneous angular velocity

$$
W=\frac{d\theta}{dt},
$$

in rad s$^{-1}$. An equivalent analytic Keplerian expression is

$$
\frac{d\nu}{dt}
=
\frac{2\pi}{P_{\rm orb}}
\frac{(1+e\cos\nu)^2}{(1-e^2)^{3/2}}.
$$

The `precision` argument controls the internal phase-grid spacing. It is a numerical resolution parameter, not a statistical uncertainty.

---

## Doppler modulation

### First-order Doppler transformation

The package uses the non-relativistic approximation

$$
D=1+\frac{v_{\rm los}}{c},
$$

where $v_{\rm los}$ is the model line-of-sight velocity. The observable transforms as

$$
E_{\rm obs}=\frac{E_0}{D},
\qquad
\lambda_{\rm obs}=\lambda_0D,
\qquad
P_{\rm obs}=P_0D.
$$

This approximation is appropriate when $|v_{\rm los}|\ll c$.

### Eccentric binary orbit

For the barycentric orbit of the emitting object,

$$
r_1=\frac{a_1(1-e^2)}{1+e\cos\nu}.
$$

The radius changes in an eccentric orbit, so the projected velocity must include both radial and tangential terms. With the package projection

$$
z=r_1\cos\theta\sin i,
$$

the complete derivative is

$$
v_{\rm orb,los}
=
\frac{dz}{dt}
=
\left(
\dot r_1\cos\theta-r_1\dot\theta\sin\theta
\right)\sin i.
$$

The radial derivative is

$$
\dot r_1
=
\frac{a_1(1-e^2)e\sin\nu}
{(1+e\cos\nu)^2}\,\dot\nu.
$$

The first term vanishes only for a circular orbit. The corrected orbital models retain both contributions.

An optional phenomenological flow term can be added:

$$
v_{\rm flow,los}=v_{\rm wind}\cos\theta\sin i.
$$

This extra term is not the same as the accelerating stellar-wind solution used in the density and absorption models.

### Logarithmic spiral

The isolated spiral is described by

$$
R(\theta)=a_{\rm s}e^{b\theta}.
$$

The spiral phase evolves as

$$
\phi_{\rm s}(t)=\phi_{{\rm s},0}+\Omega_{\rm s}(t-t_0),
$$

where `omega` is expressed in **cycles per second**. Therefore,

$$
\theta_{\rm s}=2\pi\phi_{\rm s},
\qquad
\dot\theta_{\rm s}=2\pi\Omega_{\rm s}.
$$

Because $R$ changes with time,

$$
\dot R=bR\dot\theta_{\rm s}.
$$

For the projected coordinate $z=R\cos\theta_{\rm s}\sin i_{\rm s}$, the full line-of-sight velocity is

$$
v_{\rm s,los}
=
R\dot\theta_{\rm s}
\left[b\cos\theta_{\rm s}-\sin\theta_{\rm s}\right]
\sin i_{\rm s}.
$$

The model is kinematic: it describes a trajectory but does not solve the dynamics that create the spiral.

### Spiral plus binary orbit

The combined model adds the two projected velocities,

$$
v_{\rm los}=v_{\rm orbit,los}+v_{\rm spiral,los},
$$

and applies a single Doppler transformation to the total velocity. The orbital and spiral components can have different phases and inclinations.

### Hierarchical orbit or disc model

The disc model contains two Keplerian components:

1. The centre of mass of the inner pair $(M_1+M_3)$ follows the outer binary orbit around the donor $M_2$.
2. `Mass3`, which carries the observed feature, follows an inner orbit around `Mstar1`.

The total projected velocity is

$$
v_{\rm los}=v_{\rm outer,los}+v_{\rm inner,los}.
$$

This is a model for a coherently moving emitting region. It is not a continuous accretion-disc line-profile model and does not include shear, turbulence, emissivity gradients, or radiative transfer.

---

## Stellar-wind density: CAK-like beta law

The wind implementation uses a **CAK-like beta velocity law**,

$$
v(r)=v_\infty\left(1-\frac{R_*}{r}\right)^\beta,
\qquad r>R_*.
$$

It is important to distinguish this prescription from a full solution of the Castor--Abbott--Klein line-driving equations. The package adopts the commonly used beta-law parameterization but does not solve the line force self-consistently.

For a stationary, smooth, spherically symmetric wind, mass conservation gives

$$
\dot M=4\pi r^2\rho(r)v(r),
$$

and hence

$$
\rho(r)=\frac{\dot M}{4\pi r^2v(r)}.
$$

The local density encountered by the compact object is evaluated at the relative orbital separation from the donor. No barycentric mass correction is applied to this distance.

The equivalent hydrogen number density is

$$
n_{\rm H}(r)=\frac{X_{\rm H}\rho(r)}{m_{\rm p}},
$$

where $X_{\rm H}$ is the hydrogen mass fraction. The theoretical density and ionization functions allow this value to be supplied explicitly; the current fitting kernel uses $X_{\rm H}=0.70$.

---

## Equivalent hydrogen column density

At a given orbital phase, the line of sight is parameterized by the distance $z$ from the compact object towards the observer. If $r_{\rm orb}$ is the donor--compact-object separation and $\alpha$ is the angle between the orbital radius vector and the line of sight, the distance from a point on the ray to the donor centre is

$$
x(z)=
\sqrt{r_{\rm orb}^2+z^2-2r_{\rm orb}z\cos\alpha}.
$$

The package geometry uses

$$
\cos\alpha=\cos(2\pi\phi)\sin i.
$$

The equivalent hydrogen column is

$$
N_{\rm H}(\phi)
=
\int_0^{z_{\rm max}} n_{\rm H}[x(z)]\,dz
=
\frac{X_{\rm H}}{m_{\rm p}}
\int_0^{z_{\rm max}}\rho[x(z)]\,dz.
$$

The corrected models use adaptive numerical quadrature and integrate to $z_{\rm max}=1000R_*$. Results are returned in units of $10^{22}\,{\rm cm}^{-2}$.

### Eclipse convention

If the line of sight intersects the stellar photosphere, the current package convention returns

$$
N_{\rm H}=0.
$$

This value is an occultation flag expressed through the historical output convention. It does **not** mean that the physical column is zero. A future API would be clearer if it returned a separate eclipse mask.

---

## Wind ionization

The ionization parameter is

$$
\xi=\frac{L_{\rm X}}{n_{\rm H}d^2},
$$

with units of erg cm s$^{-1}$, where $d$ is the distance to the ionizing source.

The one-dimensional line-of-sight calculation returns $\log_{10}\xi$, together with $n_{\rm H}$ and the remaining hydrogen column beyond each position,

$$
N_{\rm H}(>z)=\int_z^{z_{\rm max}}n_{\rm H}(z')\,dz'.
$$

The two-dimensional ionization map evaluates $\xi$ on a polar grid centred on the donor and includes a geometrical stellar shadow. It does not calculate attenuation of the ionizing continuum through the wind before reaching each grid cell.

!!! warning "Luminosity units"
    The one-dimensional ionization profile retains the historical input convention of luminosity in units of $10^{32}\,{\rm erg\,s^{-1}}$. `ionization_map_phase` accepts luminosity directly in ${\rm erg\,s^{-1}}$. This difference should be kept explicit in scripts until the API is unified.

The ionization routines do not solve thermal balance, ionic populations, recombination, or non-equilibrium effects. They should not be interpreted as replacements for XSTAR, Cloudy, SPEX, XSPEC photoionization tables, or another radiative-transfer calculation.

---

## Timing analysis

### Hardness and colour ratios

For hard- and soft-band values $H$ and $S$, the hardness ratio is

$$
HR=\frac{H-S}{H+S}.
$$

Assuming independent one-sigma uncertainties, error propagation gives

$$
\sigma_{HR}=
\sqrt{
\left(\frac{2S}{(H+S)^2}\sigma_H\right)^2+
\left(\frac{2H}{(H+S)^2}\sigma_S\right)^2
}.
$$

The colour ratio is

$$
CR=\frac{H}{S},
$$

with

$$
\sigma_{CR}=
\sqrt{
\left(\frac{\sigma_H}{S}\right)^2+
\left(\frac{H\sigma_S}{S^2}\right)^2
}.
$$

These expressions assume no covariance between the two bands.

### Inverse-variance rebinning

For measurements $y_i$ with uncertainties $\sigma_i$, the package combines points with weights

$$
w_i=\frac{1}{\sigma_i^2},
$$

so that

$$
\bar y=\frac{\sum_iw_iy_i}{\sum_iw_i},
\qquad
\sigma_{\bar y}=\sqrt{\frac{1}{\sum_iw_i}}.
$$

`rebin_snr` accumulates consecutive points until

$$
\frac{|\bar y|}{\sigma_{\bar y}}\geq {\rm S/N}_{\rm min}
$$

and at least `min_bins` input points have been collected. The current convention is the usual **signal divided by noise**; for example, `min_snr=5` means S/N $\geq5$.

`rebin_resolution` uses the temporal overlap between input and output bins and preserves the time-integrated quantity. For rates or fluxes,

$$
I=\sum_i y_i\Delta t_i,
\qquad
\bar y=\frac{I}{\sum_i\Delta t_i},
$$

and, for independent errors,

$$
\sigma_I=\sqrt{\sum_i(\sigma_i\Delta t_i)^2},
\qquad
\sigma_{\bar y}=\frac{\sigma_I}{\sum_i\Delta t_i}.
$$

For gapped light curves, real bin widths or exposures should be supplied rather than inferred only from bin centres.

### Pulse folding

For a trial period $P$, pulse phase is

$$
\phi_{\rm pulse}=
\frac{(t-t_0)\bmod P}{P}.
$$

The package sorts the light curve by pulse phase and can return the unbinned folded data or rebin the profile by S/N or by a fixed number of phase-sorted points. Folding does not itself account for period derivatives, binary demodulation, or a changing pulse profile.

### Lomb--Scargle period search

The sliding-window routine computes an uncertainty-weighted Lomb--Scargle periodogram inside successive time intervals. `window_sec` and `step_sec` are physical time durations in the same units as the input times.

Candidate peaks are retained when their false-alarm probability satisfies the user threshold. The frequency uncertainty is estimated from the local profile-(\chi^2) curve around each Lomb--Scargle peak. After refining the peak frequency (f_{\rm best}), the lower and upper frequency bounds, (f_-) and (f_+), are defined by

[
\chi^2(f)-\chi^2(f_{\rm best})=\Delta\chi^2,
]

with (\Delta\chi^2=1) by default. Under the usual Gaussian and single-parameter assumptions, this approximately corresponds to a local (68.3%) confidence interval. The frequency uncertainty is estimated from the local profile-(\chi^2) curve around each Lomb--Scargle peak. After refining the peak frequency (f_{\rm best}), the lower and upper frequency bounds, (f_-) and (f_+), are defined by

$$
\chi^2(f)-\chi^2(f_{\rm best})=\Delta\chi^2,
$$

with (\Delta\chi^2=1) by default. Under the usual Gaussian and single-parameter assumptions, this approximately corresponds to a local (68.3%) confidence interval. The frequency uncertainty is estimated from the local profile-(\chi^2) curve around each Lomb--Scargle peak. After refining the peak frequency (f_{\rm best}), the lower and upper frequency bounds, (f_-) and (f_+), are defined by

$$
\chi^2(f)-\chi^2(f_{\rm best})=\Delta\chi^2,
$$

with (\Delta\chi^2=1) by default. Under the usual Gaussian and single-parameter assumptions, this approximately corresponds to a local (68.3%) confidence interval. The resulting frequency uncertainties are therefore generally asymmetric,

$$
\sigma_{f,-}=f_{\rm best}-f_-,
$$

and

$$
\sigma_{f,+}=f_+-f_{\rm best}.
$$

The corresponding period bounds are obtained directly from the nonlinear transformation (P=1/f),

$$
P_{\rm lower}=\frac{1}{f_+},
$$

and

$$
P_{\rm upper}=\frac{1}{f_-}.
$$

This gives asymmetric period uncertainties,

$$
\sigma_{P,-}=P_{\rm best}-P_{\rm lower},
$$

and

$$
\sigma_{P,+}=P_{\rm upper}-P_{\rm best}.
$$

The periodogram power S/N is estimated robustly as

$$
{\rm S/N}_{\rm power}
=====================

\frac{P_{\rm peak}-\operatorname{median}(P)}
{1.4826,\operatorname{MAD}(P)}.
$$

These profile-(\chi^2) uncertainties are local confidence estimates around a selected peak. They do not account for ambiguity between aliases or multiple competing peaks. Likewise, the standard Lomb--Scargle false-alarm probability does not automatically include red-noise effects, the full number of sliding-window trials, or subsequent selection cuts.




---

## Model fitting

### Weighted objective function

For observed values $y_i$, model values $m_i(\boldsymbol\theta)$, and uncertainties $\sigma_i$, the fitting routines minimize

$$
\chi^2(\boldsymbol\theta)
=
\sum_i
\left[
\frac{y_i-m_i(\boldsymbol\theta)}{\sigma_i}
\right]^2.
$$

A standard statistical interpretation of $\chi^2$ requires appropriate, independent uncertainties and a suitable model. If no measurement errors are supplied, the package creates working weights, but the resulting absolute $\chi^2$ should not be interpreted as a calibrated goodness-of-fit probability.

### Discrete and extended observations

In `discrete` mode, each datum is represented by one time:

$$
m_i=m(t_i).
$$

In `extended` mode, each datum covers a time interval $[t_{i,0},t_{i,1}]$. The physically relevant prediction is the bin-averaged observable,

$$
\bar m_i=
\frac{1}{t_{i,1}-t_{i,0}}
\int_{t_{i,0}}^{t_{i,1}}m(t)\,dt.
$$

The implementation approximates this integral by sampling and averaging the **final observable**. Narrow bins, as defined by `extended_binsize` in orbital-phase cycles, are evaluated at the midpoint.

---

## Particle Swarm Optimization

Particle Swarm Optimization (PSO) is a derivative-free global optimization method. It is useful for the strongly non-linear, bounded, and sometimes multi-modal models in `xraybinaryorbit`.

### What is a particle?

A particle is one complete candidate solution. For a model with $D$ fitted parameters, particle $k$ has a position

$$
\mathbf{x}_k=(\theta_{1,k},\theta_{2,k},\ldots,\theta_{D,k}).
$$

A particle is **not** one parameter. For example, in `fit_orbit_ps`, one particle simultaneously contains a proposed phase, semi-major axis, period, eccentricity, periapsis, inclination, masses, wind velocity, and rest feature.

Particles are initialized within the lower and upper bounds. At every swarm update, each particle is evaluated through the full physical model and assigned an objective value, here the weighted $\chi^2$.

### How particles move

In standard PSO, particle velocity and position are updated schematically as

$$
\mathbf{v}_k^{(j+1)}=
\omega\mathbf{v}_k^{(j)}
+c_{\rm p}\mathbf{r}_{\rm p}
\left(\mathbf{p}_k-\mathbf{x}_k^{(j)}\right)
+c_{\rm g}\mathbf{r}_{\rm g}
\left(\mathbf{g}-\mathbf{x}_k^{(j)}\right),
$$

$$
\mathbf{x}_k^{(j+1)}=
\mathbf{x}_k^{(j)}+\mathbf{v}_k^{(j+1)}.
$$

Here:

- $\mathbf{p}_k$ is the best position previously found by that particle;
- $\mathbf{g}$ is the best position found by the entire swarm;
- $\omega$ controls inertia;
- $c_{\rm p}$ and $c_{\rm g}$ control attraction to the personal and global best positions;
- $\mathbf{r}_{\rm p}$ and $\mathbf{r}_{\rm g}$ are random vectors.

### Meaning of the package controls

| Argument | Meaning |
|---|---|
| `swarmsize` | Number of particles in one swarm. More particles sample parameter space more densely but increase the cost of every update. |
| `maxiter` | Maximum number of internal swarm updates in one PSO run. It is not the number of independent fits. |
| `num_iterations` | Number of complete, independent PSO runs launched by `xraybinaryorbit`. Each run starts a new swarm. |

An approximate upper scale for the number of model evaluations is

$$
N_{\rm eval}\sim
N_{\rm runs}\,N_{\rm particles}\,(N_{\rm updates}+1),
$$

although early stopping can reduce it. This explains why increasing all three controls simultaneously can make a fit much slower, especially for $N_{\rm H}$, where each objective evaluation contains numerical line-of-sight integrations.

### Which PSO result is returned?

For each independent run, the package stores the best parameter vector and its $\chi^2$. The run with the lowest $\chi^2$ is selected as the reported best solution.

The `Std` row in the returned table is

$$
\sigma_{\rm run}(\theta_j)
=
\operatorname{std}
\left(
\theta_{j,1},\ldots,\theta_{j,N_{\rm runs}}
\right),
$$

the dispersion of the best parameter from each independent run.

!!! warning "PSO `Std` is not a confidence interval"
    This value measures optimizer stability under the chosen bounds, swarm size, stopping criteria, and random initialization. It is not a posterior standard deviation, a one-sigma observational uncertainty, or a profile-likelihood interval. With only a few independent runs, it can be especially unstable.

A small PSO `Std` means that the optimizer repeatedly found similar solutions. It does not prove that the parameter is physically well constrained. Conversely, a large value may reveal multiple minima, insufficient swarm exploration, or a genuine parameter degeneracy.

### Recommended PSO checks

- Use physically justified bounds.
- Run enough independent swarms to examine stability.
- Inspect whether best-fit values lie on a bound.
- Compare the best $\chi^2$ values from all runs.
- Plot residuals and the model, not only the parameter table.
- For scientific confidence intervals, use profile likelihoods, bootstrap simulations, or posterior sampling around the identified solution.

---

## Bounded least squares and its errors

The `_ls` functions use bounded non-linear least squares through `scipy.optimize.curve_fit`. The initial point is the midpoint of the supplied bounds.

When measurement uncertainties are supplied, they are passed as `sigma` with `absolute_sigma=True`. The parameter covariance matrix is approximated locally from the Jacobian near the optimum, and the reported errors are

$$
\sigma_{\theta_j}=\sqrt{C_{jj}},
$$

where $C$ is the covariance matrix.

These errors rely on a local linearization and are most meaningful when:

- the model is smooth near the optimum;
- the optimum is not on a bound;
- the parameters are locally identifiable;
- residuals are compatible with the assumed uncertainties;
- the likelihood is approximately Gaussian near the minimum.

They can be misleading for multi-modal, highly degenerate, discontinuous, or boundary-limited problems. The LS routines return $R^2$ as a descriptive statistic,

$$
R^2=1-
\frac{\sum_i(y_i-m_i)^2}
{\sum_i(y_i-\bar y)^2},
$$

but $R^2$ is not a replacement for $\chi^2$, residual analysis, or physical validation.

### PSO and LS errors are different quantities

| Output | Meaning |
|---|---|
| PSO `Std` | Run-to-run dispersion of independently optimized best values |
| LS `Std` | Local one-sigma approximation from the covariance matrix |

They should not be compared as though they were calculated from the same statistical definition.

---

## Main physical degeneracies

### Doppler orbit

The velocity amplitude scales approximately as

$$
K\propto\frac{a_1\sin i}{P_{\rm orb}},
$$

so semi-major axis and inclination can be strongly correlated. Initial phase, periapsis, and eccentricity can also compensate for each other when orbital coverage is incomplete. The rest feature can correlate with the mean offset, while `wind_vel` can absorb part of a sinusoidal component.

### Spiral models

`semimajor_spiral`, `b`, and `omega` are correlated because they jointly control the spatial scale, expansion, phase evolution, and velocity amplitude. A short observation may not distinguish an expanding spiral from a sinusoid with a slowly changing amplitude.

### Hierarchical orbit

Two periodic components can generate aliases and exchange amplitude when their periods are similar or when only a small fraction of the outer orbit is covered. `Mass3` can remain weakly constrained if both orbital scales are also free.

### Absorption

To first order,

$$
N_{\rm H}\propto\frac{\dot M}{v_\infty},
$$

so mass-loss rate and terminal velocity are strongly degenerate. Inclination, separation, eccentricity, periapsis, and $\beta$ can also trade off. With only a few $N_{\rm H}$ measurements, fitting all wind and orbital parameters simultaneously is generally not identifiable.

---

## Scope and limitations

The orbital models assume fixed Keplerian motion and do not include apsidal precession, relativistic delays, secular period evolution, Roche distortion, or dynamical mass transfer.

The wind model assumes a smooth, stationary, spherically symmetric flow. It does not self-consistently include clumping, porosity, an accretion wake, focused streams, radiative inhibition, ionization feedback on acceleration, shocks, or magnetic structure.

The spiral and disc models are phenomenological kinematic descriptions. A good fit alone does not establish that the corresponding hydrodynamic structure exists.

The timing tools assume that the supplied uncertainties are meaningful. Lomb--Scargle peak FAP values and peak-width errors can be optimistic when the light curve contains red noise, strong non-stationarity, or extensive multiple testing.

---
