Welcome to xraybinaryorbit's documentation!
============

**xraybinaryorbit** is a user-friendly Python package designed to facilitate the analysis of orbital modulations in X-ray binaries.

**Overview**

Orbital modulations in X-ray binaries are widely known but complex to analyze due to their dependence on several parameters and geometrical considerations. To address these challenges, we have collected various functions developed over years of analyzing close X-ray binaries and combined them into a Python package useful for almost every X-ray binary analysis. The aim is to simplify the implementation of these analyses for other astronomers.

Given the complexity and the number of parameters involved, we propose a user-friendly "form" method to enhance the package's usability. These forms, once loaded and submitted, are saved for future interactions, allowing for a more streamlined workflow.

Each function is briefly described when loaded, providing guidance for its execution.

**Underlying Physics**

The primary physics models and principles used in this package include:

- The Doppler effect
- Kepler's laws
- The CAK model (Castor, Abbott, and Klein), as described in:
  
  Castor, J. I., Abbott, D. C., & Klein, R. I. (1975). *Radiation-driven winds in Of stars*. Astrophysical Journal, 195, 157-174.

**Dependencies**

To function correctly, this package requires the following Python libraries:

- `numpy`
- `pandas`
- `scipy`
- `pyswarm`
- `matplotlib`
- `os`
- `warnings`
- `tkinter`
- `astropy`
- `inspect`

These dependencies support numerical computations, data manipulation, scientific computing, optimization, visualization, GUI creation, time-series analysis, and interaction with the operating system. Please ensure these packages are installed to use **xraybinaryorbit** successfully.

In out github repository https://github.com/xragua/xraybinaryorbit some practical examples and other usefull info can be found.
