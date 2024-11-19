# __init__.py in your submodule

# Import shared resources from main.py
from ..main import *  # Import everything from main.py (make sure this is what you want)

# Import everything from submodules
from .fit_conic_orbit import _conic_orbit
from .fit_disc import _disc_in_orbit
from .fit_nh import _nh_orbit
from .fit_spiral import _spiral
from .fit_spiral_orbit import _spiral_orbit

# Define what should be accessible when using from fitting import *
__all__ = [
    '_conic_orbit',
    '_disc_in_orbit',
    '_nh_orbit',
    '_spiral',
    '_spiral_orbit'
]
