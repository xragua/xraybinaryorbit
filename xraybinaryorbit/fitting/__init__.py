# __init__.py in your submodule

# Import shared resources from main.py
from ..main import *  # Import everything from main.py (make sure this is what you want)

# Import everything from submodules
from .fit_conic_orbit import *
from .fit_disc import *
from .fit_nh import *
from .fit_spiral import *
from .fit_spiral_orbit import *

