# __init__.py in your submodule

# Import shared resources from main.py
from ..main import *  # Import everything from main.py (make sure this is what you want)

# Import everything from submodules
from .fit_conic_orbit import *
from .fit_disc import *
from .fit_nh import *
from .fit_spiral import *
from .fit_spiral_orbit import *

from ..helpers.data_helpers import (
    _define_x_y_sy,
    _copy_fields,
    _load_values_to_interface,
    _manage_parameters,
    _load_bounds_to_interface,
    _manage_bounds
)

from ..helpers.math_helpers import (
    _gaussian,
    _time_pairs,
    _interpolate_pchip,
    _chi_squared_weighted,
    _chi_squared,
    _orbital_phase_to_time,
    _orbital_time_to_phase
)




