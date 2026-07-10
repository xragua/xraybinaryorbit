import os

import pandas as pd
import pytest

from xraybinaryorbit import ionization_map_phase


def _unpack_map_output(output):
    assert len(output) in (2, 3)
    if len(output) == 3:
        xi_map, electron_density_map, area = output
        assert isinstance(electron_density_map, pd.DataFrame)
    else:
        xi_map, area = output
    return xi_map, area


def test_ionization_map_phase_default():
    xi_map, area = _unpack_map_output(ionization_map_phase(load_directly=True))
    assert isinstance(xi_map, pd.DataFrame)
    assert not xi_map.empty
    assert isinstance(area, float)
    assert area >= 0


def test_ionization_map_phase_plot():
    filename = "test_ionization_map.png"
    _unpack_map_output(ionization_map_phase(
        save_plot=True, name="test_ionization_map", load_directly=True,
    ))
    assert os.path.isfile(filename)
    os.remove(filename)


def test_ionization_map_phase_custom_size_color():
    xi_map, area = _unpack_map_output(ionization_map_phase(
        size_in_Rstar=3.0, min_color=1e-5, max_color=1e2,
        load_directly=True,
    ))
    assert isinstance(xi_map, pd.DataFrame)
    assert not xi_map.empty
    assert area >= 0
