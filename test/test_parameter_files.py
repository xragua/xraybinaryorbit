"""Verify that every tested load_directly function has a valid TXT file."""

from pathlib import Path

import pytest

TEST_DIR = Path(__file__).resolve().parent

EXPECTED_HEADERS = {
    "bounds_orbit.txt": ["iphase","semimajor","orbitalperiod","eccentricity","periapsis","inclination","Rstar","Mstar1","Mstar2","wind_vel","feature"],
    "bounds_spiral.txt": ["iphase_spiral","semimajor_spiral","b","omega","inclination_spiral","feature"],
    "bounds_spiral_orbit.txt": ["iphase_orbit","semimajor_orbit","orbitalperiod","eccentricity","periapsis","inclination","Rstar","Mstar1","Mstar2","iphase_spiral","semimajor_spiral","b","omega","inclination_spiral","feature"],
    "bounds_disc.txt": ["iphase","semimajor","orbitalperiod","eccentricity","periapsis","inclination","Rstar","Mstar1","Mstar2","iphase2","semimajor2","orbitalperiod2","eccentricity2","periapsis2","inclination2","Mass3","feature","wind_vel"],
    "bounds_nh.txt": ["iphase","semimajor","orbitalperiod","eccentricity","periapsis","inclination","Rstar","Mdot","v_inf","beta"],
    "orbit.txt": ["iphase","semimajor","orbitalperiod","eccentricity","periapsis","inclination","Rstar","Mstar1","Mstar2","wind_vel","feature"],
    "spiral.txt": ["iphase_spiral","semimajor_spiral","b","omega","inclination_spiral","feature"],
    "disc.txt": ["iphase","semimajor","orbitalperiod","eccentricity","periapsis","inclination","Rstar","Mstar1","Mstar2","iphase2","semimajor2","orbitalperiod2","eccentricity2","periapsis2","inclination2","Mass3","wind_vel","feature"],
    "spiral_in_orbit.txt": ["iphase","semimajor","orbitalperiod","eccentricity","periapsis","inclination","iphase_spiral","semimajor_spiral","b","omega","inclination_spiral","Rstar","Mstar1","Mstar2","wind_vel","feature"],
    "density_through_orbit.txt": ["semimajor","orbitalperiod","eccentricity","periapsis","Rstar","wind_infinite_velocity","Mass_loss_rate","beta"],
    "absorption_column_through_orbit.txt": ["semimajor","orbitalperiod","eccentricity","periapsis","inclination","Rstar","wind_infinite_velocity","Mass_loss_rate","beta"],
    "den_chi_orbphase.txt": ["orb_phase","luminosity","semimajor","eccentricity","periapsis","inclination","Rstar","wind_infinite_velocity","Mass_loss_rate","beta"],
    "ionization_map_phase.txt": ["phase","semimajor","eccentricity","periapsis","Rstar","wind_infinite_velocity","Mass_loss_rate","beta","luminosity","bound1","bound2"],
}


@pytest.mark.parametrize("filename,expected", EXPECTED_HEADERS.items())
def test_parameter_file_exists_and_matches_current_header(filename, expected):
    path = TEST_DIR / filename
    assert path.is_file(), f"Missing required parameter file: {filename}"
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    expected_line_count = 3 if filename.startswith("bounds_") else 2
    assert len(lines) == expected_line_count
    assert lines[0].split(",") == expected
    for values_line in lines[1:]:
        values = [float(value) for value in values_line.split(",")]
        assert len(values) == len(expected)
