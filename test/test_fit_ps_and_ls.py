"""Fast deterministic contract tests for wrappers using TXT bounds only."""

from __future__ import annotations

import inspect
from dataclasses import dataclass

import numpy as np
import pandas as pd
import pytest

import xraybinaryorbit as xbo


@dataclass(frozen=True)
class WrapperCase:
    ps_name: str
    ls_name: str | None
    model_name: str
    phase_name: str
    lower: np.ndarray
    upper: np.ndarray
    feature_index: int | None
    return_count: int


CASES = [
    WrapperCase(
        "fit_orbit_ps", "fit_orbit_ls", "conic_orbit", "_fitted_phase",
        np.array([0.0,1.0,1.0,0.0,0.0,10.0,1.0,0.1,0.1,-10.0,5.0]),
        np.array([1.0,3.0,5.0,0.8,360.0,90.0,30.0,30.0,30.0,10.0,10.0]),
        10, 4,
    ),
    WrapperCase(
        "fit_spiral_ps", "fit_spiral_ls", "spiral", "_fitted_phase",
        np.array([0.0,0.01,-0.1,1e-5,10.0,5.0]),
        np.array([1.0,3.0,0.1,1e-3,90.0,10.0]),
        5, 4,
    ),
    WrapperCase(
        "fit_spiral_in_orbit_ps", "fit_spiral_in_orbit_ls", "spiral_orbit", "_fitted_orbit_phase",
        np.array([0.0,1.0,1.0,0.0,0.0,10.0,1.0,0.1,0.1,0.0,0.01,-0.1,1e-5,10.0,5.0]),
        np.array([1.0,3.0,5.0,0.8,360.0,90.0,30.0,30.0,30.0,1.0,3.0,0.1,1e-3,90.0,10.0]),
        14, 4,
    ),
    WrapperCase(
        "fit_disc_ps", "fit_disc_ls", "disc_in_orbit", "_fitted_phases",
        np.array([0.0,1.0,1.0,0.0,0.0,10.0,1.0,0.1,0.1,0.0,0.1,0.1,0.0,0.0,10.0,0.0,5.0,-10.0]),
        np.array([1.0,3.0,5.0,0.8,360.0,90.0,30.0,30.0,30.0,1.0,2.0,2.0,0.8,360.0,90.0,2.0,10.0,10.0]),
        16, 5,
    ),
    WrapperCase(
        "fit_nh_ps", None, "nh_orbit", "_fitted_phase",
        np.array([0.0,1.1,1.0,0.0,0.0,10.0,1.0,1e-8,500.0,0.0]),
        np.array([1.0,4.0,5.0,0.8,360.0,90.0,30.0,1e-5,2500.0,2.0]),
        None, 4,
    ),
]


def _constant_model_for_case(case, midpoint):
    value = midpoint[0] if case.feature_index is None else midpoint[case.feature_index]

    def model(x, *params):
        length = np.asarray(x).shape[0]
        selected = params[0] if case.feature_index is None else params[case.feature_index]
        return np.full(length, selected, dtype=float)

    return value, model


@pytest.mark.parametrize("case", CASES, ids=lambda case: case.ps_name)
def test_pso_wrappers_load_bounds_from_txt_without_gui(monkeypatch, case):
    function = getattr(xbo, case.ps_name)
    module = inspect.getmodule(function)
    midpoint = 0.5 * (case.lower + case.upper)
    expected_value, fake_model = _constant_model_for_case(case, midpoint)

    def fake_pso(objective, lb, ub, **kwargs):
        assert np.allclose(lb, case.lower)
        assert np.allclose(ub, case.upper)
        params = 0.5 * (np.asarray(lb) + np.asarray(ub))
        return params, objective(params)

    monkeypatch.setattr(module, "pso", fake_pso)
    monkeypatch.setattr(module, case.model_name, fake_model)

    x = np.array([0.0, 1.0, 2.0])
    y = np.full(x.size, expected_value)
    sigma = np.full(x.size, 0.1)

    if case.phase_name == "_fitted_phases":
        monkeypatch.setattr(module, case.phase_name, lambda x_data, params, method: (
            np.arange(len(x_data), dtype=float),
            np.arange(len(x_data), dtype=float) + 0.5,
        ))
    else:
        monkeypatch.setattr(module, case.phase_name, lambda x_data, params, method: np.arange(len(x_data), dtype=float))

    output = function(
        x, y, y_err=sigma,
        num_iterations=1, maxiter=2, swarmsize=2,
        method_="discrete",
        load_directly=True,
    )

    assert len(output) == case.return_count
    frame, predicted, score = output[0], output[-2], output[-1]
    assert isinstance(frame, pd.DataFrame)
    assert np.allclose(frame.loc["Value"].to_numpy(dtype=float), midpoint)
    assert np.allclose(frame.loc["Std"].to_numpy(dtype=float), 0.0)
    assert np.allclose(predicted, y)
    assert score == pytest.approx(0.0)


@pytest.mark.parametrize("case", [c for c in CASES if c.ls_name], ids=lambda case: case.ls_name)
def test_least_squares_wrappers_load_bounds_from_txt_without_gui(monkeypatch, case):
    function = getattr(xbo, case.ls_name)
    module = inspect.getmodule(function)
    midpoint = 0.5 * (case.lower + case.upper)
    expected_value, fake_model = _constant_model_for_case(case, midpoint)

    def fake_curve_fit(model, x, y, p0, bounds, **kwargs):
        lower, upper = bounds
        assert np.allclose(lower, case.lower)
        assert np.allclose(upper, case.upper)
        assert np.allclose(p0, midpoint)
        assert np.asarray(model(x, *p0)).shape == np.asarray(y).shape
        return np.asarray(p0), np.eye(len(p0)) * 0.04

    monkeypatch.setattr(module, "curve_fit", fake_curve_fit)
    monkeypatch.setattr(module, case.model_name, fake_model)

    x = np.array([0.0, 1.0, 2.0])
    y = np.full(x.size, expected_value)
    sigma = np.full(x.size, 0.1)

    if case.phase_name == "_fitted_phases":
        monkeypatch.setattr(module, case.phase_name, lambda x_data, params, method: (
            np.arange(len(x_data), dtype=float),
            np.arange(len(x_data), dtype=float) + 0.5,
        ))
    else:
        monkeypatch.setattr(module, case.phase_name, lambda x_data, params, method: np.arange(len(x_data), dtype=float))

    output = function(
        x, y, y_err=sigma,
        method_="discrete",
        load_directly=True,
    )

    assert len(output) == case.return_count
    frame, predicted, r_squared = output[0], output[-2], output[-1]
    assert isinstance(frame, pd.DataFrame)
    assert np.allclose(frame.loc["Value"].to_numpy(dtype=float), midpoint)
    assert np.allclose(frame.loc["Std"].to_numpy(dtype=float), 0.2)
    assert np.allclose(predicted, y)
    assert np.isnan(r_squared)


def test_fit_data_switches_one_dimensional_extended_input_to_discrete():
    module = inspect.getmodule(xbo.fit_orbit_ps)
    x = np.array([0.0, 1.0])
    y = np.array([6.4, 6.5])
    sigma = np.array([0.1, 0.1])

    with pytest.warns(RuntimeWarning, match="discrete"):
        prepared_x, prepared_y, prepared_sigma, method = module._prepare_fit_data(x, y, sigma, "extended")

    assert method == "discrete"
    assert np.array_equal(prepared_x, x)
    assert np.array_equal(prepared_y, y)
    assert np.array_equal(prepared_sigma, sigma)
