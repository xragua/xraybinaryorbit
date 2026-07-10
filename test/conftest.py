"""Shared configuration: headless plots and no interactive parameter GUI."""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib
import pytest

TEST_DIR = Path(__file__).resolve().parent
os.chdir(TEST_DIR)
matplotlib.use("Agg")


@pytest.fixture(autouse=True)
def forbid_parameter_gui(monkeypatch):
    """Fail clearly if a missing TXT file would otherwise open Tkinter."""
    try:
        from xraybinaryorbit.helpers import data_helpers
    except (ImportError, AttributeError):
        return

    def fail_gui(*args, **kwargs):
        pytest.fail(
            "A test attempted to open the Tkinter parameter interface. "
            "Check load_directly=True and the required TXT file in test/."
        )

    monkeypatch.setattr(
        data_helpers,
        "_load_values_to_interface",
        fail_gui,
        raising=False,
    )
    monkeypatch.setattr(
        data_helpers,
        "_load_bounds_to_interface",
        fail_gui,
        raising=False,
    )
