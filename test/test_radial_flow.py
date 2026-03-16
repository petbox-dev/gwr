"""
Regression tests for infinite-acting radial flow with wellbore storage.

Uses besselk(0, sqrt(s)) with cD=1e3, skin=0.
Compares against pre-computed reference data in reference_data.npz.
Regenerate with:  python test/generate_reference.py
"""

from pathlib import Path

import numpy as np
from mpmath import mp

from gwr_inversion import gwr

REF = np.load(Path(__file__).parent / 'reference_data.npz')
ATOL = 1e-12


def _lap_rad(s, cD=1e3, skin=0.0):
    b = mp.besselk(0, mp.sqrt(s))
    return 1.0 / (s / (b + skin) + (s ** 2) * cD)


class TestRadialFlow:

    def test_pressure(self):
        result = np.array([float(v) for v in gwr(
            lambda s: _lap_rad(s, cD=1e3), REF['time_rad'], 8
        )])
        np.testing.assert_allclose(result, REF['rad_pressure'], atol=ATOL)

    def test_derivative(self):
        time = REF['time_rad']
        result = np.array([float(v) for v in (
            time * gwr(lambda s: s * _lap_rad(s, cD=1e3), time, 8)
        )])
        np.testing.assert_allclose(result, REF['rad_derivative'], atol=ATOL)

    def test_integral(self):
        time = REF['time_rad']
        result = np.array([float(v) for v in (
            1 / time * gwr(lambda s: 1 / s * _lap_rad(s, cD=1e3), time, 8)
        )])
        np.testing.assert_allclose(result, REF['rad_integral'], atol=ATOL)
