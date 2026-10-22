"""
Regression tests for unit step inversion: H(t-1) via e^(-s)/s.

Compares against pre-computed reference data in reference_data.npz.
Regenerate with:  python test/generate_reference.py
"""

from pathlib import Path

import numpy as np
from mpmath import mp

from gwr_inversion import gwr

from backend_markers import pins_gmpy2_rounding

REF = np.load(Path(__file__).parent / 'reference_data.npz')
ATOL = 1e-12


class TestUnitStep:

    @pins_gmpy2_rounding
    def test_m8(self):
        result = np.array([float(v) for v in gwr(
            lambda s: 1 / s * mp.exp(-1.0 * s), REF['time_unit'], 8
        )])
        np.testing.assert_allclose(result, REF['unit_m8'], atol=ATOL)

    def test_m128(self):
        result = np.array([float(v) for v in gwr(
            lambda s: 1 / s * mp.exp(-1.0 * s), REF['time_unit'], 128
        )])
        np.testing.assert_allclose(result, REF['unit_m128'], atol=ATOL)
