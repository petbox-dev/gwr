"""
Regression tests for ln(t) inversion: -(ln(s) + 0.577216) / s.

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


class TestLnT:

    @pins_gmpy2_rounding
    def test_m4(self):
        result = np.array([float(v) for v in gwr(
            lambda s: -mp.log(s) / s - 0.577216 / s, REF['time_ln'], 4
        )])
        np.testing.assert_allclose(result, REF['ln_m4'], atol=ATOL)
