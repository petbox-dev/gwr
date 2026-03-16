"""
Regression tests for hard oscillatory inversion: sin(t) via 1/(1+s^2).

Compares against pre-computed reference data in reference_data.npz.
Regenerate with:  python test/generate_reference.py
"""

from pathlib import Path

import numpy as np
from mpmath import mp

from gwr_inversion import gwr

REF = np.load(Path(__file__).parent / 'reference_data.npz')
ATOL = 1e-12


def _lap_sin(s, *args):
    return 1.0 / (1.0 + s * s)


class TestHardOscillatory:

    def test_m256(self):
        result = float(gwr(_lap_sin, REF['t0_hard'], 256)[0])
        assert abs(result - float(REF['sin_m256'][0])) < ATOL

    def test_m512(self):
        result = float(gwr(_lap_sin, REF['t0_hard'], 512)[0])
        assert abs(result - float(REF['sin_m512'][0])) < ATOL

    def test_m768(self):
        result = float(gwr(_lap_sin, REF['t0_hard'], 768)[0])
        assert abs(result - float(REF['sin_m768'][0])) < ATOL

    def test_m1024(self):
        result = float(gwr(_lap_sin, REF['t0_hard'], 1024)[0])
        assert abs(result - float(REF['sin_m1024'][0])) < ATOL
