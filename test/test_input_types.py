"""
Tests for input type handling (scalar, list, ndarray).

Compares against pre-computed reference data in reference_data.npz.
Regenerate with:  python test/generate_reference.py
"""

from pathlib import Path

import numpy as np
from mpmath import mp

from gwr_inversion import gwr

REF = np.load(Path(__file__).parent / 'reference_data.npz')
ATOL = 1e-12


class TestInputTypes:

    def test_scalar(self):
        result = float(gwr(lambda s: -mp.log(s) / s - 0.577216 / s, 5.0, 4))
        assert abs(result - float(REF['scalar_result'][0])) < ATOL

    def test_list(self):
        result = [float(v) for v in gwr(
            lambda s: -mp.log(s) / s - 0.577216 / s, [1.0, 5.0, 10.0], 4
        )]
        np.testing.assert_allclose(result, REF['list_result'], atol=ATOL)

    def test_ndarray(self):
        result = np.array([float(v) for v in gwr(
            lambda s: -mp.log(s) / s - 0.577216 / s,
            np.array([1.0, 5.0, 10.0]), 4
        )])
        np.testing.assert_allclose(result, REF['list_result'], atol=ATOL)

    def test_scalar_returns_float(self):
        result = gwr(lambda s: -mp.log(s) / s - 0.577216 / s, 5.0, 4)
        assert isinstance(result, float)

    def test_scalar_returns_mpf(self):
        result = gwr(lambda s: -mp.log(s) / s - 0.577216 / s, 5.0, 4, dtype="mpf")
        assert isinstance(result, mp.mpf)

    def test_list_returns_list_of_float(self):
        result = gwr(lambda s: -mp.log(s) / s - 0.577216 / s, [1.0, 5.0], 4)
        assert isinstance(result, list)
        assert all(isinstance(v, float) for v in result)

    def test_list_returns_list_of_mpf(self):
        result = gwr(lambda s: -mp.log(s) / s - 0.577216 / s, [1.0, 5.0], 4, dtype="mpf")
        assert isinstance(result, list)
        assert all(isinstance(v, mp.mpf) for v in result)

    def test_ndarray_returns_float64(self):
        result = gwr(
            lambda s: -mp.log(s) / s - 0.577216 / s,
            np.array([1.0, 5.0]), 4
        )
        assert isinstance(result, np.ndarray)
        assert result.dtype == np.float64

    def test_ndarray_returns_object_mpf(self):
        result = gwr(
            lambda s: -mp.log(s) / s - 0.577216 / s,
            np.array([1.0, 5.0]), 4, dtype="mpf"
        )
        assert isinstance(result, np.ndarray)
        assert result.dtype == object
