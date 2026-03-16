"""Tests for Fixed Talbot method wrapper."""

import math

import numpy as np
from mpmath import mp

from gwr_inversion import talbot


class TestTalbot:

    def test_ln_scalar(self):
        result = float(talbot(lambda s: -mp.log(s) / s - 0.577216 / s, 5.0, degree=32))
        assert abs(result - math.log(5.0)) < 0.01

    def test_exp_decay(self):
        """L^-1{1/(s+1)} = e^(-t)."""
        result = float(talbot(lambda s: 1 / (s + 1), 1.0, degree=32))
        assert abs(result - math.exp(-1.0)) < 1e-10

    def test_array_input(self):
        time = np.array([1.0, 2.0, 5.0])
        result = talbot(lambda s: 1 / (s + 1), time, degree=32)
        assert isinstance(result, np.ndarray)
        expected = np.exp(-time)
        np.testing.assert_allclose(
            np.array([float(v) for v in result]), expected, atol=1e-10
        )

    def test_list_input(self):
        result = talbot(lambda s: 1 / (s + 1), [1.0, 2.0], degree=32)
        assert isinstance(result, list)
        assert len(result) == 2
