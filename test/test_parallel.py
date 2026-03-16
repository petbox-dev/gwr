"""Tests for parallel time-point evaluation."""

import numpy as np
from mpmath import mp

from gwr_inversion import gwr


def _lap_ln(s):
    """Ln(t) transform — module-level for pickling."""
    return -mp.log(s) / s - 0.577216 / s


class TestParallel:

    def test_parallel_matches_sequential(self):
        time = np.array([1.0, 5.0, 10.0, 50.0, 100.0])
        seq = np.array([float(v) for v in gwr(_lap_ln, time, 4, workers=1)])
        par = np.array([float(v) for v in gwr(_lap_ln, time, 4, workers=2)])
        np.testing.assert_allclose(seq, par, atol=1e-12)

    def test_parallel_list_input(self):
        time = [1.0, 5.0, 10.0]
        seq = [float(v) for v in gwr(_lap_ln, time, 4, workers=1)]
        par = [float(v) for v in gwr(_lap_ln, time, 4, workers=2)]
        np.testing.assert_allclose(seq, par, atol=1e-12)

    def test_workers_1_is_sequential(self):
        result = gwr(_lap_ln, np.array([1.0, 5.0]), 4, workers=1)
        assert isinstance(result, np.ndarray)
        assert len(result) == 2
