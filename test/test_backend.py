"""Tests for backend abstraction and gmpy2 acceleration."""

import pytest
import numpy as np
from mpmath import mp

from gwr_inversion import gwr
from gwr_inversion._gmpy2 import HAS_GMPY2

needs_gmpy2 = pytest.mark.skipif(not HAS_GMPY2, reason="gmpy2 not installed")


class TestMpmathBackend:

    def test_ln_scalar(self):
        result = float(gwr(lambda s: -mp.log(s) / s - 0.577216 / s, 5.0, 32, backend="mpmath"))
        assert abs(result - 1.6094375773356333) < 1e-10

    def test_unit_step_array(self):
        time = np.array([0.5, 2.0, 5.0])
        result = gwr(lambda s: mp.exp(-s) / s, time, 32, backend="mpmath")
        assert isinstance(result, np.ndarray)
        assert len(result) == 3


@needs_gmpy2
class TestGmpy2Backend:

    def test_ln_scalar_matches(self):
        mp_result = float(gwr(
            lambda s: -mp.log(s) / s - 0.577216 / s, 5.0, 32, backend="mpmath"
        ))
        gm_result = float(gwr(
            lambda s: -mp.log(s) / s - 0.577216 / s, 5.0, 32, backend="gmpy2"
        ))
        assert abs(mp_result - gm_result) < 1e-10

    def test_unit_step_array_matches(self):
        time = np.array([0.5, 2.0, 5.0])
        mp_result = np.array([float(v) for v in gwr(
            lambda s: mp.exp(-s) / s, time, 32, backend="mpmath"
        )])
        gm_result = np.array([float(v) for v in gwr(
            lambda s: mp.exp(-s) / s, time, 32, backend="gmpy2"
        )])
        np.testing.assert_allclose(mp_result, gm_result, atol=1e-10)

    def test_auto_selects_gmpy2(self):
        from gwr_inversion._gmpy2 import get_backend
        assert get_backend("auto") == "gmpy2"


class TestBackendErrors:

    def test_gmpy2_not_installed_error(self):
        if HAS_GMPY2:
            pytest.skip("gmpy2 is installed")
        with pytest.raises(ImportError, match="gmpy2 is not installed"):
            gwr(lambda s: 1 / s, 1.0, 4, backend="gmpy2")
