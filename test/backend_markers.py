"""Skip markers for tests that depend on an optional backend."""

import pytest

from gwr_inversion._gmpy2 import HAS_GMPY2

# A test that needs gmpy2 to run at all.
needs_gmpy2 = pytest.mark.skipif(not HAS_GMPY2, reason="gmpy2 not installed")

# A reference-data test whose stored value is the gmpy2 backend's rounding at low M;
# the mpmath backend rounds differently and fails it
# (docs/todo/open/reference-tests-require-gmpy2.md).
pins_gmpy2_rounding = pytest.mark.skipif(
    not HAS_GMPY2, reason="reference value pins the gmpy2 backend's rounding")
