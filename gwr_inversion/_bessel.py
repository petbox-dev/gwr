"""Bessel function dispatch across backends (gmpy2, flint, mpmath)."""

from typing import Any

from mpmath import mp

from ._flint import (
    besselk as besselk_flint, besseli as besseli_flint, HAS_FLINT,
)
from ._gmpy2_bessel import (
    besselk_gmpy2 as _besselk_gmpy2, besseli_gmpy2 as _besseli_gmpy2,
    HAS_GMPY2_BESSEL,
)


def besselk(n: int, x: Any) -> Any:
    """Modified Bessel function of the second kind K_n(x).

    Tries gmpy2 MPFR first, then ARB (python-flint), then mpmath.
    """
    if HAS_GMPY2_BESSEL:
        return _besselk_gmpy2(n, x)
    if HAS_FLINT:
        return besselk_flint(n, x)
    return mp.besselk(n, x)


def besseli(n: int, x: Any) -> Any:
    """Modified Bessel function of the first kind I_n(x).

    Tries ARB (python-flint) first (fastest for I_n which has no exponential
    cancellation), then gmpy2 MPFR, then mpmath.
    """
    if HAS_FLINT:
        return besseli_flint(n, x)
    if HAS_GMPY2_BESSEL:
        return _besseli_gmpy2(n, x)
    return mp.besseli(n, x)
