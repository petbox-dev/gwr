"""Bessel function dispatch across backends (flint, gmpy2, mpmath)."""

from typing import Any

from mpmath import mp

from ._flint import (
    besselk as besselk_flint, besseli as besseli_flint, HAS_FLINT,
    besselk0_integral as besselk0_integral_flint,
    besseli0_integral as besseli0_integral_flint,
)
from ._gmpy2_bessel import (
    besselk_gmpy2 as _besselk_gmpy2, besseli_gmpy2 as _besseli_gmpy2,
    HAS_GMPY2_BESSEL,
)


def besselk(n: int, x: Any) -> Any:
    """Modified Bessel function of the second kind K_n(x).

    Tries ARB (python-flint) first, then gmpy2 MPFR, then mpmath. All three
    return a result accurate to the current mp.prec.
    """
    if HAS_FLINT:
        return besselk_flint(n, x)
    if HAS_GMPY2_BESSEL:
        return _besselk_gmpy2(n, x)
    return mp.besselk(n, x)


def besseli(n: int, x: Any) -> Any:
    """Modified Bessel function of the first kind I_n(x).

    Tries ARB (python-flint) first, then gmpy2 MPFR, then mpmath.
    """
    if HAS_FLINT:
        return besseli_flint(n, x)
    if HAS_GMPY2_BESSEL:
        return _besseli_gmpy2(n, x)
    return mp.besseli(n, x)


def besselk0_integral(z: Any) -> Any:
    """``int_0^z K_0(t) dt`` (ARB when available, else mpmath closed form)."""
    return besselk0_integral_flint(z)


def besseli0_integral(z: Any) -> Any:
    """``int_0^z I_0(t) dt`` (ARB when available, else mpmath)."""
    return besseli0_integral_flint(z)
