"""
Optional python-flint backend for fast arbitrary-precision special functions.

When python-flint is installed, functions here use ARB (via FLINT) instead of
mpmath, providing ~15x speedup for Bessel function evaluations.

Usage in a Laplace transform function::

    from gwr_inversion import besselk
    from mpmath import mp, sqrt

    def lap_rad(s):
        b = besselk(0, sqrt(s))   # uses ARB if available, else mp.besselk
        return 1.0 / (s / b + s**2 * 1e3)
"""

from typing import Any

from mpmath import mp

try:
    from flint import arb, ctx as _flint_ctx
    HAS_FLINT = True
except ImportError:
    arb = None
    _flint_ctx = None
    HAS_FLINT = False


def _arb_to_mpf(x: Any) -> Any:
    """Convert arb to mpmath mpf via the internal _mpf_ tuple (zero-copy)."""
    result = mp.mpf.__new__(mp.mpf)
    result._mpf_ = x._mpf_
    return result


def _mpf_to_arb(x: Any) -> Any:
    """Convert mpmath mpf to arb via string representation."""
    return arb(mp.nstr(x, int(mp.dps * 1.05) + 5, strip_zeros=False))


def _set_flint_prec() -> None:
    """Sync flint context precision with current mp.dps."""
    _flint_ctx.prec = int(mp.dps * 3.3219280948873626) + 10


def besselk(n: int, x: Any) -> Any:
    """
    Modified Bessel function of the second kind K_n(x).

    Uses ARB (python-flint) when available (~15x faster than mpmath).
    Falls back to mp.besselk otherwise.
    """
    if not HAS_FLINT:
        return mp.besselk(n, x)
    _set_flint_prec()
    return _arb_to_mpf(_mpf_to_arb(x).bessel_k(n))


def besseli(n: int, x: Any) -> Any:
    """
    Modified Bessel function of the first kind I_n(x).

    Uses ARB (python-flint) when available (~15x faster than mpmath).
    Falls back to mp.besseli otherwise.
    """
    if not HAS_FLINT:
        return mp.besseli(n, x)
    _set_flint_prec()
    return _arb_to_mpf(_mpf_to_arb(x).bessel_i(n))
