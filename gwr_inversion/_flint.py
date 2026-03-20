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
    from flint import arb, fmpz, ctx as _flint_ctx
    HAS_FLINT = True
except ImportError:
    arb = None
    fmpz = None
    _flint_ctx = None
    HAS_FLINT = False


def _arb_to_mpf(x: Any) -> Any:
    """Convert arb to mpmath mpf via the internal _mpf_ tuple (zero-copy)."""
    result = mp.mpf.__new__(mp.mpf)
    result._mpf_ = x._mpf_
    return result


def _mpf_to_arb(x: Any) -> Any:
    """Convert mpmath mpf to arb via the _mpf_ tuple (lossless)."""
    sign, man, exp, bc = x._mpf_
    a = arb(fmpz(int(man))) * arb(2) ** int(exp)
    return -a if sign else a


def _base_bits() -> int:
    """Return the base precision in bits from current mp.dps."""
    return int(mp.dps * 3.3219280948873626) + 10


def besselk(n: int, x: Any) -> Any:
    """
    Modified Bessel function of the second kind K_n(x).

    Uses ARB (python-flint) when available with extra working bits to
    compensate for K_n(x) ~ e^{-x} cancellation. Falls back to mpmath
    when the required overshoot would more than triple the base precision
    (crossover where mpmath's asymptotic expansion is faster anyway).
    """
    if not HAS_FLINT:
        return mp.besselk(n, x)
    bits = _base_bits()
    x_abs = float(abs(x))
    extra_bits = int(x_abs * 1.5) + 64
    total_bits = bits + extra_bits
    if total_bits > 3 * bits:
        return mp.besselk(n, x)
    _flint_ctx.prec = total_bits
    x_arb = _mpf_to_arb(x)
    return _arb_to_mpf(x_arb.bessel_k(n))


def besseli(n: int, x: Any) -> Any:
    """
    Modified Bessel function of the first kind I_n(x).

    Uses ARB (python-flint) when available. Falls back to mpmath when
    flint is not installed.
    """
    if not HAS_FLINT:
        return mp.besseli(n, x)
    bits = _base_bits()
    _flint_ctx.prec = bits
    x_arb = _mpf_to_arb(x)
    return _arb_to_mpf(x_arb.bessel_i(n))
