"""Bessel function dispatch across backends (flint, gmpy2, mpmath)."""

from typing import Any, Callable, Optional, Tuple

from mpmath import mp

from ._flint import (
    besselk as besselk_flint, besseli as besseli_flint, HAS_FLINT,
    besselk0_integral as besselk0_integral_flint,
    besseli0_integral as besseli0_integral_flint,
    bessel_k_over_i as bessel_k_over_i_flint,
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


# log2(K_n(z) / I_n(z)) <= _RATIO_BOUND_BITS - 2 z log2(e) for z >= 1, n in {0, 1}.
# Measured maximum of the left side plus 2 z log2(e) over z in [1, 4000]: 1.65
# (n = 0, approaching log2(pi) from below) and 2.98 (n = 1, at z = 1).
_RATIO_BOUND_BITS = 3.5
_TWO_LOG2_E = 2.8853900817779268
_RATIO_PROBE_PRECISION = 64


def guard_sensitivity(log2_sensitivity: Callable[[Any], float]) -> Callable[[Any], float]:
    """Treat a zero denominator in the caller's formula as "no bits can be spared"."""
    def guarded_sensitivity(ratio: Any) -> float:
        try:
            return log2_sensitivity(ratio)
        except ZeroDivisionError:
            return float("inf")
    return guarded_sensitivity


def plan_k_over_i(
    z_float: float, log2_sensitivity: Callable[[Any], float],
) -> Tuple[Optional[Any], int]:
    """
    For ``K_n(z) / I_n(z)`` at ``z = z_float``: return ``(bound, first_bits)``.

    ``bound`` is an upper bound on the ratio when the term is negligible for the
    caller (no evaluation needed), else None. ``first_bits`` is the ARB
    precision for the first evaluation: the bits the bound says are needed,
    plus guard bits, instead of full precision.
    """
    target_bits = mp.prec + 10
    full_bits = target_bits + 16
    if z_float < 1.0:
        return None, full_bits
    ratio_bound = mp.mpf(2) ** (_RATIO_BOUND_BITS - _TWO_LOG2_E * z_float)
    needed_bits = target_bits + log2_sensitivity(ratio_bound)
    if needed_bits <= 0:
        return ratio_bound, 0
    return None, int(min(max(needed_bits, 0.0) + 32, full_bits))


def bessel_k_over_i(order: int, z: Any, log2_sensitivity: Callable[[Any], float]) -> Any:
    """
    ``K_order(z) / I_order(z)`` (order 0 or 1) only as accurate as the result needs.

    In a bounded-reservoir solution the outer boundary enters through this
    ratio at ``z = sqrt(s) r_eD``, where it is about ``pi e^{-2z}``: a small
    correction that needs few bits, and none once ``z`` is large. Evaluating
    K and I there to full precision forces the K-series cancellation retries
    (``2 z log2(e)`` bits) for no gain.

    ``log2_sensitivity(r)`` returns log2 of ``|d ln(result) / d ln(r)|`` for a
    ratio value ``r``; it must be non-decreasing in ``r``. The ratio returned
    is accurate to ``mp.prec + 10 + log2_sensitivity(R)`` bits, and is
    ``mp.zero`` when that is <= 0 (dropping the term moves the result by less
    than ``2^-(mp.prec + 10)``).
    """
    guarded_sensitivity = guard_sensitivity(log2_sensitivity)
    target_bits = mp.prec + 10
    z_mpf = mp.mpf(z)
    bound, first_bits = plan_k_over_i(float(z_mpf), guarded_sensitivity)
    if bound is not None:
        return mp.zero
    if HAS_FLINT and z_mpf > 0:
        ratio = bessel_k_over_i_flint(order, z_mpf, guarded_sensitivity, first_bits)
        if ratio is not None:
            return ratio
    with mp.workprec(_RATIO_PROBE_PRECISION):
        probe = besselk(order, z_mpf) / besseli(order, z_mpf)
    needed_bits = target_bits + guarded_sensitivity(2 * abs(probe))
    if needed_bits <= 0:
        return mp.zero
    with mp.workprec(max(int(needed_bits) + 16, _RATIO_PROBE_PRECISION)):
        ratio = besselk(order, z_mpf) / besseli(order, z_mpf)
    return +ratio


def besselk0_integral(z: Any) -> Any:
    """``int_0^z K_0(t) dt`` (ARB when available, else mpmath closed form)."""
    return besselk0_integral_flint(z)


def besseli0_integral(z: Any) -> Any:
    """``int_0^z I_0(t) dt`` (ARB when available, else mpmath)."""
    return besseli0_integral_flint(z)
