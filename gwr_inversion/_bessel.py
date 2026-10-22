"""
Bessel function dispatch across backends: ARB (python-flint), then gmpy2
MPFR, then mpmath.

Each fast backend returns None for an argument it declines (complex, zero,
negative, non-finite, a precision loop that does not converge, or the
backend is not installed); the next backend takes it, and mpmath takes
everything. All three return a result accurate to the current ``mp.prec``,
for the argument rounded to ``mp.prec`` (the fast backends round an argument
that carries more bits; the transforms never pass one).
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any

from mpmath import mp

from . import _flint, _gmpy2_bessel
from ._numeric import INITIAL_GUARD_BITS, TARGET_GUARD_BITS, TWO_LOG2_E, positive_real


def besselk(n: int, x: Any) -> Any:
    """Modified Bessel function of the second kind K_n(x)."""
    result = _flint.besselk(n, x)
    if result is None:
        result = _gmpy2_bessel.besselk_gmpy2(n, x)
    return mp.besselk(n, x) if result is None else result


def besseli(n: int, x: Any) -> Any:
    """Modified Bessel function of the first kind I_n(x)."""
    result = _flint.besseli(n, x)
    if result is None:
        result = _gmpy2_bessel.besseli_gmpy2(n, x)
    return mp.besseli(n, x) if result is None else result


def besselk0_integral(z: Any) -> Any:
    """``int_0^z K_0(t) dt``: ARB, else the same modified Struve closed form in mpmath."""
    if z == 0:
        return mp.zero
    result = _flint.besselk0_integral(z)
    if result is not None:
        return result
    return mp.pi * z / 2 * (mp.besselk(0, z) * mp.struvel(-1, z)
                            + mp.besselk(1, z) * mp.struvel(0, z))


def besseli0_integral(z: Any) -> Any:
    """``int_0^z I_0(t) dt``: ARB, else mpmath's ``besseli(0, z, -1)`` (derivative
    order -1 is the integral)."""
    if z == 0:
        return mp.zero
    result = _flint.besseli0_integral(z)
    return mp.besseli(0, z, -1) if result is None else result


# log2(K_n(z) / I_n(z)) <= _RATIO_BOUND_BITS - 2 z log2(e) for z >= 1, n in {0, 1}.
# Measured maximum of the left side plus 2 z log2(e) over z in [1, 4000]: 1.65
# (n = 0, approaching log2(pi) from below) and 2.98 (n = 1, at z = 1).
_RATIO_BOUND_BITS = 3.5
# Bits added to the bound's estimate of the needed precision for the first
# evaluation, and to the probe's estimate on the mpmath path.
_PLANNED_GUARD_BITS = 32
_PROBED_GUARD_BITS = 16
# mpmath precision (bits) of the probe that measures the ratio's magnitude.
_RATIO_PROBE_PRECISION = 64


def _guard_sensitivity(log2_sensitivity: Callable[[Any], float]) -> Callable[[Any], float]:
    """
    Treat a zero denominator (or an infinite value) in the caller's formula as
    "no bits can be spared": ``mp.prec`` extra bits, a finite number so that
    the precision arithmetic downstream stays integral.
    """
    def guarded_sensitivity(ratio: Any) -> float:
        try:
            sensitivity = log2_sensitivity(ratio)
        except (ZeroDivisionError, OverflowError):  # zero denominator; exponent beyond float
            return float(mp.prec)
        return sensitivity if math.isfinite(sensitivity) else float(mp.prec)
    return guarded_sensitivity


def _plan_first_bits(z_float: float, log2_sensitivity: Callable[[Any], float]) -> int | None:
    """
    ARB precision for the first evaluation of ``K_n(z) / I_n(z)`` at
    ``z = z_float``: the bits an upper bound on the ratio says are needed, plus
    guard bits, instead of full precision. None when the bound already shows
    the term is negligible for the caller (no evaluation needed).
    """
    target_bits = mp.prec + TARGET_GUARD_BITS
    full_bits = target_bits + INITIAL_GUARD_BITS
    if z_float < 1.0:
        return int(full_bits)
    ratio_bound = mp.mpf(2) ** (_RATIO_BOUND_BITS - TWO_LOG2_E * z_float)
    needed_bits = target_bits + log2_sensitivity(ratio_bound)
    if needed_bits <= 0:
        return None
    return int(min(needed_bits + _PLANNED_GUARD_BITS, full_bits))


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
    than ``2^-(mp.prec + 10)``). A ``z`` that is not a positive finite real
    (complex on the Talbot contour) gets the full-precision ratio.
    """
    if order not in (0, 1):
        raise ValueError(f"bessel_k_over_i supports order 0 or 1, got {order!r}")
    z_mpf = positive_real(z)
    if z_mpf is None:
        return besselk(order, z) / besseli(order, z)
    guarded_sensitivity = _guard_sensitivity(log2_sensitivity)
    target_bits = mp.prec + TARGET_GUARD_BITS
    first_bits = _plan_first_bits(float(z_mpf), guarded_sensitivity)
    if first_bits is None:
        return mp.zero
    ratio = _flint.bessel_k_over_i(order, z_mpf, guarded_sensitivity, first_bits)
    if ratio is not None:
        return ratio
    with mp.workprec(_RATIO_PROBE_PRECISION):
        probe = besselk(order, z_mpf) / besseli(order, z_mpf)
    needed_bits = target_bits + guarded_sensitivity(2 * abs(probe))
    if needed_bits <= 0:
        return mp.zero
    with mp.workprec(max(int(needed_bits) + _PROBED_GUARD_BITS, _RATIO_PROBE_PRECISION)):
        ratio = besselk(order, z_mpf) / besseli(order, z_mpf)
    return +ratio
