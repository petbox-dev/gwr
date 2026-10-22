"""
ARB (python-flint) evaluation of Bessel functions and their integrals.

ARB evaluates in ball arithmetic: every result carries a rigorous error
radius. ARB does not compensate internally for the cancellation in the
K_n(x) power series (K_n(x) ~ e^{-x} while the series terms grow like e^{x}),
so at a fixed working precision the result loses about ``2 x log2(e)`` bits.
The radius reports that loss. Each function here therefore re-evaluates at a
higher working precision until the relative accuracy of the ball reaches
``mp.prec + TARGET_GUARD_BITS`` bits.

Every public function returns None when ARB declines: python-flint is not
installed, the argument is outside ``_numeric.positive_real``, or the
precision loop does not converge. ``_bessel`` owns the fallback to gmpy2 and
mpmath; users call ``gwr_inversion.besselk`` / ``besseli``.
"""
from __future__ import annotations

from collections.abc import Callable
from functools import lru_cache
from typing import Any, NamedTuple

from mpmath import mp

from ._numeric import (
    INITIAL_GUARD_BITS,
    MAX_ATTEMPTS,
    RETRY_GUARD_BITS,
    TARGET_GUARD_BITS,
    positive_real,
)

# When python-flint is absent, arb/fmpz/flint_context stay unbound; every use is
# behind a HAS_FLINT check.
try:
    import flint
    from flint import arb, fmpz
    flint_context = flint.ctx
    HAS_FLINT = True
except ImportError:
    HAS_FLINT = False


# Distinct (argument, precision) pairs kept by the order-0/1 cache. A Laplace
# transform evaluates a handful of Bessel arguments per s.
_ORDER_01_CACHE_SIZE = 32

# A deficit larger than this is treated as "no finite bound" (rel_accuracy_bits
# returns about -2^63 for a ball containing zero or infinity).
_UNBOUNDED_DEFICIT_BITS = 1 << 20

# Bits below which a ball midpoint is not trusted for the ratio's magnitude.
_RATIO_MAGNITUDE_BITS = 4


class BesselOrder01(NamedTuple):
    """K_0, K_1, I_0 and I_1 at one argument (arb balls, or mpf once converted)."""

    k0: Any
    k1: Any
    i0: Any
    i1: Any


def next_working_bits(working_bits: int, target_bits: int, accuracy_bits: int) -> int:
    """Working precision for the next attempt after reaching ``accuracy_bits``."""
    deficit_bits = target_bits - accuracy_bits
    if deficit_bits > _UNBOUNDED_DEFICIT_BITS:
        return 2 * working_bits
    return working_bits + max(deficit_bits, 0) + RETRY_GUARD_BITS


def mpf_from_raw(raw: tuple[Any, ...]) -> Any:
    """The mpf whose ``_mpf_`` tuple is ``raw``, exactly (no rounding)."""
    result = mp.mpf.__new__(mp.mpf)
    result._mpf_ = raw
    return result


def arb_to_mpf(x: Any) -> Any:
    """Convert the midpoint of an arb ball to mpmath mpf (rounded to mp.prec)."""
    return +mpf_from_raw(x._mpf_)


def mpf_to_arb(x: Any) -> Any:
    """Convert mpmath mpf to an exact arb (zero radius) at any context precision."""
    sign, mantissa, exponent, bit_count = x._mpf_
    # arb(fmpz) and the power of 2 are exact, so the product is exact once the
    # context holds the mantissa's bit count.
    with flint_context.workprec(max(bit_count, flint_context.prec)):
        value = arb(fmpz(-int(mantissa) if sign else int(mantissa))) * arb(2) ** int(exponent)
    return value


def _adaptive_balls(
    evaluate: Callable[[Any], tuple[Any, ...]], x_mpf: Any,
) -> tuple[Any, ...] | None:
    """
    Evaluate ``evaluate(x_arb)`` (a tuple of balls) and raise the ARB working
    precision until every ball's relative accuracy is at least
    ``mp.prec + TARGET_GUARD_BITS``; the least accurate ball sets the retry.

    ``x_mpf`` is a positive finite mpf (``positive_real``). Returns the balls'
    midpoints as mpf, or None when the loop does not converge.
    """
    target_bits = mp.prec + TARGET_GUARD_BITS
    working_bits = target_bits + INITIAL_GUARD_BITS
    x_arb = mpf_to_arb(x_mpf)
    for _ in range(MAX_ATTEMPTS):
        with flint_context.workprec(working_bits):
            balls = evaluate(x_arb)
        accuracy_bits = min(ball.rel_accuracy_bits() for ball in balls)
        if accuracy_bits >= target_bits:
            return tuple(arb_to_mpf(ball) for ball in balls)
        working_bits = next_working_bits(working_bits, target_bits, accuracy_bits)
    return None


def _adaptive_ball(evaluate: Callable[[Any], Any], x_mpf: Any) -> Any | None:
    """``_adaptive_balls`` for an ``evaluate`` that returns one ball."""
    values = _adaptive_balls(lambda x_arb: (evaluate(x_arb),), x_mpf)
    return None if values is None else values[0]


def bessel_order_01_arb(x_arb: Any) -> BesselOrder01:
    """
    (K_0, K_1, I_0, I_1) balls at one argument.

    K_0 comes from ARB directly. K_1 comes from the Wronskian
    ``I_0 K_1 + I_1 K_0 = 1/x``: ARB's direct K_1 has a slow path near
    the switch to its asymptotic expansion (up to minutes per call at
    M=128), and the Wronskian form has no cancellation because
    ``I_1 K_0 -> 1/(2x)``.
    """
    k0 = x_arb.bessel_k(0)
    i0 = x_arb.bessel_i(0)
    i1 = x_arb.bessel_i(1)
    k1 = (1 / x_arb - i1 * k0) / i0
    return BesselOrder01(k0, k1, i0, i1)


def _besselk_integer(order: int) -> Callable[[Any], Any]:
    """
    Return an arb evaluator for K_order(x), order >= 2: the upward recurrence
    ``K_{n+1} = K_{n-1} + (2n/x) K_n`` (stable for K) from K_0 and K_1.
    """
    def evaluate(x_arb: Any) -> Any:
        k_previous, k_current, _, _ = bessel_order_01_arb(x_arb)
        for n in range(1, order):
            k_previous, k_current = k_current, k_previous + 2 * n / x_arb * k_current
        return k_current
    return evaluate


@lru_cache(maxsize=_ORDER_01_CACHE_SIZE)
def _bessel_order_01_cached(x_key: tuple[Any, ...], precision_bits: int) -> Any | None:
    """
    (K_0, K_1, I_0, I_1) for the mpf whose ``_mpf_`` tuple is ``x_key``.

    A transform typically needs several of these four at the same argument
    (``radial_bounded_nf`` asks for K_0, K_1, I_0, I_1 at sqrt(s)); K_1 alone
    already needs the other three. One ARB evaluation serves all of them.
    ``precision_bits`` is ``mp.prec`` at the call, so a change of precision
    never returns a stale entry.
    """
    values = _adaptive_balls(bessel_order_01_arb, mpf_from_raw(x_key))
    return None if values is None else BesselOrder01(*values)


def _bessel_order_01(x_mpf: Any) -> BesselOrder01 | None:
    """Cached (K_0, K_1, I_0, I_1) at a positive finite mpf; None if not converged."""
    return _bessel_order_01_cached(x_mpf._mpf_, mp.prec)


def besselk(n: Any, x: Any) -> Any | None:
    """
    K_n(x) in ARB for integer ``n`` and positive finite real ``x``, else None.

    Orders 0 and 1 share one cached evaluation of (K_0, K_1, I_0, I_1) per
    argument.
    """
    if not HAS_FLINT:
        return None
    try:
        if int(n) != n:
            return None
    except (TypeError, ValueError, OverflowError):  # nan, inf, complex or non-numeric order
        return None
    x_mpf = positive_real(x)
    if x_mpf is None:
        return None
    order = abs(int(n))
    if order <= 1:
        values = _bessel_order_01(x_mpf)
        if values is None:
            return None
        return values.k0 if order == 0 else values.k1
    return _adaptive_ball(_besselk_integer(order), x_mpf)


def besseli(n: Any, x: Any) -> Any | None:
    """
    I_n(x) in ARB for positive finite real ``x``, else None. Orders 0 and 1
    share the cached evaluation used by ``besselk``.
    """
    if not HAS_FLINT:
        return None
    try:
        float(n)  # ARB takes a real order (nan and inf end in a non-converged loop: None)
    except (TypeError, ValueError):  # complex or non-numeric order: mpmath decides
        return None
    x_mpf = positive_real(x)
    if x_mpf is None:
        return None
    if n in (0, 1):
        values = _bessel_order_01(x_mpf)
        if values is None:
            return None
        return values.i0 if n == 0 else values.i1
    return _adaptive_ball(lambda x_arb: x_arb.bessel_i(n), x_mpf)


def cache_clear() -> None:
    """Clear the cached order-0/1 Bessel values."""
    if HAS_FLINT:
        _bessel_order_01_cached.cache_clear()


def _bessel_k_over_i_arb(z_arb: Any, order: int) -> Any:
    """K_order(z) / I_order(z) for order 0 or 1 (K_1 from the Wronskian)."""
    if order == 0:
        return z_arb.bessel_k(0) / z_arb.bessel_i(0)
    _, k1, _, i1 = bessel_order_01_arb(z_arb)
    return k1 / i1


def bessel_k_over_i_ball(
    order: int, z_arb: Any, log2_sensitivity: Callable[[Any], float], first_bits: int,
) -> tuple[Any, float] | None:
    """
    Ball for ``K_order(z) / I_order(z)`` and the bits the caller needs from it.

    ``log2_sensitivity(r)`` is log2 of ``|d ln(result) / d ln(r)|`` when the
    ratio has value ``r``, non-decreasing in ``r``; the caller needs
    ``mp.prec + 10 + log2_sensitivity(R)`` bits. The first evaluation runs at
    ``first_bits``; its midpoint gives the magnitude of R. Returns
    ``(ball, needed_bits)`` once the ball's relative accuracy covers
    ``needed_bits`` (or ``needed_bits <= 0``: the term is negligible), and
    None when the loop does not converge (for example when the radius of
    ``z_arb`` itself limits the accuracy).
    """
    target_bits = mp.prec + TARGET_GUARD_BITS
    working_bits = first_bits
    for _ in range(MAX_ATTEMPTS + 2):
        with flint_context.workprec(working_bits):
            ratio = _bessel_k_over_i_arb(z_arb, order)
        accuracy_bits = ratio.rel_accuracy_bits()
        if accuracy_bits < _RATIO_MAGNITUDE_BITS:
            working_bits = max(next_working_bits(working_bits, target_bits, accuracy_bits),
                               target_bits + INITIAL_GUARD_BITS)
            continue
        # 2 |midpoint| bounds |ratio| from above once 4+ bits are correct.
        needed_bits = target_bits + log2_sensitivity(2 * abs(arb_to_mpf(ratio)))
        if needed_bits <= 0 or accuracy_bits >= needed_bits:
            return ratio, needed_bits
        working_bits += int(needed_bits - accuracy_bits) + RETRY_GUARD_BITS
    return None


def bessel_k_over_i(
    order: int, z_mpf: Any, log2_sensitivity: Callable[[Any], float], first_bits: int,
) -> Any | None:
    """
    ``K_order(z) / I_order(z)`` at a positive finite mpf, accurate to the bits
    the caller needs (see ``bessel_k_over_i_ball``); ``mp.zero`` when the term
    is negligible. Returns None when ARB is absent or the loop does not converge.
    """
    if not HAS_FLINT:
        return None
    found = bessel_k_over_i_ball(order, mpf_to_arb(z_mpf), log2_sensitivity, first_bits)
    if found is None:
        return None
    ratio, needed_bits = found
    return mp.zero if needed_bits <= 0 else arb_to_mpf(ratio)


def _besselk0_integral_arb(z_arb: Any) -> Any:
    """
    ``int_0^z K_0(t) dt = (pi z / 2) [K_0(z) L_{-1}(z) + K_1(z) L_0(z)]``.

    L_nu are modified Struve functions, written as 1F2 series:
    ``L_{-1}(z) = (2/pi) 1F2(1; 1/2, 3/2; z^2/4)`` and
    ``L_0(z) = (2z/pi) 1F2(1; 3/2, 3/2; z^2/4)``. Both sums are positive,
    so the combination has no cancellation.
    """
    quarter_z2 = z_arb * z_arb / 4
    half = arb(1) / 2
    three_halves = arb(3) / 2
    struve_minus_one_scaled = quarter_z2.hypgeom([1], [half, three_halves])
    struve_zero_scaled = quarter_z2.hypgeom([1], [three_halves, three_halves])
    k0, k1, _, _ = bessel_order_01_arb(z_arb)
    # (pi z / 2) * (2/pi) = z, and (pi z / 2) * (2 z / pi) = z^2.
    return z_arb * (k0 * struve_minus_one_scaled + z_arb * k1 * struve_zero_scaled)


def _besseli0_integral_arb(z_arb: Any) -> Any:
    """``int_0^z I_0(t) dt = z 1F2(1/2; 1, 3/2; z^2/4)``."""
    quarter_z2 = z_arb * z_arb / 4
    return z_arb * quarter_z2.hypgeom([arb(1) / 2], [1, arb(3) / 2])


def besselk0_integral(z: Any) -> Any | None:
    """``int_0^z K_0(t) dt`` in ARB (modified Struve closed form) for positive
    finite real ``z``, else None."""
    if not HAS_FLINT:
        return None
    z_mpf = positive_real(z)
    return None if z_mpf is None else _adaptive_ball(_besselk0_integral_arb, z_mpf)


def besseli0_integral(z: Any) -> Any | None:
    """``int_0^z I_0(t) dt`` in ARB for positive finite real ``z``, else None."""
    if not HAS_FLINT:
        return None
    z_mpf = positive_real(z)
    return None if z_mpf is None else _adaptive_ball(_besseli0_integral_arb, z_mpf)
