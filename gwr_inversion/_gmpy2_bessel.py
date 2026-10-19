"""
gmpy2/MPFR backend for Bessel functions K_n(x) and I_n(x), n in {0, 1}.

Uses power series for small x and asymptotic expansion for large x,
all in gmpy2 MPFR arithmetic at the current context precision.

The K_n power series cancels: its terms grow like e^{x} while K_n(x) ~ e^{-x},
so about ``2 x log2(e)`` bits are lost. The series path adds that many guard
bits. The asymptotic expansion has a minimal-term relative error near
``e^{-2x}``, so it is used once ``2 x log2(e)`` exceeds the working precision.

``besselk_gmpy2`` / ``besseli_gmpy2`` return None when they decline (gmpy2
absent, order not 0 or 1, ``x`` outside ``_numeric.positive_real``);
``_bessel`` then falls back to mpmath.
"""
from __future__ import annotations

import math
from typing import Any, Optional, Tuple

from mpmath import mp

from ._gmpy2 import HAS_GMPY2, gm_from_mpf, gm_to_mpf, gmpy2, mpfr
from ._numeric import LOG2_E, positive_real

# Whether the Bessel dispatch uses this backend. Separate from ``_gmpy2.HAS_GMPY2``
# (which selects the GWR arithmetic) so that a benchmark can switch off the
# gmpy2 Bessel functions alone.
HAS_GMPY2_BESSEL = HAS_GMPY2

# Bits beyond the working precision at which the asymptotic expansion is used.
_ASYMPTOTIC_GUARD_BITS = 16
# Bits beyond the cancellation estimate on the series path.
_SERIES_GUARD_BITS = 32


def _uses_asymptotic(x_abs: float, base_bits: int) -> bool:
    """True when the asymptotic expansion reaches ``base_bits`` of precision."""
    return 2.0 * x_abs * LOG2_E > base_bits + _ASYMPTOTIC_GUARD_BITS


def _series_guard_bits(x_abs: float) -> int:
    """Extra bits that compensate the K_n series cancellation at x."""
    return int(2.0 * x_abs * LOG2_E + 2.0 * math.log2(x_abs + 2.0)) + _SERIES_GUARD_BITS


def _series_tolerance() -> Any:
    """Relative size below which a series term no longer changes the sum
    (10 bits below the current gmpy2 context precision)."""
    return gmpy2.exp2(mpfr(-gmpy2.get_context().precision - 10))


def _besseli0_series(x: Any) -> Any:
    """I_0(x) via power series. x is gmpy2 mpfr."""
    tolerance = _series_tolerance()
    quarter_x_squared = x * x / mpfr(4)  # (x/2)^2
    term = mpfr(1)
    series_sum = mpfr(1)
    term_index = 1
    while True:
        term = term * quarter_x_squared / mpfr(term_index * term_index)
        series_sum = series_sum + term
        if abs(term) < abs(series_sum) * tolerance:
            break
        term_index += 1
    return series_sum


def _besseli1_series(x: Any) -> Any:
    """I_1(x) via power series. x is gmpy2 mpfr."""
    tolerance = _series_tolerance()
    quarter_x_squared = x * x / mpfr(4)  # (x/2)^2
    term = x / mpfr(2)
    series_sum = term
    term_index = 1
    while True:
        term = term * quarter_x_squared / mpfr(term_index * (term_index + 1))
        series_sum = series_sum + term
        if abs(term) < abs(series_sum) * tolerance:
            break
        term_index += 1
    return series_sum


def _besselk0_series(x: Any) -> Tuple[Any, Any]:
    """(K_0(x), I_0(x)) via power series; the K_0 series sums I_0 on the way.
    x is gmpy2 mpfr."""
    tolerance = _series_tolerance()
    gamma = gmpy2.const_euler()
    quarter_x_squared = x * x / mpfr(4)  # (x/2)^2

    # I_0 terms and the harmonic-number-weighted sum H_k (x/2)^{2k} / (k!)^2.
    i0_term = mpfr(1)
    i0_sum = mpfr(1)
    harmonic_number = mpfr(0)
    harmonic_sum = mpfr(0)  # sum starts at k=1

    term_index = 1
    while True:
        i0_term = i0_term * quarter_x_squared / mpfr(term_index * term_index)
        i0_sum = i0_sum + i0_term
        harmonic_number = harmonic_number + mpfr(1) / mpfr(term_index)
        harmonic_sum = harmonic_sum + harmonic_number * i0_term
        if abs(i0_term) < abs(i0_sum) * tolerance:
            break
        term_index += 1

    return -(gmpy2.log(x / mpfr(2)) + gamma) * i0_sum + harmonic_sum, i0_sum


def _besselk1_wronskian(x: Any, i0: Any, i1: Any, k0: Any) -> Any:
    """K_1(x) from the Wronskian: I_0 K_1 + I_1 K_0 = 1/x."""
    return (mpfr(1) / x - i1 * k0) / i0


def _asymptotic_sum(n: int, x: Any, alternating: bool) -> Any:
    """``sum_k (+-1)^k a_k(n) / x^k`` with ``a_k(n) = prod_{j=1}^{k} (4n^2 - (2j-1)^2)
    / (k! 8^k)``, stopped at the smallest term. Signs alternate for I_n."""
    tolerance = _series_tolerance()
    inverse_x = mpfr(1) / x
    term = mpfr(1)
    series_sum = mpfr(1)
    four_n_squared = mpfr(4 * n * n)
    term_index = 1
    previous_magnitude = abs(term)
    while True:
        odd = mpfr(2 * term_index - 1)
        term = (term * (four_n_squared - odd * odd) / (mpfr(8) * mpfr(term_index))
                * inverse_x)
        if alternating:
            term = -term
        magnitude = abs(term)
        if magnitude >= previous_magnitude:
            break  # divergent: stop
        series_sum = series_sum + term
        if magnitude < abs(series_sum) * tolerance:
            break
        previous_magnitude = magnitude
        term_index += 1
    return series_sum


def _besselk_asymptotic(n: int, x: Any) -> Any:
    """K_n(x) via asymptotic expansion. x is gmpy2 mpfr."""
    # K_n(x) = sqrt(pi/(2x)) * e^{-x} * sum a_k(n) / x^k
    pi = gmpy2.const_pi()
    prefix = gmpy2.sqrt(pi / (mpfr(2) * x)) * gmpy2.exp(-x)
    return prefix * _asymptotic_sum(n, x, alternating=False)


def _besseli_asymptotic(n: int, x: Any) -> Any:
    """I_n(x) via asymptotic expansion. x is gmpy2 mpfr."""
    # I_n(x) = e^x / sqrt(2*pi*x) * sum (-1)^k * a_k(n) / x^k
    pi = gmpy2.const_pi()
    prefix = gmpy2.exp(x) / gmpy2.sqrt(mpfr(2) * pi * x)
    return prefix * _asymptotic_sum(n, x, alternating=True)


def _accepted_argument(n: Any, x: Any) -> Optional[Any]:
    """``x`` as a positive finite mpf when this backend evaluates order ``n`` there."""
    if not (HAS_GMPY2_BESSEL and HAS_GMPY2) or n not in (0, 1):
        return None
    return positive_real(x)


def besselk_gmpy2(n: int, x: Any) -> Optional[Any]:
    """K_n(x) for n in {0, 1} using gmpy2 MPFR arithmetic; None when declined.

    Output: mpmath mpf at mp.prec. The working precision derives from mp.prec,
    not the gmpy2 context. Uses ``2 x log2(e)`` extra working bits on the series
    path to compensate for cancellation (K_n(x) ~ e^{-x} while intermediate
    terms grow like e^x).
    """
    x_mpf = _accepted_argument(n, x)
    if x_mpf is None:
        return None

    x_abs = float(x_mpf)
    base_bits = mp.prec

    if _uses_asymptotic(x_abs, base_bits):
        with gmpy2.context(gmpy2.get_context(),
                           precision=base_bits + _ASYMPTOTIC_GUARD_BITS):
            result = _besselk_asymptotic(n, gm_from_mpf(x_mpf, 0))
        # e^{-x} underflows MPFR's exponent range near x = 8e8; mpmath has none.
        return None if result == 0 else gm_to_mpf(result, 0)

    # Series path: boost precision to handle cancellation
    with gmpy2.context(gmpy2.get_context(),
                       precision=base_bits + _series_guard_bits(x_abs)):
        x_mpfr = gm_from_mpf(x_mpf, 0)
        k0, i0 = _besselk0_series(x_mpfr)
        if n == 0:
            result = k0
        else:
            result = _besselk1_wronskian(x_mpfr, i0, _besseli1_series(x_mpfr), k0)

    return gm_to_mpf(result, 0)


def besseli_gmpy2(n: int, x: Any) -> Optional[Any]:
    """I_n(x) for n in {0, 1} using gmpy2 MPFR arithmetic; None when declined.

    Output: mpmath mpf at mp.prec. The working precision derives from mp.prec,
    not the gmpy2 context.
    """
    x_mpf = _accepted_argument(n, x)
    if x_mpf is None:
        return None

    base_bits = mp.prec
    use_asymptotic = _uses_asymptotic(float(x_mpf), base_bits)
    with gmpy2.context(gmpy2.get_context(),
                       precision=base_bits + _ASYMPTOTIC_GUARD_BITS):
        x_mpfr = gm_from_mpf(x_mpf, 0)
        if use_asymptotic:
            result = _besseli_asymptotic(n, x_mpfr)
        elif n == 0:
            result = _besseli0_series(x_mpfr)
        else:
            result = _besseli1_series(x_mpfr)

    # e^{x} overflows MPFR's exponent range near x = 8e8; mpmath has none.
    return gm_to_mpf(result, 0) if gmpy2.is_finite(result) else None
