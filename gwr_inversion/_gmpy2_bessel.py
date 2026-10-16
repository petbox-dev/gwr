"""
gmpy2/MPFR backend for Bessel functions K_n(x) and I_n(x), n in {0, 1}.

Uses power series for small x and asymptotic expansion for large x,
all in gmpy2 MPFR arithmetic at the current context precision.

The K_n power series cancels: its terms grow like e^{x} while K_n(x) ~ e^{-x},
so about ``2 x log2(e)`` bits are lost. The series path adds that many guard
bits. The asymptotic expansion has a minimal-term relative error near
``e^{-2x}``, so it is used once ``2 x log2(e)`` exceeds the working precision.
"""

import math
from typing import Any

from mpmath import mp

try:
    import gmpy2
    from gmpy2 import mpfr
    HAS_GMPY2_BESSEL = True
except ImportError:
    gmpy2 = None
    mpfr = None
    HAS_GMPY2_BESSEL = False


_LOG2_E = 1.4426950408889634
# Bits beyond the working precision at which the asymptotic expansion is used.
_ASYMPTOTIC_GUARD_BITS = 16
# Bits beyond the cancellation estimate on the series path.
_SERIES_GUARD_BITS = 32


def _mpf_to_mpfr(x: Any) -> Any:
    """Convert mpmath mpf to gmpy2 mpfr from the mantissa and exponent (exact
    when the context precision holds the mantissa)."""
    sign, mantissa, exponent, _ = mp.mpf(x)._mpf_
    value = gmpy2.mul_2exp(mpfr(int(mantissa)), int(exponent))
    return -value if sign else value


def _mpfr_to_mpf(x: Any) -> Any:
    """Convert gmpy2 mpfr to mpmath mpf from the mantissa and exponent
    (rounded to mp.prec)."""
    mantissa, exponent = x.as_mantissa_exp()
    return mp.mpf((int(mantissa), int(exponent)))


def _uses_asymptotic(x_abs: float, base_bits: int) -> bool:
    """True when the asymptotic expansion reaches ``base_bits`` of precision."""
    return 2.0 * x_abs * _LOG2_E > base_bits + _ASYMPTOTIC_GUARD_BITS


def _series_guard_bits(x_abs: float) -> int:
    """Extra bits that compensate the K_n series cancellation at x."""
    return int(2.0 * x_abs * _LOG2_E + 2.0 * math.log2(x_abs + 2.0)) + _SERIES_GUARD_BITS


def _besseli0_series(x: Any) -> Any:
    """I_0(x) via power series. x is gmpy2 mpfr."""
    x2 = x * x / mpfr(4)  # (x/2)^2
    term = mpfr(1)
    s = mpfr(1)
    k = 1
    while True:
        term = term * x2 / mpfr(k * k)
        s = s + term
        if abs(term) < abs(s) * gmpy2.exp2(mpfr(-gmpy2.get_context().precision - 10)):
            break
        k += 1
    return s


def _besseli1_series(x: Any) -> Any:
    """I_1(x) via power series. x is gmpy2 mpfr."""
    x2 = x * x / mpfr(4)  # (x/2)^2
    term = x / mpfr(2)
    s = term
    k = 1
    while True:
        term = term * x2 / mpfr(k * (k + 1))
        s = s + term
        if abs(term) < abs(s) * gmpy2.exp2(mpfr(-gmpy2.get_context().precision - 10)):
            break
        k += 1
    return s


def _besselk0_series(x: Any) -> Any:
    """K_0(x) via power series. x is gmpy2 mpfr."""
    gamma = gmpy2.const_euler()
    x2 = x * x / mpfr(4)  # (x/2)^2

    # I_0 and H_k sum simultaneously
    i0_term = mpfr(1)
    i0_sum = mpfr(1)
    hk = mpfr(0)
    hk_sum = mpfr(0)  # sum starts at k=1

    k = 1
    while True:
        i0_term = i0_term * x2 / mpfr(k * k)
        i0_sum = i0_sum + i0_term
        hk = hk + mpfr(1) / mpfr(k)
        hk_sum = hk_sum + hk * i0_term
        if abs(i0_term) < abs(i0_sum) * gmpy2.exp2(mpfr(-gmpy2.get_context().precision - 10)):
            break
        k += 1

    return -(gmpy2.log(x / mpfr(2)) + gamma) * i0_sum + hk_sum


def _besselk1_wronskian(x: Any, i0: Any, i1: Any, k0: Any) -> Any:
    """K_1(x) from the Wronskian: I_0 K_1 + I_1 K_0 = 1/x."""
    return (mpfr(1) / x - i1 * k0) / i0


def _asymptotic_coeffs(n: int, k_max: int) -> list[Any]:
    """Compute asymptotic coefficients a_k(n) for Bessel functions.

    a_0 = 1, a_k(n) = prod_{j=1}^{k} (4n^2 - (2j-1)^2) / (k! * 8^k)
    """
    coeffs = [mpfr(1)]
    four_n2 = mpfr(4 * n * n)
    for k in range(1, k_max + 1):
        j = mpfr(2 * k - 1)
        coeffs.append(coeffs[-1] * (four_n2 - j * j) / (mpfr(8) * mpfr(k)))
    return coeffs


def _besselk_asymptotic(n: int, x: Any) -> Any:
    """K_n(x) via asymptotic expansion. x is gmpy2 mpfr."""
    # K_n(x) = sqrt(pi/(2x)) * e^{-x} * sum a_k(n) / x^k
    pi = gmpy2.const_pi()
    prefix = gmpy2.sqrt(pi / (mpfr(2) * x)) * gmpy2.exp(-x)

    inv_x = mpfr(1) / x
    term = mpfr(1)
    s = mpfr(1)
    four_n2 = mpfr(4 * n * n)
    k = 1
    prev_abs = abs(term)
    while True:
        j = mpfr(2 * k - 1)
        term = term * (four_n2 - j * j) / (mpfr(8) * mpfr(k)) * inv_x
        t_abs = abs(term)
        if t_abs >= prev_abs:
            break  # divergent: stop
        s = s + term
        if t_abs < abs(s) * gmpy2.exp2(mpfr(-gmpy2.get_context().precision - 10)):
            break
        prev_abs = t_abs
        k += 1

    return prefix * s


def _besseli_asymptotic(n: int, x: Any) -> Any:
    """I_n(x) via asymptotic expansion. x is gmpy2 mpfr."""
    # I_n(x) = e^x / sqrt(2*pi*x) * sum (-1)^k * a_k(n) / x^k
    pi = gmpy2.const_pi()
    prefix = gmpy2.exp(x) / gmpy2.sqrt(mpfr(2) * pi * x)

    inv_x = mpfr(1) / x
    term = mpfr(1)
    s = mpfr(1)
    four_n2 = mpfr(4 * n * n)
    k = 1
    prev_abs = abs(term)
    while True:
        j = mpfr(2 * k - 1)
        term = term * (four_n2 - j * j) / (mpfr(8) * mpfr(k)) * inv_x
        # (-1)^k factor: negate each step
        term = -term
        t_abs = abs(term)
        if t_abs >= prev_abs:
            break  # divergent: stop
        s = s + term
        if t_abs < abs(s) * gmpy2.exp2(mpfr(-gmpy2.get_context().precision - 10)):
            break
        prev_abs = t_abs
        k += 1

    return prefix * s


def besselk_gmpy2(n: int, x: Any) -> Any:
    """K_n(x) for n in {0, 1} using gmpy2 MPFR arithmetic.

    Input: mpmath mpf. Output: mpmath mpf at mp.prec.
    The working precision derives from mp.prec, not the gmpy2 context.
    Uses ``2 x log2(e)`` extra working bits on the series path to compensate
    for cancellation (K_n(x) ~ e^{-x} while intermediate terms grow like e^x).
    """
    if n not in (0, 1) or not mp.mpf(x) > 0:
        return mp.besselk(n, x)

    x_abs = float(x)
    base_bits = mp.prec

    if _uses_asymptotic(x_abs, base_bits):
        with gmpy2.context(gmpy2.get_context(),
                           precision=base_bits + _ASYMPTOTIC_GUARD_BITS):
            result = _besselk_asymptotic(n, _mpf_to_mpfr(x))
        return _mpfr_to_mpf(result)

    # Series path: boost precision to handle cancellation
    with gmpy2.context(gmpy2.get_context(),
                       precision=base_bits + _series_guard_bits(x_abs)):
        xm = _mpf_to_mpfr(x)
        if n == 0:
            result = _besselk0_series(xm)
        else:
            i0 = _besseli0_series(xm)
            i1 = _besseli1_series(xm)
            k0 = _besselk0_series(xm)
            result = _besselk1_wronskian(xm, i0, i1, k0)

    return _mpfr_to_mpf(result)


def besseli_gmpy2(n: int, x: Any) -> Any:
    """I_n(x) for n in {0, 1} using gmpy2 MPFR arithmetic.

    Input: mpmath mpf. Output: mpmath mpf at mp.prec.
    The working precision derives from mp.prec, not the gmpy2 context.
    """
    if n not in (0, 1) or not mp.mpf(x) > 0:
        return mp.besseli(n, x)

    base_bits = mp.prec
    use_asymptotic = _uses_asymptotic(float(x), base_bits)
    with gmpy2.context(gmpy2.get_context(),
                       precision=base_bits + _ASYMPTOTIC_GUARD_BITS):
        xm = _mpf_to_mpfr(x)
        if use_asymptotic:
            result = _besseli_asymptotic(n, xm)
        elif n == 0:
            result = _besseli0_series(xm)
        else:
            result = _besseli1_series(xm)

    return _mpfr_to_mpf(result)
