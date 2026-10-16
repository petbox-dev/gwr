"""
Optional python-flint backend for fast arbitrary-precision special functions.

When python-flint is installed, functions here use ARB (via FLINT) instead of
mpmath, providing ~10-60x speedup for Bessel function evaluations.

ARB evaluates in ball arithmetic: every result carries a rigorous error
radius. ARB does not compensate internally for the cancellation in the
K_n(x) power series (K_n(x) ~ e^{-x} while the series terms grow like e^{x}),
so at a fixed working precision the result loses about ``2 x log2(e)`` bits.
The radius reports that loss. Each function here therefore re-evaluates at a
higher working precision until the relative accuracy of the ball reaches the
current mpmath precision.

Usage in a Laplace transform function::

    from gwr_inversion import besselk
    from mpmath import mp, sqrt

    def lap_rad(s):
        b = besselk(0, sqrt(s))   # uses ARB if available, else mp.besselk
        return 1.0 / (s / b + s**2 * 1e3)
"""

from typing import Any, Callable, Optional

from mpmath import mp

# When python-flint is absent, arb/fmpz/_flint_ctx stay unbound; every use is
# behind a HAS_FLINT check.
try:
    from flint import arb, fmpz, ctx as _flint_ctx
    HAS_FLINT = True
except ImportError:
    HAS_FLINT = False


# Bits of relative accuracy required above mp.prec.
_TARGET_GUARD_BITS = 10
# Bits added to the first working precision above the target.
_INITIAL_GUARD_BITS = 16
# Bits added on each retry beyond the measured deficit.
_RETRY_GUARD_BITS = 32
_MAX_ATTEMPTS = 8


def _arb_to_mpf(x: Any) -> Any:
    """Convert the midpoint of an arb ball to mpmath mpf (rounded to mp.prec)."""
    result = mp.mpf.__new__(mp.mpf)
    result._mpf_ = x._mpf_
    return +result


def _mpf_to_arb(x: Any) -> Any:
    """Convert mpmath mpf to an exact arb (zero radius)."""
    sign, mantissa, exponent, _ = x._mpf_
    value = arb(fmpz(int(mantissa))) * arb(2) ** int(exponent)
    return -value if sign else value


def _adaptive(evaluate: Callable[[Any], Any], x: Any) -> Optional[Any]:
    """
    Evaluate ``evaluate(x_arb)`` and raise the ARB working precision until
    the result's relative accuracy is at least ``mp.prec + 10`` bits.

    Returns None when the input is not a positive real number or the loop
    does not converge; the caller then falls back to mpmath.
    """
    try:
        x_mpf = mp.mpf(x)
    except TypeError:
        return None
    if not x_mpf > 0:
        return None

    target_bits = mp.prec + _TARGET_GUARD_BITS
    working_bits = target_bits + _INITIAL_GUARD_BITS
    previous_bits = _flint_ctx.prec
    try:
        # The conversion is exact when the context holds the full mantissa.
        _flint_ctx.prec = max(working_bits, mp.prec + 64)
        x_arb = _mpf_to_arb(x_mpf)
        for _ in range(_MAX_ATTEMPTS):
            _flint_ctx.prec = working_bits
            value = evaluate(x_arb)
            accuracy_bits = value.rel_accuracy_bits()
            if accuracy_bits >= target_bits:
                return _arb_to_mpf(value)
            working_bits += max(target_bits - accuracy_bits, 0) + _RETRY_GUARD_BITS
        return None
    finally:
        _flint_ctx.prec = previous_bits


def _besselk_integer(order: int) -> Callable[[Any], Any]:
    """
    Return an arb evaluator for K_order(x), order >= 0.

    K_0 comes from ARB directly. K_1 comes from the Wronskian
    ``I_0 K_1 + I_1 K_0 = 1/x``: ARB's direct K_1 has a slow path near
    the switch to its asymptotic expansion (up to minutes per call at
    M=128), and the Wronskian form has no cancellation because
    ``I_1 K_0 -> 1/(2x)``. Higher orders use the upward recurrence
    ``K_{n+1} = K_{n-1} + (2n/x) K_n``, which is stable for K.
    """
    def evaluate(x_arb: Any) -> Any:
        k_previous = x_arb.bessel_k(0)
        if order == 0:
            return k_previous
        k_current = (1 / x_arb - x_arb.bessel_i(1) * k_previous) / x_arb.bessel_i(0)
        for n in range(1, order):
            k_previous, k_current = k_current, k_previous + 2 * n / x_arb * k_current
        return k_current
    return evaluate


def besselk(n: int, x: Any) -> Any:
    """
    Modified Bessel function of the second kind K_n(x).

    Uses ARB (python-flint) with adaptive working precision for integer
    ``n`` and real ``x > 0``. Falls back to mpmath otherwise.
    """
    if not HAS_FLINT or int(n) != n:
        return mp.besselk(n, x)
    result = _adaptive(_besselk_integer(abs(int(n))), x)
    return mp.besselk(n, x) if result is None else result


def besseli(n: int, x: Any) -> Any:
    """
    Modified Bessel function of the first kind I_n(x).

    Uses ARB (python-flint) with adaptive working precision for real
    ``x > 0``. Falls back to mpmath otherwise.
    """
    if not HAS_FLINT:
        return mp.besseli(n, x)
    result = _adaptive(lambda x_arb: x_arb.bessel_i(n), x)
    return mp.besseli(n, x) if result is None else result


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
    k0 = _besselk_integer(0)(z_arb)
    k1 = _besselk_integer(1)(z_arb)
    # (pi z / 2) * (2/pi) = z, and (pi z / 2) * (2 z / pi) = z^2.
    return z_arb * (k0 * struve_minus_one_scaled + z_arb * k1 * struve_zero_scaled)


def _besseli0_integral_arb(z_arb: Any) -> Any:
    """``int_0^z I_0(t) dt = z 1F2(1/2; 1, 3/2; z^2/4)``."""
    quarter_z2 = z_arb * z_arb / 4
    return z_arb * quarter_z2.hypgeom([arb(1) / 2], [1, arb(3) / 2])


def besselk0_integral(z: Any) -> Any:
    """
    ``int_0^z K_0(t) dt`` via the modified Struve closed form.

    Uses ARB with adaptive working precision when available, else the same
    closed form in mpmath.
    """
    if mp.mpf(z) == 0:
        return mp.zero
    if HAS_FLINT:
        result = _adaptive(_besselk0_integral_arb, z)
        if result is not None:
            return result
    return mp.pi * z / 2 * (mp.besselk(0, z) * mp.struvel(-1, z)
                            + mp.besselk(1, z) * mp.struvel(0, z))


def besseli0_integral(z: Any) -> Any:
    """
    ``int_0^z I_0(t) dt``.

    Uses ARB with adaptive working precision when available, else mpmath's
    ``besseli(0, z, -1)`` (derivative order -1 is the integral).
    """
    if mp.mpf(z) == 0:
        return mp.zero
    if HAS_FLINT:
        result = _adaptive(_besseli0_integral_arb, z)
        if result is not None:
            return result
    return mp.besseli(0, z, -1)
