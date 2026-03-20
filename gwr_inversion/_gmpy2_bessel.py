"""
gmpy2/MPFR backend for Bessel functions K_n(x) and I_n(x), n in {0, 1}.

Uses power series for small x and asymptotic expansion for large x,
all in gmpy2 MPFR arithmetic at the current context precision.
"""

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


def _mpf_to_mpfr(x: Any) -> Any:
    """Convert mpmath mpf to gmpy2 mpfr via string."""
    return mpfr(str(x))


def _mpfr_to_mpf(x: Any) -> Any:
    """Convert gmpy2 mpfr to mpmath mpf via string."""
    dps = mp.dps
    # Use enough digits to round-trip
    mp.dps = int(gmpy2.get_context().precision / 3.32) + 5
    result = mp.mpf(str(x))
    mp.dps = dps
    return result


def _prec_digits() -> int:
    """Return current gmpy2 precision in decimal digits."""
    return int(gmpy2.get_context().precision / 3.3219280948873626)


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

    Input: mpmath mpf. Output: mpmath mpf.
    Uses extra working bits for the series path to compensate for
    cancellation (K_n(x) ~ e^{-x} while intermediate terms grow like e^x).
    """
    if n not in (0, 1):
        return mp.besselk(n, x)

    prec_digits = _prec_digits()
    x_abs = float(abs(x))

    if x_abs > 1.5 * prec_digits:
        xm = _mpf_to_mpfr(x)
        result = _besselk_asymptotic(n, xm)
        return _mpfr_to_mpf(result)

    # Series path: boost precision to handle cancellation
    base_prec = gmpy2.get_context().precision
    extra_bits = int(x_abs * 1.5) + 64
    with gmpy2.context(gmpy2.get_context(), precision=base_prec + extra_bits):
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

    Input: mpmath mpf. Output: mpmath mpf.
    """
    if n not in (0, 1):
        return mp.besseli(n, x)

    xm = _mpf_to_mpfr(x)
    prec_digits = _prec_digits()

    if float(abs(x)) > 1.5 * prec_digits:
        result = _besseli_asymptotic(n, xm)
    elif n == 0:
        result = _besseli0_series(xm)
    else:
        result = _besseli1_series(xm)

    return _mpfr_to_mpf(result)
