"""
gmpy2/MPFR backend for multiprecision arithmetic.

Provides mpmath ops by default; uses gmpy2.mpfr when available for ~10x speedup.
"""
from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Callable, Dict, Generator, Tuple

from mpmath import mp

try:
    import gmpy2
    from gmpy2 import mpfr
    HAS_GMPY2 = True
except ImportError:
    gmpy2 = None
    mpfr = None
    HAS_GMPY2 = False


def _gmpy2_bits(precin: int) -> int:
    """Convert decimal digits of precision to bits."""
    return int(precin * 3.3219280948873626) + 10


@contextmanager
def gmpy2_precision(precin: int) -> Generator[None, None, None]:
    """Context manager: set gmpy2 context precision for all arithmetic within the block."""
    if not HAS_GMPY2:
        yield
        return
    with gmpy2.context(gmpy2.get_context(), precision=_gmpy2_bits(precin)):
        yield


def get_backend(name: str = "auto") -> str:
    """Return a backend name ('mpmath' or 'gmpy2')."""
    if name == "auto":
        return "gmpy2" if HAS_GMPY2 else "mpmath"
    if name == "gmpy2" and not HAS_GMPY2:
        raise ImportError("gmpy2 is not installed. Install with: pip install gmpy2")
    return name


# --- mpmath backend ops ---

def mp_log2(precin: int) -> Any:
    return mp.log(2)


def mp_mpf(val: Any, precin: int) -> Any:
    return mp.mpf(val)


def mp_fac(n: int, precin: int) -> Any:
    return mp.fac(n)


def mp_binomial(n: int, i: int, precin: int) -> Any:
    return mp.binomial(n, i)


def mp_to_mpf(val: Any, precin: int) -> Any:
    return val


def mp_from_mpf(val: Any, precin: int) -> Any:
    return val


# --- gmpy2 backend ops ---

def gm_log2(precin: int) -> Any:
    return gmpy2.log(mpfr(2))


def gm_mpf(val: Any, precin: int) -> Any:
    return mpfr(str(val))


def gm_fac(n: int, precin: int) -> Any:
    return mpfr(gmpy2.fac(n))


def gm_binomial(n: int, i: int, precin: int) -> Any:
    return mpfr(gmpy2.comb(n, i))


def gm_to_mpf(val: Any, precin: int) -> Any:
    """Convert gmpy2 mpfr to mpmath mpf for public API consistency."""
    dps = mp.dps
    mp.dps = precin
    if gmpy2.is_finite(val):
        mantissa, exponent = val.as_mantissa_exp()
        result = mp.mpf((int(mantissa), int(exponent)))
    else:
        result = mp.mpf(str(val))  # inf and nan carry no mantissa
    mp.dps = dps
    return result


def gm_from_mpf(val: Any, precin: int) -> Any:
    """Convert mpmath mpf to gmpy2 mpfr from the mantissa and exponent.

    Other types (Python int/float, mpc) keep the string conversion.
    """
    if not isinstance(val, mp.mpf):
        return mpfr(str(val))
    sign, mantissa, exponent, _ = val._mpf_
    if mantissa == 0:
        return mpfr(str(val))  # zero, inf and nan carry no mantissa
    value = gmpy2.mul_2exp(mpfr(int(mantissa)), int(exponent))
    return -value if sign else value


# --- Backend dispatch ---

OpsDict = Dict[str, Callable[..., Any]]

BACKENDS: Dict[str, OpsDict] = {
    "mpmath": {
        "log2": mp_log2,
        "mpf": mp_mpf,
        "fac": mp_fac,
        "binomial": mp_binomial,
        "to_mpf": mp_to_mpf,
        "from_mpf": mp_from_mpf,
    },
    "gmpy2": {
        "log2": gm_log2,
        "mpf": gm_mpf,
        "fac": gm_fac,
        "binomial": gm_binomial,
        "to_mpf": gm_to_mpf,
        "from_mpf": gm_from_mpf,
    },
}


def get_ops(backend_name: str) -> Tuple[OpsDict, str]:
    """Return the ops dict and resolved backend name."""
    name = get_backend(backend_name)
    return BACKENDS[name], name
