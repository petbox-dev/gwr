"""
GWR Algorithm Numerical Laplace Inversion.

Provides gwr() for Gaver-Wynn-Rho inversion and talbot() for Fixed Talbot.
"""

from typing import Any, Callable, Iterable, Optional, Union

import numpy as np
from mpmath import mp

from ._gmpy2 import get_ops, HAS_GMPY2, gm_to_mpf, gmpy2_precision
from ._flint import besselk, besseli, HAS_FLINT
from ._gwr import gwr_single


def gwr(fn: Union[Callable[[float], Any], Callable[[float, int], Any]],
        time: Union[float, Iterable[float], "np.ndarray[Any, Any]"],
        M: int = 32,
        precin: Optional[int] = None,
        backend: str = "auto",
        workers: int = 1) -> Any:
    """
    Inverse Laplace transform via the Gaver-Wynn-Rho algorithm.

    Parameters
    ----------
    fn : callable
        Laplace transform function using mpmath.
    time : float, list, or ndarray
        Time value(s) at which to evaluate.
    M : int
        Number of Gaver functional terms (default 32).
    precin : int or None
        Decimal digits of precision (default: round(2.1 * M)).
    backend : str
        'auto' (use gmpy2 if available), 'gmpy2', or 'mpmath'.
    workers : int
        Number of parallel workers for array inputs (default 1 = sequential).
    """
    dps = mp.dps
    prec = round(21 * M / 10) if precin is None else precin
    mp.dps = prec

    ops, backend_name = get_ops(backend)

    try:
        with gmpy2_precision(prec):
            if not isinstance(time, Iterable):
                return _to_output(gwr_single(fn, time, M, prec, ops), backend_name, prec)

            if isinstance(time, np.ndarray):
                if time.ndim >= 2:
                    time = np.squeeze(time)
                    if time.ndim >= 2:
                        raise ValueError(f'Expected ndim < 2, but got {time.ndim}')
                if time.ndim < 1:
                    val = gwr_single(fn, time.item(), M, prec, ops)
                    return np.array([_to_output(val, backend_name, prec)], dtype=object)

                times = list(time)
                if workers > 1:
                    from ._parallel import parallel_gwr
                    results = parallel_gwr(fn, times, M, prec, backend, workers)
                else:
                    results = [gwr_single(fn, t, M, prec, ops) for t in times]
                return np.array(
                    [_to_output(r, backend_name, prec) for r in results], dtype=object
                )

            # list or other Iterable
            times = list(time)
            if workers > 1:
                from ._parallel import parallel_gwr
                results = parallel_gwr(fn, times, M, prec, backend, workers)
            else:
                results = [gwr_single(fn, t, M, prec, ops) for t in times]
            return [_to_output(r, backend_name, prec) for r in results]

    finally:
        mp.dps = dps


def _to_output(val: Any, backend_name: str, precin: int) -> Any:
    """Convert backend result to mpmath mpf for public API consistency."""
    if backend_name == "gmpy2":
        return gm_to_mpf(val, precin)
    return val


def talbot(fn: Callable[[float], Any],
           time: Union[float, Iterable[float], "np.ndarray[Any, Any]"],
           degree: int = 32) -> Any:
    """
    Inverse Laplace transform via mpmath's Fixed Talbot method.

    Parameters
    ----------
    fn : callable
        Laplace transform function using mpmath.
    time : float, list, or ndarray
        Time value(s).
    degree : int
        Number of terms (default 32).
    """
    from ._talbot import talbot_single

    if not isinstance(time, Iterable):
        return talbot_single(fn, time, degree)

    if isinstance(time, np.ndarray):
        return np.array([talbot_single(fn, float(t), degree) for t in time], dtype=object)

    return [talbot_single(fn, t, degree) for t in time]


def cache_clear() -> None:
    """Clear all internal caches."""
    pass
