"""GWR algorithm implementation and public gwr() entry point."""
from __future__ import annotations

from functools import lru_cache
from inspect import signature
from math import factorial
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple, Union

import numpy as np
from mpmath import mp

from . import _flint
from ._gmpy2 import get_ops, gmpy2_precision

# ---------------------------------------------------------------------------
# Core algorithm (single time point)
# ---------------------------------------------------------------------------


@lru_cache(maxsize=16)
def _gaver_integers(M: int) -> Tuple[List[int], List[List[int]]]:
    """
    Exact integers of the Gaver functional for ``n = 1..M``: the coefficient
    ``(2n)! / (n ((n-1)!)^2) = n C(2n, n)`` and the binomial row ``C(n, i)``.
    Both mpmath mpf and gmpy2 mpfr multiply exactly by a Python int, so the
    sum needs no per-term conversion and each coefficient one rounding.
    """
    # Factorials give n C(2n, n) exactly.
    factorials = [factorial(k) for k in range(2 * M + 1)]
    coefficients = [0] * (M + 1)
    rows: List[List[int]] = [[1] for _ in range(M + 1)]
    for n in range(1, M + 1):
        coefficients[n] = factorials[2 * n] // (n * factorials[n - 1] ** 2)
        # Pascal's rule: one big-int addition per entry instead of a division.
        previous_row = rows[n - 1]
        rows[n] = [1] + [previous_row[i - 1] + previous_row[i] for i in range(1, n)] + [1]
    return coefficients, rows



def cache_clear() -> None:
    """Clear every internal cache: the order-0/1 Bessel values and the Gaver
    integer tables (about 54 MB at M = 1024, per process)."""
    _flint.cache_clear()
    _gaver_integers.cache_clear()


def gwr_single(
    fn: Callable[..., Any],
    time: float,
    M: int,
    precin: int,
    ops: Dict[str, Callable[..., Any]],
) -> Any:
    """
    GWR algorithm for a single time point.

    Parameters
    ----------
    fn : callable
        Laplace transform function (uses mpmath internally).
    time : float
        Time value.
    M : int
        Number of Gaver functional terms.
    precin : int
        Decimal digits of precision.
    ops : dict
        Backend operations dict from _gmpy2.get_ops().
    """
    tau = ops["log2"](precin) / ops["mpf"](time, precin)

    sig = signature(fn).parameters
    n_params = len(sig)

    # fn uses mpmath, so convert tau to mpf for fn calls
    tau_mpf = ops["to_mpf"](tau, precin)

    fni: List[Any] = [0] * (2 * M + 1)
    for i in range(1, 2 * M + 1):
        s_mpf = i * tau_mpf
        if n_params == 1:
            fni[i] = ops["from_mpf"](fn(s_mpf), precin)
        elif n_params == 2:
            fni[i] = ops["from_mpf"](fn(s_mpf, precin), precin)
        else:
            raise TypeError(
                f'Too many arguments for Laplace transform. Expected 1 or 2, got '
                f'{n_params}. Function signature:\n{sig}'
            )

    zero = ops["mpf"](0, precin)
    G0: List[Any] = [zero] * M
    Gp: List[Any] = [zero] * M

    gaver_coefficients, binomial_rows = _gaver_integers(M)
    M1 = M
    for n in range(1, M + 1):
        try:
            row = binomial_rows[n]
            s = zero
            for i in range(n + 1):
                if i % 2 == 1:
                    s = s - row[i] * fni[n + i]
                else:
                    s = s + row[i] * fni[n + i]

            G0[n - 1] = tau * gaver_coefficients[n] * s

        except Exception as e:
            if n == 1:
                raise e
            M1 = n - 1
            break

    best = G0[M1 - 1]
    Gm: List[Any] = [zero] * M1

    broken = False
    for k in range(M1 - 1):
        for n in range(M1 - 2 - k, -1, -1):
            try:
                expr = G0[n + 1] - G0[n]
            except:
                broken = True
                break

            Gp[n] = Gm[n + 1] + (k + 1) / expr
            if k % 2 == 1 and n == M1 - 2 - k:
                best = Gp[n]

        if broken:
            break

        for n in range(M1 - k):
            Gm[n] = G0[n]
            G0[n] = Gp[n]

    return best


# ---------------------------------------------------------------------------
# Core evaluation (single M, no refinement)
# ---------------------------------------------------------------------------


def _gwr_core(
    fn: Union[Callable[[float], Any], Callable[[float, int], Any]],
    time: Union[float, Iterable[float], "np.ndarray[Any, Any]"],
    M: int,
    precin: Optional[int],
    backend: str,
    workers: int,
    dtype: str,
) -> Any:
    """Core GWR evaluation without refinement."""
    # Decimal digits of working precision (mp.dps and the backend context).
    precision_digits = round(21 * M / 10) if precin is None else precin
    ops = get_ops(backend)

    with mp.workdps(precision_digits), gmpy2_precision(precision_digits):
        if not isinstance(time, Iterable):
            value = ops["to_mpf"](gwr_single(fn, time, M, precision_digits, ops),
                                  precision_digits)
            return float(value) if dtype == "float" else value

        is_array = isinstance(time, np.ndarray)
        parallel = workers > 1
        if isinstance(time, np.ndarray) and time.ndim < 1:
            time = [time.item()]
            parallel = False

        # ndarray, list or other Iterable
        times = list(time)
        if parallel:
            from ._parallel import parallel_gwr
            results = parallel_gwr(fn, times, M, precision_digits, backend, workers)
        else:
            results = [gwr_single(fn, t, M, precision_digits, ops) for t in times]
        values = [ops["to_mpf"](result, precision_digits) for result in results]
        if dtype == "float":
            values = [float(value) for value in values]
        if not is_array:
            return values
        return np.array(values, dtype=np.float64 if dtype == "float" else object)


def _materialize_time(
    time: Union[float, Iterable[float], "np.ndarray[Any, Any]"],
) -> Union[float, List[float], "np.ndarray[Any, Any]"]:
    """
    Time points in the form every later step reads: an Iterable becomes a list
    (a generator would otherwise be consumed by the first evaluation and the
    refinement would see nothing), and an ndarray loses its length-1 axes.
    """
    if isinstance(time, np.ndarray):
        if time.ndim >= 2:
            time = np.squeeze(time)
            if time.ndim >= 2:
                raise ValueError(f'Expected ndim < 2, but got {time.ndim}')
        return time
    if isinstance(time, Iterable):
        return list(time)
    return time


# ---------------------------------------------------------------------------
# Adaptive refinement
# ---------------------------------------------------------------------------


def _refine_low_signal(
    fn: Any,
    time: Any,
    result: Any,
    M: int,
    precin: Optional[int],
    backend: str,
    workers: int,
    dtype: str,
    threshold: float,
) -> Any:
    """
    Re-evaluate low-magnitude points at 2*M for improved precision.

    ``time`` is the materialized time points the core evaluation used
    (``_materialize_time``). Non-finite values neither set the peak nor get refined.
    """
    is_ndarray = isinstance(result, np.ndarray)
    magnitudes = np.abs(np.asarray(result, dtype=np.float64))
    finite_magnitudes = magnitudes[np.isfinite(magnitudes)]
    if finite_magnitudes.size == 0:
        return result
    peak = finite_magnitudes.max()
    if peak == 0:
        return result

    mask = magnitudes < peak * threshold
    if not np.any(mask):
        return result

    refine_times = [time_point for time_point, refine_point
                    in zip(np.atleast_1d(np.asarray(time, dtype=object)), mask) if refine_point]
    M_high = 2 * M
    # GWR needs about 2.1 M digits; an explicit precin sized for M is too few for 2 M.
    precin_high = None if precin is None else max(precin, round(21 * M_high / 10))
    refined = _gwr_core(fn, refine_times, M_high, precin_high, backend, workers, dtype)

    # Merge refined values back
    out = result.copy() if is_ndarray else list(result)
    for i, value in zip(np.flatnonzero(mask), refined):
        out[i] = value
    return out


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def gwr(
    fn: Union[Callable[[float], Any], Callable[[float, int], Any]],
    time: Union[float, Iterable[float], "np.ndarray[Any, Any]"],
    M: int = 32,
    precin: Optional[int] = None,
    backend: str = "auto",
    workers: int = 1,
    dtype: str = "float",
    refine: bool = False,
    refine_threshold: float = 1e-4,
) -> Any:
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
        The first parallel call starts a process pool (about 1.5 s on Windows,
        where every worker re-imports the package); later calls with the same
        pool size ``min(workers, os.cpu_count())`` reuse it, and
        ``shutdown_pool()`` releases it. Each worker
        gets one interleaved chunk of the time points. ``fn`` must be
        picklable (a module-level function, not a lambda).
    dtype : str
        Output type. ``'float'`` (default) returns Python float / float64
        ndarray. ``'mpf'`` returns mpmath mpf objects preserving arbitrary
        precision.
    refine : bool
        When True and input is array-like, re-evaluate low-magnitude points
        at ``2 * M`` to improve precision where the solution approaches zero
        (an explicit ``precin`` is raised to ``round(2.1 * 2 * M)`` for them).
        Default False.
    refine_threshold : float
        Fraction of peak magnitude below which points are refined.
        Default ``1e-4``.
    """
    if dtype not in ("float", "mpf"):
        raise ValueError(f"dtype must be 'float' or 'mpf', got {dtype!r}")

    time = _materialize_time(time)
    result = _gwr_core(fn, time, M, precin, backend, workers, dtype)

    if refine and isinstance(result, (np.ndarray, list)):
        result = _refine_low_signal(
            fn, time, result, M, precin, backend, workers, dtype,
            refine_threshold)

    return result


def talbot(
    fn: Callable[[float], Any],
    time: Union[float, Iterable[float], "np.ndarray[Any, Any]"],
    degree: int = 32,
) -> Any:
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
