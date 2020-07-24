"""
``GWR(fn , time, M, precin)`` gives the inverse of the Laplace transform function
named ``fn`` for a given array of ``time``. The method involves the calculation
of ``M`` terms of the Gaver-functional. The obtained series is accelerated using
the Wynn-Rho convergence acceleration scheme. The precision of internal
calculations is set to ``precin``.

``GWR(F, t, M)`` does the same, but the precision of the internal calculations is
selected automatically:  ``precin`` = 2.1 M).

GWR(F, t) uses ``M`` = 32 terms and ``precin`` = 67 as defaults. It should give
reasonable results for many problems.

Important note: The Laplace transform should be defined as a function of one
argument. It can involve anything from a simple Mathematica expression to a
sophisticated Module. Since the Laplace transform will be evaluated with
non-standard (multiple) precision, approximate numbers (with decimal point) or
Mathematica functions starting with the letter ``N`` are not allowed in the
function definition.

Example usage:

def fun(s):
    (1/s) * mp.exp(-mp.sqrt( (s ** 2 + (37 / 100) s + 1) / (s ** 2 + s + Pi)))

t0 = 100
GWR(fun, t0)
"""

from inspect import signature
from functools import lru_cache

import numpy as np
from mpmath import mp  # type: ignore

from typing import List, Dict, Tuple, Any, Callable, Union, Iterable, Optional
from typing import cast


MACHINE_PRECISION = 15


LOG2 = mp.log(2.0)


def gwr(fn: Union[Callable[[float], Any], Callable[[float, int], Any]],
        time: Union[float, Iterable[float], np.ndarray],
        M: int = 32, precin: Optional[int] = None) -> Any:
    """
    Gives the inverse of the Laplace transform function ``fn`` for a given array
    of ``time``. The method involves the calculation of ``M`` terms of the
    Gaver functional. The obtained series is accelerated using the Wynn-Rho
    convergence acceleration scheme.

    Returns a ``mp.mpf`` arbitrary-precision float number, a Sequence of
    ``mp.mpf``, or an ``np.ndarray`` of ``mp.mpf`` depending upon the type of
    the ``time`` argument.

    Parameters
    ----------

    fn: Union[Callable[[float], Any], Callable[[float, int], Any]]
        The Laplace transformation to invert. Must take the Laplace parameter
        as the first argument, and optionally ``precin`` as the second argument.

        The precision is necessary for any functions that are memoized,
        otherwise the cached result will not necessarily match the input
        ``precin``.

    time: Union[float, Iterable[float], np.ndarray]
        The array of time at which to evalute ``fn``.

    M: int = 32
        The number of terms of the Gaver functional.

    precin: Optional[int] = None
        The digits of precision to use. If None (default), automatically set to
        ``round(2.1 * M)``.

    Returns
    -------
    result: Union[float, Iterable[float], np.ndarray]
        The inverted result. The type corresponds to the type of ``time``, but
        is typed as ``Any`` to avoid the requirement for type checking in
        calling code.
    """

    dps = mp.dps
    # mp.dps = int(2.1 * M) if precin is None else precin
    mp.dps = round(21 * M / 10) if precin is None else precin
    # mp.dps = max(mp.dps, MACHINE_PRECISION)

    if not isinstance(time, Iterable):
        # should be a float, but make it a catch-all for any non-Iterable
        try:
            return _gwr(fn, time, M, mp.dps)
        except Exception as e:
            raise e
        finally:
            mp.dps = dps

    if not isinstance(time, np.ndarray):
        # evaluate any Iterable that is not an np.ndarray
        try:
            return [_gwr(fn, t, M, mp.dps) for t in time]
        except Exception as e:
            raise e
        finally:
            mp.dps = dps

    # must be an ndarray or else... !!!
    assert isinstance(time, np.ndarray), f'Unknown type for time: {type(time)}'
    if time.ndim < 1:
        # to iterate over an np.ndarray it must be a vector
        try:
            return np.array([_gwr(fn, time.item(), M, mp.dps)], dtype=object)
        except Exception as e:
            raise e
        finally:
            mp.dps = dps

    if time.ndim >= 2:
        # remove single-dimensional entries from any matrix
        np.squeeze(time)
        if time.ndim >= 2:
            # cannot iterate over a matrix
            mp.dps = dps
            raise ValueError(f'Expected ndim < 2, but got {time.ndim}')

    try:
        return np.array([_gwr(fn, t, M, mp.dps) for t in time], dtype=object)
    except Exception as e:
        raise e
    finally:
        mp.dps = dps


@lru_cache(maxsize=None)
def binomial(n: int, i: int, precin: int) -> float:
    return mp.binomial(n, i)


@lru_cache(maxsize=None)
def binomial_sum(n: int, i: int, precin: int) -> float:
    if i % 2 == 1:
        return -binomial(n, i, precin)
    else:
        return binomial(n, i, precin)


@lru_cache(maxsize=None)
def fac(n: int, precin: int) -> float:
    return mp.fac(n)


@lru_cache(maxsize=None)
def fac_prod(n: int, tau: float, precin: int) -> float:
    n_fac = fac(n - 1, precin)
    return tau * fac(2 * n, precin) / (n * n_fac * n_fac)


@lru_cache(maxsize=None)
def _gwr(fn: Union[Callable[[float], Any], Callable[[float, int], Any]],
         time: float, M: int, precin: int) -> float:
    """
    GWR alorithm with memoization.
    """
    tau = mp.log(2.0) / mp.mpf(time)

    # mypy can't type check the Callable at runtime, we must do it ourselves
    sig = signature(fn).parameters
    n_params = len(sig)

    fni: List[float]
    if n_params == 1:
        fni = [fn(i * tau) if i > 0 else 0 for i in range(2 * M + 1)]  # type: ignore
    elif n_params == 2:
        fni = [fn(i * tau, precin) if i > 0 else 0 for i in range(2 * M + 1)]  # type: ignore
    else:
        raise TypeError('Too many arguments for Laplace transform. Expected 1 or 2, got '
                        f'{n_params}. Function signature:\n{sig}')

    G0: List[float] = [0.0] * M
    Gp: List[float] = [0.0] * M

    M1 = M
    for n in range(1, M + 1):
        try:
            G0[n - 1] = fac_prod(n, tau, precin) \
                * sum(binomial_sum(n, i, precin) * fni[n + i] for i in range(n + 1))

        except Exception as e:
            if n == 1:
                # we didn't perform a single iteration... something must be broken
                raise e

            M1 = n - 1
            break

    best = G0[M1 - 1]
    Gm: List[float] = [0.0] * M1

    broken = False
    for k in range(M1 - 1):
        for n in range(M1 - 1 - k)[::-1]:
            try:
                expr = G0[n + 1] - G0[n]
            except:
                # expr = 0.0
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


def _gwr_no_memo(fn: Callable[[float], Any], time: float, M: int = 32, precin: int = 0) -> float:
    """
    GWR alorithm without memoization. This is a near 1:1 translation from
    Mathematica.
    """
    tau = mp.log(2.0) / mp.mpf(time)

    fni: List[float] = [0.0] * M
    for i, n in enumerate(fni):
        if i == 0:
            continue
        fni[i] = fn(n * tau)

    G0: List[float] = [0.0] * M
    Gp: List[float] = [0.0] * M

    M1 = M
    for n in range(1, M + 1):
        try:
            n_fac = mp.fac(n - 1)
            G0[n - 1] = tau * mp.fac(2 * n) / (n * n_fac * n_fac)

            s = 0.0
            for i in range(n + 1):
                s += mp.binomial(n, i) * (-1) ** i * fni[n + i]

            G0[n - 1] *= s

        except:
            M1 = n - 1
            break

    best = G0[M1 - 1]
    Gm: List[float] = [0.0] * M1

    broken = False
    for k in range(M1 - 1):
        for n in range(M1 - 1 - k)[::-1]:
            try:
                expr = G0[n + 1] - G0[n]
            except:
                expr = 0.0
                broken = True
                break

            expr = Gm[n + 1] + (k + 1) / expr
            Gp[n] = expr
            if k % 2 == 1 and n == M1 - 2 - k:
                best = expr

        if broken:
            break

        for n in range(M1 - k):
            Gm[n] = G0[n]
            G0[n] = Gp[n]

    return best


def cache_clear() -> None:
    binomial.cache_clear()
    binomial_sum.cache_clear()
    fac.cache_clear()
    fac_prod.cache_clear()
    _gwr.cache_clear()
