"""
GWR(F, t, M, precin) gives the inverse of the Laplace transform function named
'F' for a given time point 't'. The method involves the calculation of 'M' terms
of the Gaver-functional. The obtained series is accelerated using the Wynn-Rho
convergence acceleration scheme. The precision of internal calculations is set
to 'precin'.

GWR(F, t, M) does the same, but the precision of the internal calculations is
selected automatically:  precin = 2.1 M).

GWR(F, t) uses M = 32 terms and precin = 67 as defaults. It should give
reasonable results for many problems.

Important note: The Laplace transform should be degfined as a function of one
argument. It can involve anything from a simple Mathematica expression to a
sophisticated Module. Since the Laplace transform will be evaluated with
non-standard (multiple) precision, approximate numbers (with decimal point) or
Mathematica functions starting with the letter 'N' are not allowed in the
function definition.

Example usage:

def fun(s):
    (1/s) * mp.exp(-mp.sqrt( (s ** 2 + (37 / 100) s + 1) / (s ** 2 + s + Pi)))

t0 = 100
GWR(fun, t0)
"""
import numpy as np
import mpmath  # type: ignore
from mpmath import mp

from typing import List, Dict, Tuple, Any, Callable, Union


MACHINE_PRECISION = 15


LOG2 = mpmath.log(2.0)


def gwr(fn: Callable[[float], Any], time: np.ndarray, M: int = 32, precin: int = 0) -> np.ndarray:
    time = np.atleast_1d(time)
    return np.array([_gwr(fn, t, M, precin) for t in time], dtype=object)


def _gwr(fn: Callable[[float], Any], time: float, M: int = 32, precin: int = 0) -> float:
    tau = LOG2 / mpmath.mpf(time)
    mp.dps = max(21 * M / 10 if precin <= 0 else precin, MACHINE_PRECISION)
    broken = False

    fni = np.arange(0, 2 * M + 1, dtype=object)
    for i, n in enumerate(fni):
        if i == 0:
            continue
        fni[i] = fn(n * tau)

    G0 = np.empty(M, dtype=object)
    Gp = np.empty(M, dtype=object)

    M1 = M
    for n in range(1, M + 1):
        try:
            n_fac = mpmath.fac(n - 1)
            G0[n - 1] = tau * mpmath.fac(2 * n) / (n * n_fac * n_fac)

            s = 0.0
            for i in range(n + 1):
                s += mpmath.binomial(n, i) * (-1) ** i * fni[n + i]
            # s = mpmath.fsum(
            #     lambda i: mpmath.binomial(n, i) * (-1) ** i * fni[n + int(i)], [0, n + 1])

            G0[n - 1] *= s

        except:
            M1 = n - 1
            break

    best = G0[M1 - 1]
    Gm = np.full(M1, 0.0, dtype=object)

    for k in range(M1 - 1):
        for n in range(M1 - 1 - k)[::-1]:
            try:
                expr = G0[n + 1] - G0[n]
            except:
                broken = True

            if broken or expr == 0.0:
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

    return mpmath.mpf(best)
