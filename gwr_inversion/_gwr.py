"""GWR algorithm implementation using pluggable backend."""

from inspect import signature
from typing import Any, Callable, Dict, List


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

    M1 = M
    for n in range(1, M + 1):
        try:
            n_fac = ops["fac"](n - 1, precin)
            coeff = tau * ops["fac"](2 * n, precin) / (n * n_fac * n_fac)

            s = zero
            for i in range(n + 1):
                binom = ops["binomial"](n, i, precin)
                if i % 2 == 1:
                    s = s - binom * fni[n + i]
                else:
                    s = s + binom * fni[n + i]

            G0[n - 1] = coeff * s

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
