"""Parallel time-point evaluation using concurrent.futures."""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from typing import Any, Callable, List


def _worker(time: float, fn: Callable[..., Any], M: int, precin: int, backend: str) -> Any:
    """Worker function for parallel evaluation. Must be top-level for pickling."""
    from mpmath import mp
    from ._gmpy2 import get_ops, gmpy2_precision
    from ._gwr import gwr_single

    mp.dps = precin
    ops, _ = get_ops(backend)
    with gmpy2_precision(precin):
        return gwr_single(fn, time, M, precin, ops)


def parallel_gwr(
    fn: Callable[..., Any],
    times: List[float],
    M: int,
    precin: int,
    backend: str,
    workers: int,
) -> List[Any]:
    """
    Evaluate gwr_single across time points in parallel.

    Note: fn must be picklable (module-level function, not a lambda).
    """
    max_workers = min(workers, len(times), os.cpu_count() or 1)
    worker = partial(_worker, fn=fn, M=M, precin=precin, backend=backend)
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        return list(pool.map(worker, times))
