"""Parallel time-point evaluation using a persistent process pool."""
from __future__ import annotations

import atexit
import os
import threading
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from typing import Any, Callable, List, Optional

# One pool per process, created on first use and reused by later calls. Starting
# a pool costs about 1-1.5 s on Windows (spawn re-imports numpy, mpmath, gmpy2,
# python-flint and this package in every worker), so a pool per call made short
# curves slower than sequential evaluation.
_pool: Optional[ProcessPoolExecutor] = None
_pool_workers = 0
_pool_lock = threading.Lock()


def _evaluate_chunk(
    times: List[float], fn: Callable[..., Any], M: int, precin: int, backend: str,
) -> List[Any]:
    """Evaluate gwr_single at every time in ``times``. Top-level for pickling."""
    from mpmath import mp

    from ._gmpy2 import get_ops, gmpy2_precision
    from ._gwr import gwr_single

    mp.dps = precin
    ops, _ = get_ops(backend)
    with gmpy2_precision(precin):
        return [gwr_single(fn, time, M, precin, ops) for time in times]


def _get_pool(workers: int) -> ProcessPoolExecutor:
    """Return the shared pool, (re)creating it when the worker count changes."""
    global _pool, _pool_workers
    with _pool_lock:
        if _pool is None or _pool_workers != workers:
            if _pool is not None:
                _pool.shutdown(wait=True)
            _pool = ProcessPoolExecutor(max_workers=workers)
            _pool_workers = workers
        return _pool


def shutdown_pool() -> None:
    """Shut down the shared worker pool. The next ``workers > 1`` call starts a new one."""
    global _pool, _pool_workers
    with _pool_lock:
        if _pool is not None:
            _pool.shutdown(wait=True)
        _pool = None
        _pool_workers = 0


atexit.register(shutdown_pool)


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

    The time points are split into one interleaved chunk per worker
    (``times[i::n]``), so cheap and expensive parts of the curve are spread
    evenly, and each worker receives ``fn`` once per call.

    Note: fn must be picklable (module-level function, not a lambda).
    """
    pool_size = min(workers, os.cpu_count() or 1)
    chunk_count = min(pool_size, len(times))
    chunks = [times[i::chunk_count] for i in range(chunk_count)]
    pool = _get_pool(pool_size)
    try:
        futures = [pool.submit(_evaluate_chunk, chunk, fn, M, precin, backend)
                   for chunk in chunks]
        parts = [future.result() for future in futures]
    except BrokenProcessPool:
        # A worker died; drop the pool so the next call starts a fresh one.
        shutdown_pool()
        raise
    results: List[Any] = [None] * len(times)
    for i, part in enumerate(parts):
        results[i::chunk_count] = part
    return results
