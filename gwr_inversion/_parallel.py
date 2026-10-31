"""Parallel time-point evaluation using a persistent process pool."""
from __future__ import annotations

import atexit
import os
import threading
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from typing import Any

# One pool per process, created on first use and reused by later calls. Starting
# a pool costs about 1-1.5 s on Windows (spawn re-imports numpy, mpmath, gmpy2,
# python-flint and this package in every worker), so a pool per call made short
# curves slower than sequential evaluation.
_pool: ProcessPoolExecutor | None = None
_pool_workers = 0
_pool_lock = threading.Lock()


def _evaluate_chunk(
    times: list[float], fn: Callable[..., Any], M: int, precin: int, backend: str,
) -> list[Any]:
    """Evaluate gwr_single at every time in ``times``. Top-level for pickling."""
    from mpmath import mp

    from ._gmpy2 import get_ops, gmpy2_precision
    from ._gwr import gwr_single

    ops = get_ops(backend)
    with mp.workdps(precin), gmpy2_precision(precin):
        return [gwr_single(fn, time, M, precin, ops) for time in times]


def _evaluate_batch_chunk(
    times: list[float], fns: list[Callable[..., Any]], M: int, precin: int, backend: str,
) -> list[list[Any]]:
    """
    gwr_single for every function at every time in ``times``, time-major: all
    functions at one time before the next, so they share that time's 2M Bessel
    nodes through the order-0/1 cache while the entries are fresh. One list of
    per-function results per time. Top-level for pickling.
    """
    from mpmath import mp

    from ._gmpy2 import get_ops, gmpy2_precision
    from ._gwr import gwr_single

    ops = get_ops(backend)
    with mp.workdps(precin), gmpy2_precision(precin):
        return [[gwr_single(fn, time, M, precin, ops) for fn in fns] for time in times]


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


def _interleaved_chunks(times: list[float], pool_size: int) -> list[list[float]]:
    """
    One interleaved chunk per worker (``times[i::n]``), so cheap and expensive
    parts of the curve spread evenly; never more chunks than time points, so no
    worker gets an empty chunk.
    """
    chunk_count = min(pool_size, len(times))
    return [times[i::chunk_count] for i in range(chunk_count)]


def _map_interleaved(
    evaluate_chunk: Callable[..., list[Any]],
    times: list[float],
    arguments: tuple[Any, ...],
    workers: int,
) -> list[Any]:
    """
    ``evaluate_chunk(chunk, *arguments)`` on one interleaved chunk of ``times``
    per worker, in the shared pool; the per-time results back in time order.

    Interleaving (``times[i::n]``) spreads cheap and expensive parts of the
    curve evenly. Grouping times whose GWR nodes coincide (frequent on evenly
    spaced grids) cuts duplicate Bessel work but measured slower: the load
    balance matters more (docs/todo/closed/parallel-chunks-ignore-shared-nodes.md).
    """
    pool_size = min(workers, os.cpu_count() or 1)
    chunks = _interleaved_chunks(times, pool_size)
    chunk_count = len(chunks)
    pool = _get_pool(pool_size)
    try:
        futures = [pool.submit(evaluate_chunk, chunk, *arguments) for chunk in chunks]
        parts = [future.result() for future in futures]
    except BrokenProcessPool:
        # A worker died; drop the pool so the next call starts a fresh one.
        shutdown_pool()
        raise
    results: list[Any] = [None] * len(times)
    for i, part in enumerate(parts):
        results[i::chunk_count] = part
    return results


def parallel_gwr(
    fn: Callable[..., Any],
    times: list[float],
    M: int,
    precin: int,
    backend: str,
    workers: int,
) -> list[Any]:
    """
    Evaluate gwr_single across time points in parallel.

    The time points are split into one interleaved chunk per worker
    (``times[i::n]``), so cheap and expensive parts of the curve are spread
    evenly, and each worker receives ``fn`` once per call.

    Note: fn must be picklable (module-level function, not a lambda).
    """
    return _map_interleaved(_evaluate_chunk, times, (fn, M, precin, backend), workers)


def parallel_gwr_batch(
    fns: list[Callable[..., Any]],
    times: list[float],
    M: int,
    precin: int,
    backend: str,
    workers: int,
) -> list[list[Any]]:
    """
    ``parallel_gwr`` for several functions in one dispatch: each worker takes
    one interleaved chunk of the times and evaluates every function at each of
    them (``_evaluate_batch_chunk``). Returns one list of per-function results
    per time. Every function must be picklable.
    """
    return _map_interleaved(_evaluate_batch_chunk, times, (fns, M, precin, backend), workers)
