# A raising chunk leaves the other chunks running in the persistent pool; untested on fork platforms

Severity: low
Opened: 2026-10-19

## What happens now

`parallel_gwr` submits one chunk per worker to a persistent
`ProcessPoolExecutor` and collects `future.result()` in order. If one chunk
raises (or the user presses Ctrl-C), the remaining futures are not cancelled.
They keep running, and the next call's chunks queue behind them. Before the
persistent pool (BASE v1.2.0) a per-call `with` block waited for them instead.

## Evidence

- Executed in the `/code-fix` Phase 3 plumbing audit: a raising `fn` in a
  worker gives `RuntimeError` in the caller; the pool recovered in every case
  tried (raise, unpicklable lambda, `os._exit` giving `BrokenProcessPool` then a
  fresh pool).
- `_get_pool` returns the pool after releasing its lock, so a concurrent
  `shutdown_pool` can make `submit` raise "cannot schedule new futures after
  shutdown" (reproduced by hand-interleaving). `mp.workdps` already mutates
  global mpmath state, so concurrent `gwr` calls from threads are unsafe anyway.
- A call with `workers > 1` and a single time point still goes through the pool
  (1-2 s cold start, needs a picklable `fn`); this was the same at BASE.

## What is still unverified

With the fork start method (the Linux default to Python 3.13; Python 3.14 and
later use forkserver, macOS has used spawn since 3.8) workers are forked at first use. A
function defined after that first parallel call (REPL, notebook) is unknown to
the workers and may fail to unpickle; a fork after pool creation inherits
`_pool`. Neither case is tested on any platform, and the forkserver behaviour
(Linux, Python 3.14+) for them is not examined.

## What would unblock it

Cancel the remaining futures in an `except BaseException` block in
`parallel_gwr`; take the sequential path when `len(times) < 2`; run the
`test_parallel.py` suite and a REPL scenario on Linux.

## Files

- `gwr_inversion/_parallel.py`: `parallel_gwr`, `_get_pool`, `shutdown_pool`
