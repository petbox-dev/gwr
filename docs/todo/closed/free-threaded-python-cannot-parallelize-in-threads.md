# Free-threaded Python cannot parallelize gwr in threads

Severity: low
Opened: 2026-10-31
Closed: 2026-10-31 -- settled by: measured; threads corrupt results (process-global precision in mpmath and python-flint) and subinterpreters cannot load numpy, gmpy2 or python-flint, so `workers` stays process-based.
If this was wrong, you will see: mpmath and python-flint gain per-thread precision contexts (the probe below prints PER THREAD for both), or numpy, gmpy2 and python-flint load in `InterpreterPoolExecutor`.

## Question

Can free-threaded CPython (3.14t) replace the process pool with threads,
or 3.14's `InterpreterPoolExecutor` replace it with subinterpreters?

## Evidence

Environments: CPython 3.14.6 and 3.14.6t; numpy 2.5.4,
mpmath 1.4.1, gmpy2 2.3.2, python-flint 0.9.0, scipy 1.18.1. The GIL stays off
after importing the package on the free-threaded build.

Precision set inside `with` blocks in two overlapping threads (3.14t):

| Library | Result |
|---|---|
| mpmath `mp.workdps` | shared: a thread that set 90 digits ran at 15 |
| python-flint `ctx.workprec` | shared: a thread that set 300 bits ran at 53 |
| gmpy2 `context` | per thread |

`_evaluate_batch_chunk` on 4 threads against sequential `evaluate`, 6 x 21
RadialNoFlow + storage table:

| Build | M | Values that differ | Worst relative error |
|---|---|---|---|
| 3.14t | 7 | 0-1 of 126 | 4e-11 |
| 3.14t | 32 | 7-17 of 126 | 1.3e-3 |

Our code changes precision mid-evaluation (the ARB adaptive loop, the boundary
ratio, `gwr`'s own `workdps`), and so can a user's transform. Threads were also
only 1.3-1.4x faster than sequential on 3.14t.

Subinterpreters (`InterpreterPoolExecutor`, 3.14.6): mpmath
imports; numpy, gmpy2 and python-flint raise ImportError ("does not support
loading in subinterpreters").

Process pool, same table, best of runs (seconds; one process / 4 workers warm):

| Interpreter | M = 7 | M = 32 |
|---|---|---|
| 3.13.16 | 0.30 / 0.10 | 2.62 / 0.77-0.87 |
| 3.14.6 | 0.30 / 0.09 | 2.63 / 0.71 |
| 3.14.6t | 0.30 / 0.12 | 2.63 / 0.92 |

The first parallel call is about 0.25 s slower from 3.14 on Linux: the default
start method became forkserver, which starts workers fresh.

## Ruled out

- Per-thread mpmath contexts: the package and users' transforms use the module
  global `mp`; a thread-local proxy would patch mpmath.
- Holding one precision for a whole parallel call: the adaptive loops must
  change it per argument.

## Files

- `gwr_inversion/_parallel.py`, `.github/workflows/ci.yml` (`test-free-threaded`)
