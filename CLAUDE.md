# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`gwr-inversion` — Python package for numerical inversion of Laplace transforms. Implements the Gaver-Wynn-Rho algorithm (Valkó & Abate, 2002) and wraps mpmath's Fixed Talbot method.

Public API: `gwr(fn, time, M=32, precin=None, backend="auto", workers=1)`, `talbot(fn, time, degree=32)`, `besselk(n, x)`, `besseli(n, x)`.

## Build & Install

```bash
pip install -e .
pip install -e ".[fast]"   # gmpy2: ~10x speedup on GWR internals
pip install -e ".[flint]"  # python-flint: ~15x speedup on Bessel functions
pip install -e ".[all]"    # both
```

Dependencies: `numpy>=1.17`, `mpmath>=1.1.0`. Optional: `gmpy2>=2.1.0`, `python-flint>=0.6.0`. Requires Python >=3.7.

## Testing

```bash
python -m pytest test/ -v --ignore=test/rust-bench
```

Regression tests compare against pre-computed reference data in `test/reference_data.npz`. Regenerate with `PYTHONPATH=. python test/generate_reference.py`. The f0_hard tests (M=256-1024) take ~10 minutes.

Linting/type checking:
```bash
ruff check gwr_inversion/
mypy gwr_inversion/
```

## Architecture

```
gwr_inversion/
├── __init__.py    # Public API: gwr(), talbot(), besselk(), besseli(), cache_clear()
├── _gmpy2.py      # gmpy2/MPFR backend: ops dict + gmpy2_precision context manager
├── _flint.py      # python-flint/ARB backend: besselk(), besseli() helpers
├── _gwr.py        # Core GWR algorithm (gwr_single)
├── _talbot.py     # Fixed Talbot wrapper (via mpmath.invertlaplace)
├── _parallel.py   # Parallel time-point dispatch (ProcessPoolExecutor)
```

- **`__init__.py`** — Entry points. `gwr()` handles input polymorphism (float/list/ndarray), precision context, backend selection, and parallel dispatch. Returns `mp.mpf` regardless of backend.
- **`_gmpy2.py`** — Provides `get_ops(backend)` returning a dict of math operations (`log2`, `mpf`, `fac`, `binomial`, `to_mpf`, `from_mpf`). When gmpy2 is installed and `backend="auto"`, uses MPFR via gmpy2. Also provides `gmpy2_precision(precin)` context manager to set gmpy2 context precision for all arithmetic within a block.
- **`_flint.py`** — Provides `besselk(n, x)` and `besseli(n, x)` using ARB (python-flint) when available (~15x faster than mpmath). Falls back to mpmath. Users call these in their Laplace transform functions.
- **`_gwr.py`** — `gwr_single(fn, time, M, precin, ops)` implements the GWR algorithm for one time point using the ops dict. User's `fn` always receives mpmath types; backend conversion happens transparently.
- **`_parallel.py`** — `parallel_gwr()` dispatches time points across `ProcessPoolExecutor` workers. Requires picklable `fn` (module-level functions, not lambdas).
- **`_talbot.py`** — Thin wrapper around `mpmath.invertlaplace(method='talbot')`.

Key patterns:
- `inspect.signature()` detects if `fn` accepts a `precin` parameter
- Graceful degradation: reduces M if precision is insufficient
- Backend ops dict avoids import-time gmpy2 dependency

## Code Standards

- Max line length: 100 (ruff)
- Strict mypy enabled
- Reference Mathematica implementation: `gwr.m`
