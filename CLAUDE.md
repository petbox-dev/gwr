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

Dependencies: `numpy>=1.17`, `mpmath>=1.1.0`. Optional: `gmpy2>=2.1.0`, `python-flint>=0.9.0`. Requires Python >=3.7.

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
├── _bessel.py     # Bessel dispatch: ARB -> gmpy2 -> mpmath
├── _flint.py      # python-flint/ARB backend: besselk(), besseli(), K0/I0 integrals
├── _gmpy2_bessel.py  # gmpy2/MPFR K_n, I_n (n in {0, 1}): series + asymptotic
├── transforms.py  # Composable Laplace-domain geometries, modifiers, outputs
├── _gwr.py        # Core GWR algorithm (gwr_single)
├── _talbot.py     # Fixed Talbot wrapper (via mpmath.invertlaplace)
├── _parallel.py   # Parallel time-point dispatch (ProcessPoolExecutor)
```

- **`__init__.py`** — Entry points. `gwr()` handles input polymorphism (float/list/ndarray), precision context, backend selection, and parallel dispatch. Returns `float` / float64 `ndarray` by default (`dtype="float"`); `dtype="mpf"` returns `mp.mpf` regardless of backend.
- **`_gmpy2.py`** — Provides `get_ops(backend)` returning a dict of math operations (`log2`, `mpf`, `fac`, `binomial`, `to_mpf`, `from_mpf`). When gmpy2 is installed and `backend="auto"`, uses MPFR via gmpy2. Also provides `gmpy2_precision(precin)` context manager to set gmpy2 context precision for all arithmetic within a block.
- **`_flint.py`** — Provides `besselk(n, x)` and `besseli(n, x)` using ARB (python-flint) when available. Raises ARB working precision until `rel_accuracy_bits()` reaches `mp.prec + 10` (ARB does not compensate the K_n series cancellation of ~`2x log2(e)` bits). `K_1` comes from the Wronskian (ARB's direct `bessel_k(1)` has a slow path near x = 50-300). Also `besselk0_integral` / `besseli0_integral` via Struve/1F2 closed forms. Falls back to mpmath.
- **`_bessel.py`** — `besselk`/`besseli` dispatch (ARB, then gmpy2, then mpmath); users call these in their Laplace transform functions.
- **`_gmpy2_bessel.py`** — gmpy2 K_0, K_1, I_0, I_1. Series path adds `2x log2(e)` guard bits; asymptotic expansion once `2x log2(e) > mp.prec + 16`. Working precision derives from `mp.prec`, not the gmpy2 context.
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
