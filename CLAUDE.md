# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`gwr-inversion` — Python package for numerical inversion of Laplace transforms. Implements the Gaver-Wynn-Rho algorithm (Valkó & Abate, 2004) and wraps mpmath's Fixed Talbot method.

Public API: `gwr(fn, time, M=32, precin=None, backend="auto", workers=1, dtype="float", refine=False, refine_threshold=1e-4)`, `talbot(fn, time, degree=32)`, `besselk(n, x)`, `besseli(n, x)`, `Model(geometry, *effects)` (see `models.py`).

## Build & Install

```bash
pip install -e .
pip install -e ".[fast]"   # gmpy2: ~3x on GWR arithmetic (M=32; README tables)
pip install -e ".[flint]"  # python-flint: ~15x speedup on Bessel functions
pip install -e ".[hybrid]" # scipy: needed by hybrid_rate / Model.hybrid_rate
pip install -e ".[all]"    # all three
```

Dependencies: `numpy>=1.21.3`, `mpmath>=1.1.0`. Optional: `gmpy2>=2.2.0`, `python-flint>=0.9.0`, `scipy>=1.7.2` (floors tested by the `minimum-versions` CI job). Requires Python >=3.10 (python-flint 0.9 needs it). CI: `.github/workflows/ci.yml`.

## Testing

```bash
python -m pytest test/ -v
```

Regression tests compare against pre-computed reference data in `test/reference_data.npz`. Regenerate with `PYTHONPATH=. python test/generate_reference.py`. The f0_hard tests (M=256-1024) take ~10 s with gmpy2 and ~1.5-2.5 minutes with gmpy2 absent (92 s and 135 s in two runs, 2026-10-21). The reference data was generated with the gmpy2 backend; five tests pin its rounding and skip without gmpy2 (`pins_gmpy2_rounding` in `test/backend_markers.py`; see `docs/todo/open/reference-tests-require-gmpy2.md`).

Linting/type checking:
```bash
ruff check gwr_inversion/
mypy gwr_inversion/
```

## Architecture

```
gwr_inversion/
├── __init__.py    # Public API re-exports (full list: docs/api.md)
├── _gmpy2.py      # gmpy2/MPFR backend: ops dict + gmpy2_precision context manager
├── _bessel.py     # Bessel dispatch: ARB -> gmpy2 -> mpmath
├── _flint.py      # python-flint/ARB backend: besselk(), besseli(), K0/I0 integrals
├── _numeric.py    # Shared numeric policy: positive_real(), guard bits, log2(e)
├── _gmpy2_bessel.py  # gmpy2/MPFR K_n, I_n (n in {0, 1}): series + asymptotic
├── transforms.py  # Composable Laplace-domain geometries, outputs, compositor
├── _modifiers.py  # Physics modifiers (dual porosity, skin, WBS); mpf and ARB alike
├── models.py      # Composable Model objects: geometries, effects, outputs
├── _arb_transforms.py  # Ball-arithmetic lap_transform for radial_inf + linear geometries
├── _gwr.py        # gwr() / talbot() entry points and the GWR algorithm (gwr_single)
├── _talbot.py     # Fixed Talbot wrapper (via mpmath.invertlaplace)
├── _parallel.py   # Parallel time-point dispatch (ProcessPoolExecutor)
```

- **`__init__.py`** — Re-exports the public API (full list: `docs/api.md`). `gwr()` (in `_gwr.py`) handles input polymorphism (float/list/ndarray), precision context, backend selection, and parallel dispatch. Returns `float` / float64 `ndarray` by default (`dtype="float"`); `dtype="mpf"` returns `mp.mpf` regardless of backend.
- **`_gmpy2.py`** — Provides `get_ops(backend)` returning a dict of math operations (`log2`, `mpf`, `to_mpf`, `from_mpf`); the Gaver coefficients and binomials are exact ints cached in `_gwr._gaver_integers`. When gmpy2 is installed and `backend="auto"`, uses MPFR via gmpy2. Also provides `gmpy2_precision(precin)` context manager to set gmpy2 context precision for all arithmetic within a block.
- **`_flint.py`** — Provides `besselk(n, x)` and `besseli(n, x)` using ARB (python-flint) when available. Raises ARB working precision until `rel_accuracy_bits()` reaches `mp.prec + 10` (ARB does not compensate the K_n series cancellation of ~`2x log2(e)` bits). `K_1` comes from the Wronskian (ARB's direct `bessel_k(1)` has a slow path near x = 50-300). Also `besselk0_integral` / `besseli0_integral` via Struve/1F2 closed forms. Pure ARB: every function returns None when it declines (flint absent, argument outside `positive_real`, no convergence).
- **`_bessel.py`** — `besselk`/`besseli`, the K0/I0 integrals and `bessel_k_over_i`: the one place that owns the fallback chain (ARB, then gmpy2, then mpmath; a backend returning None passes the argument on). Users call `besselk`/`besseli` in their Laplace transform functions. Complex arguments (Talbot) always reach mpmath.
- **`_numeric.py`** — `positive_real(x)` (the domain the fast backends accept), `TARGET_GUARD_BITS` and the other precision-loop constants, `LOG2_E`. Switch a backend off at run time with `_flint.HAS_FLINT` or `_gmpy2_bessel.HAS_GMPY2_BESSEL` (read at call time; see `test/benchmark.py`).
- **`_gmpy2_bessel.py`** — gmpy2 K_0, K_1, I_0, I_1. Series path adds at least `2x log2(e)` guard bits; asymptotic expansion once `2x log2(e) > mp.prec + 16`. Working precision derives from `mp.prec`, not the gmpy2 context.
- **`_gwr.py`** — `gwr_single(fn, time, M, precin, ops)` implements the GWR algorithm for one time point using the ops dict. User's `fn` always receives mpmath types; backend conversion happens transparently.
- **`_parallel.py`** — `parallel_gwr()` dispatches one interleaved chunk of time points per worker to a persistent `ProcessPoolExecutor` (created on first use, reused while the pool size `min(workers, cpu_count)` is unchanged, `shutdown_pool()` / atexit to release). Requires picklable `fn` (module-level functions, not lambdas).
- **`_arb_transforms.py`** — `lap_transform` with `radial_inf` or a linear geometry runs entirely in ARB balls (one conversion in, one out). Bounded radial and fractured geometries stay on the mpf path (measured no gain in ARB); their outer boundary enters through `_bessel.bessel_k_over_i`, evaluated only to the bits it contributes.
- **`_modifiers.py`** — `apply_dual_porosity`, `apply_skin`, `apply_wbs`, `apply_variable_wbs`, `effective_s` (step 1) and `WellboreEffects` (skin then storage, steps 3-4). Plain `+ - * / **`, so the mpf compositor, the ARB compositor and the generic `Model` path share one copy; `transforms` re-exports them.
- **`models.py`** — `Model(geometry, *effects)` frozen dataclass; geometries declare `regions`; effects have a level (region: `DualPorosity`; wellbore: `Skin`, `WellboreStorage`, `VariableWellboreStorage`). Built-in geometries delegate to `transforms.compose` (keeps the ARB path; the `WellboreEffects` bundle is built once per model); user `Geometry` subclasses use the generic path. `model.laplace(output)` returns a picklable `LaplaceFunction`.
- **`_talbot.py`** — Thin wrapper around `mpmath.invertlaplace(method='talbot')`.

Key patterns:
- `inspect.signature()` detects if `fn` accepts a `precin` parameter
- Graceful degradation: reduces M if precision is insufficient
- Backend ops dict avoids import-time gmpy2 dependency

## Outstanding work

Deferred work, unexplained anomalies and open decisions live in `docs/todo/` (one
document per item, with its evidence). Read `docs/todo/README.md` before adding an
item; the row in its open table is part of adding one. Never leave `# TODO:` in code.

The public API is documented in `docs/api.md`; the README's API table links into it.
Keep both in step with `gwr_inversion/__init__.py` when an export changes.

## Code Standards

- Max line length: 100 (ruff)
- Strict mypy enabled
- Reference Mathematica implementation: `gwr.m`
