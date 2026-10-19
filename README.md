# gwr-inversion: Laplace inversion and well-test solutions

`gwr-inversion` began as a Python reproduction of the Mathematica package
`NumericalLaplaceInversion.m`, which provides the Gaver-Wynn-Rho (GWR) function
(<https://library.wolfram.com/infocenter/MathSource/4738/>).

Version 2.0 widens the scope. The package now has three parts:

- **Numerical Laplace inversion** — `gwr()` (Gaver-Wynn-Rho with arbitrary
  precision, the method that converges on hard oscillatory transforms) and
  `talbot()` (mpmath's Fixed Talbot).
- **Reservoir solutions in the Laplace domain** — radial, linear and
  fractured-well geometries, with dual porosity, skin and constant or variable
  wellbore storage, composed into a `Model` that returns pressure, rate, the
  Bourdet derivative, cumulative and the hybrid-rate diagnostics in time.
- **Fast multiprecision special functions** — `besselk` / `besseli` and the
  fracture integrals on python-flint (ARB ball arithmetic) or gmpy2 (MPFR),
  with adaptive precision, falling back to mpmath.

A Laplace transform you write yourself is a function of the Laplace variable
`s` that uses `mpmath` (or this package's `besselk` / `besseli`). The `math`
module and `numpy` do not support multiprecision arithmetic and give invalid
results inside a transform.

| Category | Members |
|---|---|
| Inversion | [gwr](docs/api.md#gwr), [talbot](docs/api.md#talbot) |
| Models | [Model](docs/api.md#model), [LaplaceFunction](docs/api.md#laplacefunction), [OUTPUTS](docs/api.md#outputs), [Geometry](docs/api.md#geometry) |
| Model outputs | [pressure](docs/api.md#model), [rate](docs/api.md#model), [derivative](docs/api.md#model), [cumulative](docs/api.md#model), [hybrid_rate](docs/api.md#model), [laplace](docs/api.md#model), [laplace_pressure](docs/api.md#model) |
| Geometries | [RadialInfinite](docs/api.md#radialinfinite), [RadialNoFlow](docs/api.md#radialnoflow), [RadialConstantPressure](docs/api.md#radialconstantpressure), [LinearInfinite](docs/api.md#linearinfinite), [LinearNoFlow](docs/api.md#linearnoflow), [LinearConstantPressure](docs/api.md#linearconstantpressure), [FractureInfinite](docs/api.md#fractureinfinite), [FractureNoFlow](docs/api.md#fracturenoflow), [FractureConstantPressure](docs/api.md#fractureconstantpressure) |
| Effects | [DualPorosity](docs/api.md#dualporosity), [Skin](docs/api.md#skin), [WellboreStorage](docs/api.md#wellborestorage), [VariableWellboreStorage](docs/api.md#variablewellborestorage) |
| Compositor | [lap_transform](docs/api.md#lap_transform) |
| Geometry Functions | [radial_inf](docs/api.md#geometry-functions), [radial_bounded_nf](docs/api.md#geometry-functions), [radial_bounded_cp](docs/api.md#geometry-functions), [linear_inf](docs/api.md#geometry-functions), [linear_bounded_nf](docs/api.md#geometry-functions), [linear_bounded_cp](docs/api.md#geometry-functions), [frac_inf](docs/api.md#geometry-functions), [frac_bounded_nf](docs/api.md#geometry-functions), [frac_bounded_cp](docs/api.md#geometry-functions) |
| Physics Modifiers | [apply_dual_porosity](docs/api.md#physics-modifiers), [apply_skin](docs/api.md#physics-modifiers), [apply_wbs](docs/api.md#physics-modifiers), [apply_variable_wbs](docs/api.md#physics-modifiers) |
| Output Transforms | [output_pressure](docs/api.md#output-transforms), [output_rate](docs/api.md#output-transforms), [output_derivative](docs/api.md#output-transforms), [output_cumulative](docs/api.md#output-transforms) |
| Rate Analysis | [hybrid_rate](docs/api.md#hybrid_rate) |
| Special Functions | [besselk](docs/api.md#besselk--besseli), [besseli](docs/api.md#besselk--besseli) |
| Runtime | [cache_clear](docs/api.md#cache_clear), [shutdown_pool](docs/api.md#shutdown_pool), [HAS_GMPY2, HAS_FLINT, HAS_GMPY2_BESSEL](docs/api.md#has_gmpy2--has_flint--has_gmpy2_bessel) |

## Getting Started

```shell
pip install gwr_inversion                # base install
pip install gwr_inversion[fast]          # + gmpy2: faster GWR arithmetic
pip install gwr_inversion[flint]         # + python-flint: fast Bessel functions
pip install gwr_inversion[hybrid]        # + scipy: hybrid_rate / Model.hybrid_rate
pip install gwr_inversion[all]           # all three
```

Python 3.10 or later. The required dependencies are `numpy` and `mpmath`. The
optional ones are independent and complementary:

- **gmpy2** (`[fast]`) runs the GWR arithmetic (the Gaver sums and the
  Wynn-rho recursion) in MPFR, and provides a gmpy2 Bessel backend.
- **python-flint** (`[flint]`) provides ARB Bessel functions and fracture
  integrals, and evaluates `radial_inf` and the linear geometries entirely in
  ball arithmetic. Most of the speedup on reservoir solutions comes from here.
- **scipy** (`[hybrid]`) is needed only by `hybrid_rate`.

`[all]` is the recommended install.

## Examples

### Inverting a transform

```python
import math
from mpmath import mp
from gwr_inversion import gwr

def lap_log(s):
    """Laplace transform of ln(t) + 0.577216 (Euler's constant)."""
    return -mp.log(s) / s - 0.577216 / s

gwr(lap_log, time=5.0, M=32)
# 1.6094375773356333

math.log(5.0)
# 1.6094379124341003
```

`gwr` takes a scalar, a list or other iterable, or an ndarray, and returns the
same kind of container. It returns `float` values by default;
`dtype="mpf"` returns `mp.mpf` values that keep the full working precision.

### Composable models

A model is a geometry plus effects. The composition order is fixed: region
effects (dual porosity), then the geometry, then skin, then wellbore storage.

```python
import numpy as np
from gwr_inversion import Model, RadialNoFlow, DualPorosity, Skin, WellboreStorage

model = Model(RadialNoFlow(outer_radius=1000.0),
              DualPorosity(interporosity_flow=1e-6, storativity_ratio=0.1),
              Skin(3.0), WellboreStorage(100.0))
time = 10 ** np.linspace(1, 7, 61)

if __name__ == "__main__":  # required for workers > 1 on Windows and macOS (spawn)
    pressure = model.pressure(time, M=16)
    derivative = model.derivative(time, M=16)          # Bourdet derivative t dp/dt
    rate = model.rate(time, M=16, workers=4)           # constant-pressure rate
```

Geometries and effects check their parameters (`RadialNoFlow(outer_radius=0.5)`
raises `ValueError`). `model.laplace("rate")` returns the Laplace-domain
function itself, for `gwr()` or `talbot()`. A user geometry subclasses
`Geometry` and implements `laplace_pressure(s, s_eff_by_region)`; see
[docs/api.md](docs/api.md#geometry).

### Building a transform by hand

The functional layer under `Model` composes the same pieces:

```python
from gwr_inversion import gwr, lap_transform, output_rate, frac_bounded_nf

def fractured_rate(s):
    return output_rate(lap_transform(s, frac_bounded_nf, reD=50.0, skin=0.5, cD=1e-2), s)

rate = gwr(fractured_rate, [1e-2, 1e0, 1e2], M=16)
```

### Bessel functions in your own transform

Replace `mp.besselk` / `mp.besseli` with `gwr_inversion.besselk` /
`gwr_inversion.besseli`. They give the same values at the current precision and
use ARB or gmpy2 when installed:

```python
from mpmath import mp
from gwr_inversion import gwr, besselk

def line_source_with_storage(s):
    return 1.0 / (s / besselk(0, mp.sqrt(s)) + s**2 * 1e3)

gwr(line_source_with_storage, [10.0, 100.0], M=32)
```

### Parallel evaluation

```python
result = gwr(fn, time, M=32, workers=4)
```

`workers > 1` evaluates the time points in a process pool that starts on first
use (about 1.5 s on Windows) and is reused by later calls; `shutdown_pool()`
releases it. `fn` must be picklable: a module-level function or a
`model.laplace(...)` result, not a lambda. Parallel results are bit-identical
to sequential ones.

### Fixed Talbot

```python
from gwr_inversion import talbot

talbot(lambda s: 1 / (s + 1), time=1.0, degree=32)
```

Talbot converges quickly on smooth transforms but not on hard oscillatory ones,
where GWR is the method that converges. `talbot` returns `mp.mpf` values.

## Performance

Measured with `test/benchmark.py` on 2026-10-19 (31 time points, tD = 1e1 to
1e7; fractured well 21 points, tD = 1e-2 to 1e2), one run on a laptop CPU
with other load present, so single entries vary by up to about 50% between
runs. Seconds per curve. `mpmath` and `gmpy2` run the library as if
python-flint and the gmpy2 Bessel code were absent; `all` = gmpy2 +
python-flint. Results agree with pure mpmath to <= 2.4e-8 at M=8 and
<= 2e-16 at M=32.

```
                                         M=8                    M=32
Transform                         mpmath  gmpy2    all   mpmath  gmpy2    all
------------------------------------------------------------------------------
Radial inf + WBS                    0.48   0.45   0.02     4.92   5.17   0.31
Radial inf + variable WBS           0.43   0.51   0.06     5.43   6.02   0.33
Radial inf + dual porosity + WBS    1.03   1.07   0.09     4.57   4.14   0.27
Radial no-flow (reD=1000), rate     2.67   2.01   0.11    20.12  25.62   1.05
Radial const-P (reD=1000)           1.90   1.79   0.14    19.64  17.13   1.24
Radial no-flow + dual porosity      2.36   2.49   0.17    23.55  16.52   0.97
Linear no-flow (no Bessel)          0.02   0.01   0.03     0.25   0.09   0.07
Fractured well no-flow (reD=100)    2.68   2.48   0.20    35.95  33.78   1.38
```

gmpy2 alone measured 2.8x at M=32 on a transform without Bessel functions, and
about 1x on Bessel-heavy transforms, where the Bessel evaluations dominate.

## How the Bessel functions keep their precision

`besselk` and `besseli` dispatch across three backends in priority order. Each
returns a result accurate to the current `mp.prec`, and declines (passes the
argument on) what it does not cover:

1. **python-flint / ARB** — `arb.bessel_k` / `arb.bessel_i` in ball arithmetic,
   with adaptive working precision (below). Real `x > 0`.
2. **gmpy2 / MPFR** — power series for small `x`, asymptotic expansion for large
   `x`, built from MPFR operations (gmpy2 has no Bessel functions). Orders 0
   and 1, real `x > 0`.
3. **mpmath** — everything else, including complex arguments on the Talbot
   contour.

`K_n(x)` decays as `e^{-x}` while its power-series terms grow as `e^{+x}`, so
the series loses about `2 x log2(e)` bits (about 2.9 bits per unit of x) to
cancellation. Neither ARB nor the gmpy2 series compensates for this on its own.
Measured at the base precision, ARB `bessel_k(0, x)` returns no correct digits
at x = 100 (M=32) or x = 200 (M=64).

- **ARB: adaptive precision.** ARB results are balls with a rigorous error
  radius, so the loss is measurable. Evaluation starts at `mp.prec + 26` bits,
  reads `rel_accuracy_bits()`, and, if the result is short of `mp.prec + 10`
  bits, re-evaluates with the deficit plus 32 bits added.
- **ARB: K_1 from the Wronskian.** ARB's direct `bessel_k(1, x)` has a slow
  path near the switch to its asymptotic expansion (30 s or more per call at
  M=128, x = 200, python-flint 0.9.0). `K_1` comes from
  `I_0(x) K_1(x) + I_1(x) K_0(x) = 1/x` instead, which has no cancellation
  because `I_1 K_0 -> 1/(2x)`. Orders n >= 2 use the upward recurrence
  `K_{n+1} = K_{n-1} + (2n/x) K_n`.
- **gmpy2: guard bits and switch point.** The series path adds
  `2 x log2(e) + 2 log2(x + 2) + 32` bits. The asymptotic expansion has a
  minimal-term relative error near `e^{-2x}`, so it is used once
  `2 x log2(e) > mp.prec + 16`.
- **Bounded reservoirs.** The outer boundary enters through `K_n / I_n` at
  `sqrt(s) r_eD`, about `pi e^{-2z}`. It is evaluated only to the bits it
  contributes to the result, and dropped when it contributes none.

The fractured-well solutions need `int_0^z K_0(t) dt` and `int_0^z I_0(t) dt`.
These use closed forms, not quadrature: modified Struve functions as `1F2`
series for `K_0`, and `z 1F2(1/2; 1, 3/2; z^2/4)` for `I_0`.

## Changes in 2.0

- `gwr()` returns `float` / float64 arrays by default (1.x returned `mp.mpf`).
  Pass `dtype="mpf"` for the 1.x behaviour. `talbot()` still returns `mp.mpf`.
- New: `Model` and the geometry / effect classes, `lap_transform` and the
  geometry functions, physics modifiers, output transforms, `hybrid_rate`.
- New `gwr()` options: `refine` / `refine_threshold`; `workers > 1` now reuses
  a persistent pool (`shutdown_pool()`).
- `besselk` / `besseli` gained the gmpy2 backend and adaptive ARB precision;
  1.x evaluated ARB `K_n` at a fixed precision and could lose digits at large `x`.
- `cache_clear()` now empties the Bessel and Gaver caches.
- Internal names no longer exported from the package root: `get_ops`,
  `gm_to_mpf`, `gmpy2_precision`, `gwr_single`.
- Python 3.10 or later is required (1.x declared 3.7); the `[flint]` extra
  requires python-flint >= 0.9.

## Development

```shell
git clone https://github.com/petbox-dev/gwr.git
cd gwr
pip install -e ".[all]"
python -m pytest test/ -v --ignore=test/rust-bench
ruff check gwr_inversion/
mypy gwr_inversion/
```

Regenerate reference data after an algorithm change:

```shell
PYTHONPATH=. python test/generate_reference.py
```

`test/benchmark.py` produces the Performance table, and
`test/generate_figures.py` validates every geometry (M=8 against M=16) and
writes the diagnostic figures in `test/figures/`. Open work and known issues
are tracked in [docs/todo/](docs/todo/README.md).

## Benchmarks against other inversion methods

These tables come from `test/compare_methods.py` and were measured before the
2.0 GWR core, which replaced per-term factorial and binomial arithmetic with
cached exact integers. Error columns still hold; the timings are now lower and
the gmpy2 speedups smaller (see
[docs/todo/open/readme-speed-tables-predate-lean-core.md](docs/todo/open/readme-speed-tables-predate-lean-core.md)).

### Unit step: H(t-1) at t=2.0

Talbot and de Hoog converge faster than GWR on this smooth transform. Cohen
diverges. Stehfest oscillates.

```
Method                           M/degree           Error   Time (s)
---------------------------------------------------------------------
GWR (mpmath)                            8       2.533e-02      0.002
GWR (mpmath)                           32       1.487e-04      0.017
GWR (mpmath)                           64       1.191e-07      0.061
GWR (mpmath)                          128       6.220e-12      0.242
GWR (gmpy2)                             8       2.533e-02      0.001
GWR (gmpy2)                            32       1.487e-04      0.002
GWR (gmpy2)                            64       1.191e-07      0.007
GWR (gmpy2)                           128       6.220e-12      0.026
mpmath talbot                          32       0.000e+00      0.003
mpmath stehfest                        32       6.925e-03      0.004
mpmath dehoog                          32       0.000e+00      0.039
mpmath cohen                           32       6.066e+02      0.001
```

### Ln(t) across a 31-point time array

All methods give the same accuracy. For a transform this cheap the process
pool's start-up cost outweighs the parallel gain.

```
Method                                         Max Error   Time (s)
----------------------------------------------------------------------
GWR M=32 mpmath sequential                     3.351e-07      0.582
GWR M=32 mpmath parallel (4 workers)           3.351e-07      0.777
GWR M=32 gmpy2 sequential                      3.351e-07      0.126
GWR M=32 gmpy2 parallel (4 workers)            3.351e-07      0.615
Talbot degree=32                               3.351e-07      0.193
```

### Hard oscillatory: sin(t) at t=1000

GWR is the only method that converges on this problem; the others return
near-zero or diverge. GWR needs M = 1024 for full accuracy here (M = 768 is
within 1.7e-2).

```
Method                           M/degree           Error   Time (s)
---------------------------------------------------------------------
GWR (mpmath)                           32       8.269e-01      0.016
GWR (mpmath)                          128       8.269e-01      0.245
GWR (mpmath)                          512       8.269e-01      5.007
GWR (mpmath)                          768       1.671e-02    246.048
GWR (mpmath)                         1024       0.000e+00    653.559
GWR (gmpy2)                            32       8.269e-01      0.002
GWR (gmpy2)                           128       8.269e-01      0.024
GWR (gmpy2)                           512       8.269e-01      0.645
GWR (gmpy2)                           768       1.671e-02      2.226
GWR (gmpy2)                          1024       0.000e+00      5.940
mpmath talbot                         256       8.269e-01      0.032
mpmath stehfest                       256       8.269e-01      0.206
mpmath dehoog                         256       8.269e-01      2.472
mpmath cohen                          256       4.817e+18      0.008
```

## References

The method is described in: Valkó, P.P., and Abate J. 2002. Comparison of
Sequence Accelerators for the Gaver Method of Numerical Laplace Transform
Inversion. *Computers and Mathematics with Applications* **48** (3): 629–636.
<https://doi.org/10.1016/j.camwa.2002.10.017>.

More information on multi-precision inversion can be found in: Valkó, P.P. and
Vajda, S. 2002. Inversion of Noise-Free Laplace Transforms: Towards a
Standardized Set of Test Problems. *Inverse Problems in Engineering* **10** (5):
467-483. <https://doi.org/10.1080/10682760290004294>.

The reservoir solutions follow Ozkan and Raghavan (1991), Gringarten, Ramey and
Raghavan (1974), Warren and Root (1963), van Everdingen and Hurst (1949),
Hegeman, Hallford and Joseph (1993) and Fulford and Blasingame (2013); full
citations are in [docs/api.md](docs/api.md#references).
