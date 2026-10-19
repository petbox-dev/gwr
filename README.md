# Gaver-Wynn-Rho Algorithm

This is a Python reproduction of the Mathematica package that provides the GWR
function, `NumericalLaplaceInversion.m`.

<https://library.wolfram.com/infocenter/MathSource/4738/>

This package provides numerical inverse Laplace transform via the Gaver-Wynn-Rho
(GWR) algorithm and mpmath's Fixed Talbot method. The GWR algorithm uses
arbitrary-precision arithmetic and is the only known method that converges on
hard oscillatory transforms.

The Laplace transform should be provided as a function that uses the `mpmath`
library for a scalar value of the Laplace parameter.  The `math` library and
`numpy` functions do not support multiprecision math and will return invalid
results if they are used.

## Getting Started

```shell
pip install gwr_inversion                # base install
pip install gwr_inversion[fast]          # + gmpy2: faster GWR internals
pip install gwr_inversion[flint]         # + python-flint: fast Bessel functions
pip install gwr_inversion[all]           # both (maximum speed)
```

The two optional dependencies accelerate different parts of the computation and
are complementary:

- **gmpy2** (`[fast]`) — replaces the GWR algorithm's internal arithmetic
  (log2, multiprecision multiply/add) with MPFR. Measured 2.8x at M=32 on a
  transform without Bessel functions; about 1x on Bessel-heavy transforms,
  where the Bessel evaluations dominate.

- **python-flint** (`[flint]`) — provides fast ARB-based Bessel functions via
  `gwr_inversion.besselk` and `gwr_inversion.besseli`, closed-form fracture
  integrals, and ball-arithmetic evaluation of `radial_inf` and the linear
  geometries. This is where most of the speedup on reservoir transforms comes
  from.

Measured with `test/benchmark.py` on 2026-10-19 (31 time points, tD = 1e1 to
1e7; fractured well 21 points, tD = 1e-2 to 1e2), one run on a laptop CPU
with other load present, so single entries vary by up to about 50% between
runs. Seconds per curve. `mpmath` and `gmpy2` run the library as if
python-flint and the gmpy2 Bessel code were absent; `all` = gmpy2 +
python-flint. Results agree with pure mpmath to <= 2.3e-8 at M=8 and
<= 1.8e-16 at M=32.

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

`[all]` is the recommended install for production use.

## Development

```shell
git clone https://github.com/petbox-dev/gwr.git
cd gwr
pip install -e ".[all]"                  # editable install with all optional deps
```

Run the test suite:

```shell
python -m pytest test/ -v --ignore=test/rust-bench
```

Regenerate reference data (required after algorithm changes):

```shell
PYTHONPATH=. python test/generate_reference.py
```

Lint and type-check:

```shell
ruff check gwr_inversion/
mypy gwr_inversion/
```

## Simple Example

```python
import math
from gwr_inversion import gwr
from mpmath import mp

def lap_log_fn(s: float):
    # log function
    return -mp.log(s) / s - 0.577216 / s

gwr(lap_log_fn, time=5.0, M=32)
# mpf('1.6094375773356333')

math.log(5.0)
# 1.6094379124341003
```

## Optional Speedups

### gmpy2 backend — transparent, no code changes

The gmpy2 backend accelerates GWR's internal arithmetic (factorial, binomial,
log2) using MPFR. It is selected automatically when installed and requires no
changes to user code:

```python
# These are all equivalent once gmpy2 is installed:
gwr(fn, time, M=32)                    # auto-selects gmpy2
gwr(fn, time, M=32, backend="auto")    # explicit auto
gwr(fn, time, M=32, backend="gmpy2")   # force gmpy2
gwr(fn, time, M=32, backend="mpmath")  # force mpmath
```

### python-flint — fast Bessel functions

For transforms that call Bessel functions, replace `mp.besselk` / `mp.besseli`
with `gwr_inversion.besselk` / `gwr_inversion.besseli`. These use ARB when
python-flint is installed and fall back to mpmath automatically:

```python
from mpmath import mp
from gwr_inversion import gwr, besselk

# Before (mpmath besselk):
def lap_rad_slow(s):
    b = mp.besselk(0, mp.sqrt(s))
    return 1.0 / (s / b + s**2 * 1e3)

# After (ARB besselk; see the measured table under Getting Started):
def lap_rad_fast(s):
    b = besselk(0, mp.sqrt(s))
    return 1.0 / (s / b + s**2 * 1e3)

result = gwr(lap_rad_fast, time, M=32)
```

### Backend selection

```python
# Force a specific backend (overrides auto-select):
gwr(fn, time, backend="mpmath")   # pure mpmath
gwr(fn, time, backend="gmpy2")    # require gmpy2 (raises ImportError if absent)
gwr(fn, time, backend="auto")     # default: use gmpy2 if installed
```

### Bessel function implementation

`gwr_inversion.besselk` and `gwr_inversion.besseli` dispatch across three
backends in priority order. Each backend returns a result accurate to the
current `mp.prec`:

1. **python-flint/ARB** (`_flint.py`) — calls `arb.bessel_k` / `arb.bessel_i`
   from the FLINT library's ball arithmetic, with adaptive working precision
   (see below).

2. **gmpy2/MPFR** (`_gmpy2_bessel.py`) — custom implementations composed from
   fundamental MPFR operations since gmpy2 does not provide Bessel functions
   natively. Uses power series for small x and asymptotic expansion for large x.

3. **mpmath** — fallback when neither optional backend is available.

`K_n(x)` needs careful precision management: it decays as `e^{-x}` while its
power-series terms grow as `e^{+x}`, so the series loses about `2 x log2(e)`
bits (about 2.9 bits per unit of x) to cancellation. Neither ARB nor the gmpy2
series compensates for this on its own. Measured at the base precision, ARB
`bessel_k(0, x)` returns no correct digits at x = 100 (M=32) or x = 200 (M=64).

- **ARB: adaptive precision.** ARB results are balls with a rigorous error
  radius, so the loss is measurable. `_flint.py` evaluates at `mp.prec + 26`
  bits, reads `rel_accuracy_bits()`, and, if the result is short of
  `mp.prec + 10` bits, re-evaluates with the deficit plus 32 bits added. One
  retry is typical for 5 < x < the asymptotic switch point; none outside it.
- **ARB: K_1 from the Wronskian.** ARB's direct `bessel_k(1, x)` has a slow
  path near the switch to its asymptotic expansion (30 s or more per call at
  M=128, x = 200, python-flint 0.9.0). `K_1` is computed from
  `I_0(x) K_1(x) + I_1(x) K_0(x) = 1/x` instead, which has no cancellation
  because `I_1 K_0 -> 1/(2x)`. Orders n >= 2 use the upward recurrence
  `K_{n+1} = K_{n-1} + (2n/x) K_n`.
- **gmpy2: guard bits and switch point.** The series path adds
  `2 x log2(e) + 2 log2(x + 2) + 32` bits. The asymptotic expansion
  `K_n(x) ~ sqrt(pi/2x) * e^{-x} * sum(a_k / x^k)` has a minimal-term relative
  error near `e^{-2x}`, so it is used once `2 x log2(e) > mp.prec + 16`. The
  gmpy2 backend also computes `K_1` via the Wronskian.

`besseli` does not suffer from cancellation (all series terms are positive for
real x > 0) and needs no precision boost.

The fractured-well transforms need `int_0^z K_0(t) dt` and `int_0^z I_0(t) dt`.
These use closed forms, not numerical quadrature:

- `int_0^z K_0 = (pi z / 2) [K_0(z) L_{-1}(z) + K_1(z) L_0(z)]`, with the
  modified Struve functions `L_nu` evaluated as `1F2` series in ARB
  (`arb.hypgeom`), or `mp.struvel` without python-flint.
- `int_0^z I_0 = z 1F2(1/2; 1, 3/2; z^2/4)` in ARB, or `mp.besseli(0, z, -1)`
  (derivative order -1 is the integral) without python-flint.

## Parallel Evaluation

For array inputs, time points can be evaluated in parallel:

```python
import numpy as np

# Must use a module-level function (not a lambda) for parallel
def my_transform(s):
    return mp.exp(-s) / s

time = np.logspace(-1, 3, 100)
result = gwr(my_transform, time, M=32, workers=4)
```

Note: `workers > 1` requires `fn` to be picklable (module-level function or
class method, not a lambda or closure).

## Fixed Talbot Method

For well-behaved transforms, the Fixed Talbot method is also available:

```python
from gwr_inversion import talbot

result = talbot(lambda s: 1 / (s + 1), time=1.0, degree=32)
```

## Benchmarks

The gmpy2 backend (`pip install gwr_inversion[fast]`) provides ~10x speedup on
GWR algorithm internals by using MPFR instead of mpmath. For transforms that
call mpmath functions internally (e.g. `mp.besselk`), the speedup is smaller
since those calls are not accelerated.

### Unit step: H(t-1) at t=2.0

Talbot and de Hoog converge faster than GWR on this smooth transform. Cohen
diverges. Stehfest oscillates. GWR (gmpy2) matches mpmath accuracy at ~10x
the speed.

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

### Ln(t) across 31-point time array

All methods produce identical accuracy. gmpy2 provides ~4.6x sequential
speedup; parallel evaluation provides additional speedup for large arrays.

```
Method                                         Max Error   Time (s)
----------------------------------------------------------------------
GWR M=32 mpmath sequential                     3.351e-07      0.582
GWR M=32 mpmath parallel (4 workers)           3.351e-07      0.777
GWR M=32 gmpy2 sequential                      3.351e-07      0.126
GWR M=32 gmpy2 parallel (4 workers)            3.351e-07      0.615
Talbot degree=32                               3.351e-07      0.193
```

### Radial flow with wellbore storage (besselk)

The primary petroleum engineering test case. No closed-form solution — GWR and
Talbot agree to full precision. Using `gwr_inversion.besselk` (ARB via
python-flint) instead of `mp.besselk` gives 5.2x total speedup with gmpy2 on
this transform (measured 2026-10-19).

```python
from gwr_inversion import gwr, besselk  # besselk uses ARB if python-flint installed

def lap_rad(s):
    b = besselk(0, mp.sqrt(s))
    return 1.0 / (s / b + s**2 * 1e3)
```

```
Method                              Time (s)   Max Rel Diff vs reference
------------------------------------------------------------------------
GWR M=32 mpmath besselk (ref)          2.115   (reference)
GWR M=32 ARB besselk + mpmath          0.499   0.000e+00
GWR M=8 gmpy2                          0.380   1.423e-06
GWR M=16 gmpy2                         0.754   1.888e-11
GWR M=32 gmpy2                         2.116   1.355e-16
GWR M=32 ARB besselk + gmpy2           0.407   1.355e-16
Talbot degree=8                        0.258   9.191e-07
Talbot degree=16                       0.616   1.536e-11
Talbot degree=32                       1.426   0.000e+00
Talbot degree=64                       3.969   0.000e+00
```

### Hard oscillatory: sin(t) at t=1000

GWR is the **only method that converges** on this problem. All others return
near-zero or diverge. The gmpy2 backend provides ~10-100x speedup here.

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

GWR requires M≥768 to converge on this problem. The gmpy2 backend reduces
the M=1024 run from ~11 minutes to 6 seconds (~110x speedup).

Run `python test/compare_methods.py` to reproduce all benchmarks.

## References

The method is described in: Valkó, P.P., and Abate J. 2002. Comparison of
Sequence Accelerators for the Gaver Method of Numerical Laplace Transform
Inversion. *Computers and Mathematics with Application* **48** (3): 629–636.
<https://doi.org/10.1016/j.camwa.2002.10.017>.

More information on multi-precision inversion can be found in: Valkó, P.P. and
Vajda, S. 2002. Inversion of Noise-Free Laplace Transforms: Towards a
Standardized Set of Test Problems. *Inverse Problems in Engineering* **10** (5):
467-483. <https://doi.org/10.1080/10682760290004294>.
