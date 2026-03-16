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
pip install gwr_inversion[fast]          # + gmpy2: ~10x speedup on GWR internals
pip install gwr_inversion[flint]         # + python-flint: ~15x speedup on Bessel functions
pip install gwr_inversion[all]           # both (maximum speed)
```

The two optional dependencies accelerate different parts of the computation and
are complementary:

- **gmpy2** (`[fast]`) — replaces the GWR algorithm's internal arithmetic
  (factorial, binomial, log2, multiprecision multiply/add) with MPFR. Provides
  ~6x speedup on transforms that do not call Bessel functions.

- **python-flint** (`[flint]`) — provides fast ARB-based Bessel functions via
  `gwr_inversion.besselk` and `gwr_inversion.besseli`. Provides ~3x speedup on
  Bessel-heavy transforms even when gmpy2 is also installed.

Speedup summary on a 31-point time array (M=32):

```
                              mpmath    gmpy2    flint    gmpy2+flint
--------------------------------------------------------------------
ln(t)  — no Bessel            0.53s    0.09s    0.53s       0.09s
Radial flow — Bessel-heavy    1.83s    1.34s    0.62s       0.18s
```

`[all]` is always at least as fast as either alone and is the recommended
install for production use.

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

# After (ARB besselk, ~15x faster):
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
python-flint) instead of `mp.besselk` gives ~15x total speedup with gmpy2.

```python
from gwr_inversion import gwr, besselk  # besselk uses ARB if python-flint installed

def lap_rad(s):
    b = besselk(0, mp.sqrt(s))
    return 1.0 / (s / b + s**2 * 1e3)
```

```
Method                              Time (s)   Max Rel Diff vs reference
------------------------------------------------------------------------
GWR M=32 mpmath besselk (ref)          2.738   (reference)
GWR M=32 ARB besselk + mpmath          0.610   0.000e+00
GWR M=32 ARB besselk + gmpy2           0.181   1.355e-16
Talbot degree=8                        0.168   9.191e-07
Talbot degree=16                       0.382   1.536e-11
Talbot degree=32                       0.963   0.000e+00
Talbot degree=64                       2.626   0.000e+00
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
