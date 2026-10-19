# API reference

Everything below is importable from the package root:

```python
from gwr_inversion import gwr, Model, RadialNoFlow, WellboreStorage  # and so on
```

All Laplace-domain functions take and return `mpmath` numbers (`mp.mpf`, or
`mp.mpc` on the Talbot contour). Time and every reservoir quantity are
dimensionless (`t_D`, `p_wD`, `q_D`, `r_eD`, `C_D`, ...), in the usual
well-test definitions of the cited papers.

- [Inversion](#inversion): [gwr](#gwr), [talbot](#talbot)
- [Composable models](#composable-models): [Model](#model), [LaplaceFunction](#laplacefunction), [OUTPUTS](#outputs), [Geometry](#geometry)
- [Geometries](#geometries): [RadialInfinite](#radialinfinite), [RadialNoFlow](#radialnoflow), [RadialConstantPressure](#radialconstantpressure), [LinearInfinite](#linearinfinite), [LinearNoFlow](#linearnoflow), [LinearConstantPressure](#linearconstantpressure), [FractureInfinite](#fractureinfinite), [FractureNoFlow](#fracturenoflow), [FractureConstantPressure](#fractureconstantpressure)
- [Effects](#effects): [DualPorosity](#dualporosity), [Skin](#skin), [WellboreStorage](#wellborestorage), [VariableWellboreStorage](#variablewellborestorage)
- [Laplace-domain building blocks](#laplace-domain-building-blocks): [lap_transform](#lap_transform), [geometry functions](#geometry-functions), [physics modifiers](#physics-modifiers), [output transforms](#output-transforms)
- [Rate analysis](#rate-analysis): [hybrid_rate](#hybrid_rate)
- [Special functions](#special-functions): [besselk / besseli](#besselk--besseli)
- [Runtime](#runtime): [cache_clear](#cache_clear), [shutdown_pool](#shutdown_pool), [HAS_GMPY2 / HAS_FLINT / HAS_GMPY2_BESSEL](#has_gmpy2--has_flint--has_gmpy2_bessel)

---

## Inversion

### gwr

```python
gwr(fn, time, M=32, precin=None, backend="auto", workers=1,
    dtype="float", refine=False, refine_threshold=1e-4)
```

Inverse Laplace transform by the Gaver-Wynn-Rho algorithm (Valkó and Abate 2002).

| Parameter | Meaning |
|---|---|
| `fn` | Laplace-domain function of `s` (an `mp.mpf`). It may take a second parameter, which receives `precin`. Must use `mpmath` (or this package's `besselk` / `besseli`), not `math` or `numpy`. |
| `time` | A scalar, any iterable (list, tuple, range, generator), or an ndarray. Length-1 axes are dropped; an array with two axes longer than 1 raises `ValueError`. |
| `M` | Number of Gaver terms. Accuracy and cost both grow with M; 8 to 32 covers most reservoir curves, and hard oscillatory transforms need several hundred. |
| `precin` | Decimal digits of working precision. Default `round(2.1 * M)`. |
| `backend` | `"auto"` (gmpy2 when installed), `"gmpy2"` (raises `ImportError` if absent), or `"mpmath"`. Selects the GWR arithmetic only; Bessel calls inside `fn` choose their own backend. |
| `workers` | Processes for array input. `> 1` uses a persistent process pool (see [shutdown_pool](#shutdown_pool)); `fn` must then be picklable (a module-level function or a [LaplaceFunction](#laplacefunction), not a lambda). On Windows and macOS the calling script needs an `if __name__ == "__main__":` guard. |
| `dtype` | `"float"` (default): Python `float`, float64 ndarray, or list of float. `"mpf"`: `mp.mpf` values (object ndarray for array input). |
| `refine` | Re-evaluate points whose magnitude is below `refine_threshold` times the peak at `2 M` (an explicit `precin` is raised to `round(4.2 M)` for them). Useful where a rate decays toward zero. |
| `refine_threshold` | Fraction of the peak magnitude below which `refine` acts. |

Returns the same container kind as `time`: scalar to scalar, list-like to list,
ndarray to 1-D ndarray (a 0-d array gives shape `(1,)`).

```python
from mpmath import mp
from gwr_inversion import gwr

gwr(lambda s: 1 / (s + 1), 1.0, M=32)        # e^{-1} = 0.36787944117144233
```

### talbot

```python
talbot(fn, time, degree=32)
```

Fixed Talbot inversion through `mpmath.invertlaplace(method="talbot")`. `fn`
is evaluated at complex `s`; this package's Bessel functions and transforms
accept complex arguments and fall back to `mpmath` for them. Returns `mp.mpf`
(an object ndarray for ndarray input, a list for other iterables). Converges
fast on smooth transforms; it does not converge on hard oscillatory ones
(GWR does).

---

## Composable models

### Model

```python
Model(geometry, *effects)
```

A geometry composed with effects. The composition order is fixed, whatever the
order of `effects`: region effects (dual porosity) change the Laplace variable
of a region, then the geometry gives the sandface pressure, then skin, then
wellbore storage.

Validation at construction (`ValueError` or `TypeError`): the geometry is a
[Geometry](#geometry); at most one `Skin`; at most one storage effect
(`WellboreStorage` or `VariableWellboreStorage`); at most one `DualPorosity`
per region, and the region must exist in the geometry. Parameter domains are
checked by each geometry and effect.

Models are frozen and pickle, so they work with `workers > 1`.

| Method | Returns |
|---|---|
| `pressure(time, M=32, **gwr_options)` | `p_wD(t)` |
| `rate(time, M=32, **gwr_options)` | `q_D(t)` at constant wellbore pressure (`1 / (s^2 p̄_wD)`) |
| `derivative(time, M=32, **gwr_options)` | Bourdet derivative `t dp_wD/dt`. `dtype="float"` only. |
| `cumulative(time, M=32, **gwr_options)` | `∫_0^t p_wD dτ` (`p̄_wD / s`). For cumulative production apply the same division to the rate. |
| `hybrid_rate(time, M=32, *, b=0.5, **gwr_options)` | [hybrid_rate](#hybrid_rate) applied to `rate(time)`. Needs scipy. |
| `laplace(output="pressure")` | A picklable [LaplaceFunction](#laplacefunction) for `gwr` or `talbot`. |
| `laplace_pressure(s)` | `p̄_wD(s)` at one Laplace point. |

`gwr_options` are passed to [gwr](#gwr) (`precin`, `backend`, `workers`,
`dtype`, `refine`, `refine_threshold`).

```python
import numpy as np
from gwr_inversion import Model, RadialNoFlow, DualPorosity, Skin, WellboreStorage

model = Model(RadialNoFlow(outer_radius=1000.0),
              DualPorosity(interporosity_flow=1e-6, storativity_ratio=0.1),
              Skin(3.0), WellboreStorage(100.0))
time = 10 ** np.linspace(1, 7, 61)
pressure = model.pressure(time, M=16)
derivative = model.derivative(time, M=16)
```

### LaplaceFunction

```python
LaplaceFunction(model, output)
```

`s -> output(model.laplace_pressure(s), s)`, as one picklable callable. Made by
`Model.laplace(output)`; constructing one with an `output` not in
[OUTPUTS](#outputs) raises `ValueError`.

### OUTPUTS

`("pressure", "rate", "derivative", "cumulative")`: the names
`Model.laplace` accepts.

### Geometry

Base class for a user geometry. Set the class attribute `regions` (default
`("reservoir",)`) and implement:

```python
def laplace_pressure(self, s, s_eff_by_region) -> mpf
```

`s_eff_by_region` maps each region name to its effective Laplace variable (`s`
itself where no region effect applies). The model applies skin and storage
after it. User geometries take the generic path; the built-in geometries also
have a ball-arithmetic fast path.

---

## Geometries

All are frozen dataclasses. Each rejects parameters outside the domain shown
(`ValueError`, which also rejects NaN and infinity).

| Class | Parameters (default) | Domain | Laplace function |
|---|---|---|---|
| <a id="radialinfinite"></a>`RadialInfinite()` | none | | [radial_inf](#geometry-functions) |
| <a id="radialnoflow"></a>`RadialNoFlow(outer_radius)` | `outer_radius` r_eD = r_e / r_w (1000) | > 1 | [radial_bounded_nf](#geometry-functions) |
| <a id="radialconstantpressure"></a>`RadialConstantPressure(outer_radius)` | `outer_radius` (1000) | > 1 | [radial_bounded_cp](#geometry-functions) |
| <a id="linearinfinite"></a>`LinearInfinite(distance)` | `distance` x_D from the well (0) | >= 0 | [linear_inf](#geometry-functions) |
| <a id="linearnoflow"></a>`LinearNoFlow(distance)` | `distance` x_D (0) | [0, 1] | [linear_bounded_nf](#geometry-functions) |
| <a id="linearconstantpressure"></a>`LinearConstantPressure(distance)` | `distance` x_D (0) | [0, 1] | [linear_bounded_cp](#geometry-functions) |
| <a id="fractureinfinite"></a>`FractureInfinite(fracture_position)` | `fracture_position` x_D along the half-length (0.732) | [0, 1] | [frac_inf](#geometry-functions) |
| <a id="fracturenoflow"></a>`FractureNoFlow(outer_radius, fracture_position)` | r_eD = r_e / x_f (100), x_D (0.732) | > 1, [0, 1] | [frac_bounded_nf](#geometry-functions) |
| <a id="fractureconstantpressure"></a>`FractureConstantPressure(outer_radius, fracture_position)` | r_eD (100), x_D (0.732) | > 1, [0, 1] | [frac_bounded_cp](#geometry-functions) |

`fracture_position = 0.732` evaluates the uniform-flux solution at the point
where it equals the infinite-conductivity fracture (Gringarten et al. 1974);
`0` gives the uniform-flux wellbore pressure.

## Effects

| Class | Parameters | Domain | Level |
|---|---|---|---|
| <a id="dualporosity"></a>`DualPorosity(interporosity_flow, storativity_ratio, region="reservoir")` | λ, ω (Warren and Root 1963, pseudo-steady state) | λ > 0, 0 < ω < 1 | region |
| <a id="skin"></a>`Skin(skin)` | S | finite | wellbore |
| <a id="wellborestorage"></a>`WellboreStorage(storage)` | C_D | > 0 | wellbore |
| <a id="variablewellborestorage"></a>`VariableWellboreStorage(storage, storage_change_ratio, change_rate)` | C_D, C_Dd, a_D (Hegeman et al. 1993) | C_D > 0; others finite | wellbore |

The signs of C_Dd and a_D are not restricted; see
[docs/todo/open/variable-wellbore-storage-sign-domain.md](todo/open/variable-wellbore-storage-sign-domain.md).

---

## Laplace-domain building blocks

The functional layer under `Model`, for writing transforms by hand.

### lap_transform

```python
lap_transform(s, geometry_fn, cD=0.0, cDd=None, aD=None, lam=None, omega=None,
              skin=0.0, **geom_kwargs)
```

`p̄_wD(s)` composed in the standard order: `s_eff` from dual porosity (when
`lam` and `omega` are given), `geometry_fn(s, s_eff, **geom_kwargs)`, then skin,
then constant (`cD`) or variable (`cD`, `cDd`, `aD`) wellbore storage. For
`radial_inf` and the linear geometries the whole chain runs in ARB ball
arithmetic when python-flint is installed (same value, faster).

Raises `ValueError` when only one of `lam`/`omega` or of `cDd`/`aD` is given,
when `lam <= 0` or `omega` is outside (0, 1), when `cD` is negative or NaN, or
when variable storage is given with `cD = 0`. Geometry parameters (`reD`, `xD`)
are not checked here; `Model` checks them.

```python
from gwr_inversion import gwr, lap_transform, output_rate, radial_bounded_nf

rate = gwr(lambda s: output_rate(lap_transform(s, radial_bounded_nf, reD=500.0, cD=100.0), s),
           [1e2, 1e4, 1e6], M=16)
```

### Geometry functions

Signature `(s, s_eff, **parameters) -> p̄_D`, the sandface pressure before skin
and storage. `s_eff = s` without dual porosity.

| Function | Parameters | Solution |
|---|---|---|
| `radial_inf` | | Finite wellbore, infinite reservoir: `K_0(√s_eff) / (s √s_eff K_1(√s_eff))` |
| `radial_bounded_nf` | `reD=1000.0` | Closed (no-flow) circular boundary (Ozkan and Raghavan 1991) |
| `radial_bounded_cp` | `reD=1000.0` | Constant-pressure circular boundary |
| `linear_inf` | `xD=0.0` | `exp(-√s_eff x_D) / (s √s_eff)` |
| `linear_bounded_nf` | `xD=0.0` | `cosh(√s_eff (1 - x_D)) / (s √s_eff sinh √s_eff)` |
| `linear_bounded_cp` | `xD=0.0` | `sinh(√s_eff (1 - x_D)) / (s √s_eff cosh √s_eff)` |
| `frac_inf` | `xD=0.732` | Uniform-flux vertical fracture (Gringarten et al. 1974) |
| `frac_bounded_nf` | `reD=100.0, xD=0.732` | Fracture in a closed circular reservoir |
| `frac_bounded_cp` | `reD=100.0, xD=0.732` | Fracture with a constant-pressure circular boundary |

The bounded solutions evaluate the outer-boundary ratio `K_n / I_n` at
`√s_eff r_eD` only to the bits it contributes to the result, and drop it when it
contributes none.

### Physics modifiers

| Function | Result |
|---|---|
| `apply_dual_porosity(s, lam, omega)` | `s_eff = s (λ + ω(1-ω) s) / (λ + (1-ω) s)` |
| `apply_skin(psd, s, skin)` | `psd + S / s` |
| `apply_wbs(psd, s, cD)` | `1 / (1/psd + s² C_D)` |
| `apply_variable_wbs(psd, s, cD, cDd, aD)` | `W(s) / (1/psd + s² C_D)`, `W = 1 + C_D C_Dd s / (1 + a_D s)` |

They use only `+ - * / **`, so they also work on ARB balls.

### Output transforms

Signature `(pwd, s) -> value`.

| Function | Laplace value | Inverse |
|---|---|---|
| `output_pressure` | `p̄_wD` | `p_wD(t)` |
| `output_rate` | `1 / (s² p̄_wD)` | `q_D(t)` at constant wellbore pressure |
| `output_derivative` | `s p̄_wD` | `dp_wD/dt`; multiply by `t` for the Bourdet derivative |
| `output_cumulative` | `p̄_wD / s` | `∫_0^t p_wD dτ` |

---

## Rate analysis

### hybrid_rate

```python
qD_hybrid, qDd, qDdr, qDdi, qDdid = hybrid_rate(qD, time, b=0.5)
```

Arps-corrected hybrid rate and its diagnostic functions (Fulford and
Blasingame 2013), in the real domain. `qD` and `time` are 1-D arrays of the
same length (at least 2; time positive and strictly increasing).

| Output | Definition |
|---|---|
| `qD_hybrid` | `q_D + 1/(1 + b t)^(1/b) - e^(-t)` |
| `qDd` | `t |dq/dt|` (numpy gradient) |
| `qDdr` | `1 / (t d(1/q)/dt)` |
| `qDdi` | `(1/t) ∫_0^t q dτ`; the integral over `[0, t_0]` extrapolates `q ~ tⁿ` through the first two samples |
| `qDdid` | `qDdi - qD_hybrid` |

Raises `ValueError` for `b <= 0`, mismatched or short arrays, a time grid that
is not positive and increasing, first rates that are not positive, or first
rates that fall faster than `1/t` (start the grid earlier, in transient flow).
Needs scipy (`pip install gwr-inversion[hybrid]`); raises `ImportError` naming
the extra otherwise.

---

## Special functions

### besselk / besseli

```python
besselk(n, x)    # K_n(x)
besseli(n, x)    # I_n(x)
```

Drop-in replacements for `mp.besselk` / `mp.besseli` in a transform. Accurate to
the current `mp.prec`. Backends in order, each declining what it does not cover:

1. ARB (python-flint): real `x > 0`; integer `n` for K, real `n` for I. Adaptive
   working precision against the K series cancellation; K_1 from the Wronskian.
2. gmpy2 MPFR: real `x > 0`, `n` in {0, 1}.
3. mpmath: everything else (complex `x`, other orders).

---

## Runtime

### cache_clear

`cache_clear()` empties the internal caches: the order-0/1 Bessel values and the
Gaver integer tables (about 54 MB at M = 1024, per process).

### shutdown_pool

`shutdown_pool()` stops the worker pool that `gwr(..., workers > 1)` starts on
first use and reuses while the worker count is unchanged. The next parallel
call starts a new pool. Also registered with `atexit`.

### HAS_GMPY2 / HAS_FLINT / HAS_GMPY2_BESSEL

Booleans: whether gmpy2 and python-flint imported, and whether the gmpy2 Bessel
backend is available. Read-only indicators, fixed at import.

---

## References

- Valkó, P. P. and Abate, J. 2002. Comparison of Sequence Accelerators for the
  Gaver Method of Numerical Laplace Transform Inversion. *Computers and
  Mathematics with Applications* 48 (3): 629-636.
  <https://doi.org/10.1016/j.camwa.2002.10.017>
- Ozkan, E. and Raghavan, R. 1991. New Solutions for Well-Test-Analysis
  Problems: Parts 1 and 2. *SPE Formation Evaluation* 6 (3): 359-378.
  <https://doi.org/10.2118/18615-PA>, <https://doi.org/10.2118/18616-PA>
- Gringarten, A. C., Ramey, H. J., Jr., and Raghavan, R. 1974. Unsteady-State
  Pressure Distributions Created by a Well With a Single Infinite-Conductivity
  Vertical Fracture. *SPE Journal* 14 (4): 347-360.
  <https://doi.org/10.2118/4051-PA>
- Warren, J. E. and Root, P. J. 1963. The Behavior of Naturally Fractured
  Reservoirs. *SPE Journal* 3 (3): 245-255. <https://doi.org/10.2118/426-PA>
- van Everdingen, A. F. and Hurst, W. 1949. The Application of the Laplace
  Transformation to Flow Problems in Reservoirs. *Journal of Petroleum
  Technology* 1 (12): 305-324. <https://doi.org/10.2118/949305-G>
- Hegeman, P. S., Hallford, D. L., and Joseph, J. A. 1993. Well-Test Analysis
  With Changing Wellbore Storage. *SPE Formation Evaluation* 8 (3): 201-207.
  <https://doi.org/10.2118/21829-PA>
- Fulford, D. S. and Blasingame, T. A. 2013. Evaluation of Time-Rate
  Performance of Shale Wells using the Transient Hyperbolic Relation.
  SPE-167242-MS. <https://doi.org/10.2118/167242-MS>
