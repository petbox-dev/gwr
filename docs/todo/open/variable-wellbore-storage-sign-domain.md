# Negative cDd or aD gives negative early pressure or a pole on the GWR contour

Severity: medium
Opened: 2026-10-19

## What happens now

`apply_variable_wbs` evaluates `W(s) = 1 + C_D C_Dd s / (1 + a_D s)` (the exact
algebraic form of `1 + s^2 C_D C_Dd (1/s - 1/(s + 1/a_D))`). Neither
`lap_transform` nor `VariableWellboreStorage` restricts the signs of `cDd` or
`aD`; both accept any finite value.

## Evidence

Measured by the Phase 4 review of the `/code-fix` run (commit f9590a2),
`radial_inf`, GWR M = 16, p_wD at t = 1e-6:

| cD | cDd | aD | early dp/dt | early effective storage |
|---|---|---|---|---|
| 100 | 10 | 0.01 | 999.95 | 1e-3 |
| 100 | -0.5 | 1.0 | -0.49 | negative |
| 100 | 1.0 | -0.575 | -1.73 | negative |

- Early-time effective storage is `C_D / (1 + C_D C_Dd / a_D)`. p_wD is negative
  at early time whenever `1 + C_D C_Dd / a_D <= 0`.
- `aD < 0` puts a pole at `s = 1/|aD|` on the positive real axis, where GWR
  samples. At `s = 1/0.575` exactly the ARB path returned -2.5e16 and the mpf
  path raised ZeroDivisionError.
- `Model(RadialInfinite(), VariableWellboreStorage(10.0, 1.0, -0.575)).pressure([1.0, 10.0], 8)`
  returned `[-4.26, 1.27]`.
- The docstring defines `C_Dd = (C_Ds - C_DL) / C_DL`. Under that reading the
  early storage would be `C_DL (1 + C_Dd)` (1100 for the first row), which the
  formula does not produce. Either the docstring or the formula disagrees with
  the paper.

## Where the repo already uses these values

- `test/benchmark.py` panel 4: `cD=1e7, cDd=1.0, aD=-0.575` (t >= 1e2 keeps
  s <= 0.44, below the pole).
- `test/generate_figures.py`: `CD_VWBS, CDD_VWBS, AD_VWBS = 1e3, -5.75, 1e7`;
  the curve differs from constant WBS by at most 0.06%, so the "variable WBS"
  figure shows the constant-WBS curve.

## What is still unverified

- SPE-21829-PA (Hegeman, Hallford and Joseph 1993) was not available to any
  reviewer. Whether its `C_Dd` is a ratio or an absolute change, and whether
  `a_D` is a rate or a time constant, is unknown.

## What would unblock it

Read the paper's definition of W(s), C_Dd and a_D. Then either correct the
docstring to the real early/late limits, or reject `aD <= 0` and
`1 + cD*cDd/aD <= 0` in `lap_transform` and `VariableWellboreStorage`, and fix
the benchmark and figure parameters.

## Files

- `gwr_inversion/_modifiers.py`: `apply_variable_wbs`, `WellboreEffects.apply`
- `gwr_inversion/models.py`: `VariableWellboreStorage.__post_init__`
- `gwr_inversion/transforms.py`: `lap_transform` validation
- `test/benchmark.py`, `test/generate_figures.py`: VWBS parameters
