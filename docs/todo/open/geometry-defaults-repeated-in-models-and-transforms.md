# Geometry default values are repeated in models.py and transforms.py with nothing keeping them equal

Severity: low
Opened: 2026-10-21

## What happens now

Three defaults are unnamed literals, written once in each `transforms` function
and again in the matching `models` dataclass:

| value | transforms.py | models.py |
|---|---|---|
| x_D = 0.732 (fracture position) | `frac_inf`, `frac_bounded_nf`, `frac_bounded_cp` | `FractureInfinite`, `FractureNoFlow`, `FractureConstantPressure` |
| r_eD = 100.0 (fracture outer radius) | `frac_bounded_nf`, `frac_bounded_cp` | `FractureNoFlow`, `FractureConstantPressure` |
| r_eD = 1000.0 (radial outer radius) | `radial_bounded_nf`, `radial_bounded_cp` | `RadialNoFlow`, `RadialConstantPressure` |

`Model` always passes the parameters explicitly, so the two sets are
independent: changing one leaves the other silently different.

## Evidence

- `grep` of `0.732`, `= 100.0`, `= 1000.0` in `gwr_inversion/` (2026-10-21).
- Only `FractureInfinite()` is pinned, by `test/test_models.py`
  (`{"xD": 0.732}`), and that test passes `xD` explicitly rather than reading
  the transform's default.

## What would unblock it

Name the three values once (for example in `transforms.py`), use the names in
both modules, and add one test comparing each dataclass default with
`inspect.signature(transforms.<function>).parameters[...].default`.

## Files

- `gwr_inversion/transforms.py`, `gwr_inversion/models.py`
- `test/test_models.py`
