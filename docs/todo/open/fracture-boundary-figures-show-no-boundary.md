# frac_pressure.png and frac_rate.png plot three identical curves: the time range never reaches the boundary

Severity: low
Opened: 2026-10-19

## What happens now

`test/generate_figures.py` draws a "boundary comparison" for the fractured
well: `frac_inf`, `frac_bounded_nf` and `frac_bounded_cp` on `FRAC_TIME`
(1e-2 to 1e3) with `FRAC_RED = 100`. The outer boundary is felt only near
t_D ~ r_eD^2 = 1e4, so all three panels show the infinite-acting curve.

## Evidence

- In the validation run (M = 8 vs M = 16) the three fractured rows had
  identical errors: 1.32e-06 for pressure and 1.93e-05 for rate.
- The committed `test/figures/frac_pressure.png` was opened in review: all three
  panels are the same curve. `frac_rate.png` has the same cause.

## Blast radius

Diagnostic figures only. A reader concludes the bounded fractured solutions
equal the infinite one.

## What would unblock it

Extend `FRAC_TIME` past 1e4, or set `FRAC_RED` to about 10, then regenerate the
figures with `python test/generate_figures.py` and commit the new PNGs (they are
deliberately tracked).

## Files

- `test/generate_figures.py`: `FRAC_TIME`, `FRAC_RED`, `_boundary_figures`
- `test/figures/frac_pressure.png`, `test/figures/frac_rate.png`
