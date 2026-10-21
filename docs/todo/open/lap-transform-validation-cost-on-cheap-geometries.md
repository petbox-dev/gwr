# Direct lap_transform on linear_bounded_nf is 1.0-1.7 us/call (10-15%) slower than v2.0 HEAD c14c218

Severity: low
Opened: 2026-10-19

## What happens now

The `/code-fix` run (f9590a2) added parameter validation to `lap_transform`
(lam/omega pairing and domain, cDd/aD pairing, `cD >= 0` including NaN, `cD > 0`
when storage is variable), builds a `WellboreEffects` per call, and routes
through `compose`. On the cheapest geometry (no Bessel functions) this is a
measurable share of the call.

## Evidence

`hot_path_timing.py` (scratch script, not kept; s-grid not recorded), 200
points, dps 17, best of 7,
back-to-back with a `git archive` of c14c218, three rounds:

| case | c14c218 (us/call) | f9590a2 (us/call) |
|---|---|---|
| lap linear_bounded_nf + WBS | 10.2-11.4 | 11.3-12.8 |
| lap radial_inf + WBS | 33.8-39.3 | 38.3-40.2 (noise) |
| lap radial_bounded_nf | 93.9-100.7 | 90.2-99.2 (noise) |
| Model LinearNoFlow + WBS | 12.1-13.3 | 11.7-12.8 |

- An earlier per-call mpf-parameter scan cost about 1 us more; it was replaced
  by an exception path inside the ARB loop.
- `Model.laplace_pressure` does not pay the validation (the Model validates
  once at construction) and is at parity or faster.
- At the `gwr_single` level a `radial_bounded_nf` call is about 97 us per s at
  M = 8, so there the overhead is about 1-2%. On `linear_bounded_nf` a call is
  26-31 us per s (2026-10-21 re-run), so the 1-2 us overhead is about 3-7%.
- A 2026-10-21 re-run with a log-spaced grid reproduced the linear delta
  (+1.5 to +2 us); its `radial_inf` absolutes (53-67 us) differ, so the grid
  differs from the original.

## What is still unverified

Whether the cost matters for any real workload; no end-to-end benchmark was
re-run.

## What would unblock it

If it matters: validate once in a small object built per transform (as `Model`
does) rather than per call, or move the checks behind a module-level flag.
Otherwise close it as accepted.

## Files

- `gwr_inversion/transforms.py`: `lap_transform`, `compose`
- `gwr_inversion/_arb_transforms.py`: `lap_transform`
