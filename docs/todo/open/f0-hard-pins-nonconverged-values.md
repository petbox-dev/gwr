# test_f0_hard pins wrong values at M=768 and vacuous ones at M=256/512, and never reads sin_expected

Severity: low
Opened: 2026-10-19

## What happens now

`test/test_f0_hard.py` inverts `1/(s^2 + 1)` (sin t) at t = 1000 for
M = 256, 512, 768 and 1024 and compares each result with the stored
`sin_m<M>` value at `ATOL = 1e-12`.

## Evidence

Values in `test/reference_data.npz` (read in the `/code-fix` run):

| key | value | truth sin(1000) = 0.82687954 |
|---|---|---|
| sin_m256 | -8.1e-97 | not converged |
| sin_m512 | 3.6e-71 | not converged |
| sin_m768 | 0.84358923 | off by 1.67e-2 |
| sin_m1024 | 0.82687954 | correct |
| sin_expected | 0.82687954 | stored, never read |

- M = 256 and 512 are vacuous at ATOL 1e-12: any result below 1e-12 passes.
- M = 768 pins a wrong value as correct: an accuracy improvement would fail the
  suite.
- The file and the pinned values are unchanged since BASE v1.2.0 (cdb3c36).

## What would unblock it

Assert `abs(gwr(..., M=1024) - sin_expected) < tol`; for M = 256-768 either
drop the cases or assert the documented non-convergence
(`abs(result - sin_expected) > 1e-3`). The README's statement "GWR requires
M >= 768" should then read M >= 1024 for full accuracy.

## Files

- `test/test_f0_hard.py`, `test/generate_reference.py` (lines writing `sin_*`)
