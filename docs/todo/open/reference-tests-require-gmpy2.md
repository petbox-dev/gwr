# Regression tests pin gmpy2 rounding at ATOL 1e-12 and fail on an install without gmpy2

Severity: medium
Opened: 2026-10-19

## What happens now

`test_ln.py`, `test_unit_step.py`, `test_radial_flow.py` and
`test_input_types.py` compare `gwr()` with the default backend against
`test/reference_data.npz` at `ATOL = 1e-12`. The reference was generated with
the gmpy2 backend. gmpy2 is an optional dependency (`[fast]`), so
`pip install gwr-inversion` followed by pytest runs the mpmath backend.

## Evidence

Executed in the `/code-fix` run (f9590a2), gmpy2 blocked, same inputs:

| reference key | M | mpmath vs reference |
|---|---|---|
| ln_m4 | 4 | 3.5e-5 |
| unit_m8 | 8 | 5.4e-6 |
| list_result | 4 | 1.6e-6 |
| rad_pressure | 8 | 1.1e-8 |

- 5 tests fail at HEAD c14c218 with gmpy2 blocked; 6 failed at BASE v1.2.0
  (cdb3c36), including `rad_derivative`.
- Cause: `precision_digits = round(2.1 M)` is 8 digits at M = 4. The two
  backends round differently (`_gmpy2_bits` adds 10 guard bits; `mp_mpf` adds
  none). At `ln_m4` mpmath gives -4.3892e-4, gmpy2 -4.3735e-4; both are off the
  converged value -4.3453e-4, which both reach at `precin >= 64`.
- So the tests check bitwise reproduction of one backend's rounding noise, not
  accuracy. `test_backend.py` compares the two backends only at M = 32, atol 1e-10.
- The reference was regenerated twice in this branch (76212ca, 1013d3a); the
  BASE file differs from the current one by up to 5.2e-4 (ln_m4).

## Blast radius

Anyone running the suite without gmpy2 (CI without the extra, a minimal
install) sees failures that are not regressions.

## What would unblock it

Choose one: compare against a high-precision value (`precin=256`) with an
honest tolerance; or mark the reference tests
`skipif(not HAS_GMPY2)` and say they pin the gmpy2 backend; or store a
per-backend reference. Separately decide whether `mp_mpf` should carry the same
+10 guard bits as `_gmpy2_bits`.

## Files

- `test/test_ln.py`, `test/test_unit_step.py`, `test/test_radial_flow.py`,
  `test/test_input_types.py`, `test/generate_reference.py`
- `gwr_inversion/_gmpy2.py`: `_gmpy2_bits`, `mp_mpf`
