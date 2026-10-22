# Regression tests pin gmpy2 rounding at ATOL 1e-12 and fail on an install without gmpy2

Severity: low
Opened: 2026-10-19

## What happens now

`test_ln.py`, `test_unit_step.py`, `test_radial_flow.py` and
`test_input_types.py` compare `gwr()` with the default backend against
`test/reference_data.npz` with `assert_allclose(..., atol=1e-12)` (its default
`rtol=1e-7` also applies; `test_input_types.py::test_scalar` uses a bare
`< 1e-12`). The reference was generated with
the gmpy2 backend. gmpy2 is an optional dependency (`[fast]`), so
`pip install gwr-inversion` followed by pytest runs the mpmath backend.

Since 2026-10-21 the 5 tests that fail on mpmath carry
`pins_gmpy2_rounding` (`test/backend_markers.py`) and skip without gmpy2. They
still check one backend's rounding at M = 4 and 8, not accuracy.

## Evidence

Executed in the `/code-fix` run (f9590a2), gmpy2 blocked, same inputs:

| reference key | M | mpmath vs reference |
|---|---|---|
| ln_m4 | 4 | 3.5e-5 |
| unit_m8 | 8 | 5.4e-6 |
| list_result | 4 | 1.6e-6 |
| rad_pressure | 8 | 1.1e-8 |

- 5 tests fail at HEAD c14c218 with gmpy2 blocked; 6 failed at BASE v1.2.0
  (cdb3c36), including `rad_derivative`. Re-run 2026-10-21 (gmpy2, flint and
  scipy blocked): the 5 are `test_ln::test_m4`, `test_unit_step::test_m8` and
  `test_input_types::test_scalar` / `test_list` / `test_ndarray`.
  `rad_pressure` differs by 1.1e-8 but passes on the `rtol` term. The other 13
  tests in those files and all 4 in `test_f0_hard.py` pass.
- Cause: `precision_digits = round(2.1 M)` is 8 digits at M = 4. The two
  backends round differently: `_gmpy2_bits` adds 10 guard bits to the
  digits-to-bits conversion and mpmath's own conversion adds about 4, a net 6
  bits (dps 8: 36 vs 30 bits). At `list_result[0]` (t = 1, M = 4) mpmath gives
  -4.3892e-4, gmpy2 -4.3735e-4; both are off -4.3453e-4, the M = 4 value at
  `precin >= 64` (the true value is -3.35e-7).
- So the tests check bitwise reproduction of one backend's rounding noise, not
  accuracy. `test_backend.py` compares the two backends only at M = 32, atol 1e-10.
- The reference was regenerated twice in this branch (76212ca, 1013d3a); the
  BASE file differs from the current one by up to 5.2e-4 (ln_m4).

## Blast radius

Before the marker: anyone running the suite without gmpy2 saw failures that
are not regressions. Now: those 5 checks do not run on the mpmath backend at
all, so an mpmath-path regression at low M would pass unseen there.

## What would unblock it

The skip marker is in place (one of the three options). Still open: replace
the pinned rounding with a comparison against a high-precision value
(`precin=256`) at an honest tolerance, or store a per-backend reference, so the
tests check accuracy on both backends. Separately decide whether `mp_mpf` should
carry the same guard bits as `_gmpy2_bits`.

## Files

- `test/test_ln.py`, `test/test_unit_step.py`, `test/test_radial_flow.py`,
  `test/test_input_types.py`, `test/generate_reference.py`
- `gwr_inversion/_gmpy2.py`: `_gmpy2_bits`, `mp_mpf`
- `test/test_f0_hard.py` also loads `reference_data.npz`; it passes without
  gmpy2 (135 s).
- `test/backend_markers.py`: `pins_gmpy2_rounding`. Remove it from the 5 tests
  when this item closes.
