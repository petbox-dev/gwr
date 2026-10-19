# gm_mpf builds mpfr(str(t)): the gmpy2 backend inverts the decimal string, mpmath the binary double

Severity: low
Opened: 2026-10-19

## What happens now

`_gmpy2.gm_mpf(val, precin)` returns `mpfr(str(val))`. `gwr_single` uses it for
`tau = ln 2 / t`. The mpmath backend's `mp_mpf` uses `mp.mpf(val)`, the exact
binary value of the float. For a non-dyadic float such as 0.1 the two backends
invert slightly different times.

## Evidence

Measured in the `/code-fix` run (f9590a2), radial WBS cD = 100, M = 32,
relative difference mpmath vs gmpy2:

| t | relative difference |
|---|---|
| 0.1 | 5.5e-17 |
| 0.3 | 3.7e-17 |
| 7.7 | 2.3e-17 |
| 0.5 (dyadic) | 2.4e-21 |

- A second probe with `f(t) = t^5` (slope 5) measured 2.78e-16 at t = 0.1.
- The two backends also accept different input types: `Fraction(1, 3)` gives
  `ValueError: invalid digits` on gmpy2 and `TypeError` on mpmath;
  `np.float32` and `np.int64` work on gmpy2 and raise on mpmath.
- `gm_mpf` is unchanged since BASE v1.2.0 (cdb3c36). `gm_from_mpf` already has
  the exact mantissa/exponent path.

## What is still unverified

Whether any test or user depends on the decimal-string behaviour. The
reference tests (`reference-tests-require-gmpy2`) pin gmpy2 output at 1e-12,
well above these differences.

## What would unblock it

`return gm_from_mpf(mp.mpf(val), precin)` in `gm_mpf`, then decide whether to
regenerate `reference_data.npz`.

## Files

- `gwr_inversion/_gmpy2.py`: `gm_mpf`, `gm_from_mpf`, `mp_mpf`
- `gwr_inversion/_gwr.py`: `gwr_single` (`tau`)
