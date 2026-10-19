# A constant Gaver sequence gives nan on gmpy2 and ZeroDivisionError on mpmath; a bare except hides the class

Severity: low
Opened: 2026-10-19

## What happens now

In `gwr_single`, the Gaver loop catches `Exception` to shorten `M` when a term
fails, and the Wynn-rho loop wraps only `expr = G0[n + 1] - G0[n]` in a bare
`except:`. The division `Gp[n] = Gm[n + 1] + (k + 1) / expr` sits outside the
`try`, so a zero difference reaches it.

## Evidence

Executed in the `/code-fix` Phase 3 numeric audit:

- `fn = 0 * s` and `fn = 1 + 0 * s`: the gmpy2 backend returns `nan`; the
  mpmath backend raises `ZeroDivisionError` (MPFR divides by zero to inf/nan,
  mpmath raises).
- `t = 0` gives the same split.
- `t < 0` silently returns e^{+1} = 2.718 for `1/(s + 1)` on both backends.
- The bare `except:` would also catch `KeyboardInterrupt`.
- The loop is unchanged since BASE v1.2.0 (cdb3c36) apart from the coefficient
  source; it follows the Mathematica reference `gwr.m`.

## What was ruled out

The coefficient change in this branch (exact integers, Pascal's rule) does not
affect it: the same split occurs at BASE.

## What would unblock it

Decide the contract: raise `ValueError` for `t <= 0` at the `gwr()` entry, and
make a degenerate Wynn-rho step behave the same on both backends (stop
extrapolation and return the last good estimate, as `gwr.m` does). Narrow the
bare `except:` to the arithmetic errors it is meant to catch.

## Files

- `gwr_inversion/_gwr.py`: `gwr_single`, `gwr`
