# README timing tables for unit step, Ln(t) and sin(t) predate the 2.0 GWR core

Severity: low
Opened: 2026-10-19

## What happens now

The README's "Benchmarks" section (unit step, Ln(t), sin(t) at t = 1000) was
measured before commit 1013d3a ("Lean GWR core"), which
replaced per-term factorial and binomial arithmetic with exact cached integers.
The "Performance" table near the top was re-measured on 2026-10-19 and is
current.

## Evidence

Measured in the `/code-fix` Phase 3 plumbing audit (one run each, loaded laptop):

| case | README mpmath / gmpy2 (s) | measured mpmath / gmpy2 (s) |
|---|---|---|
| Ln(t), 31 points, M = 32 | 0.582 / 0.126 ("~4.6x") | 0.126 / 0.056 (2.2x) |
| Unit step, M = 128 | 0.242 / 0.026 ("~10x") | 0.0576 / 0.0203 (2.8x) |
| sin(t) at t = 1000, M = 512 | 5.007 / 0.645 | 1.511 / 0.581 (2.6x) |

- The "~10x", "~10-100x", "~4.6x" and "~11 minutes to 6 seconds (~110x)" prose
  claims (README at f9590a2) were removed by the 2.0 README rewrite (d4d3359);
  the tables they summarised are still in the README.
- `test_f0_hard.py` (M = 256-1024) took 10.3 s with gmpy2.
- The README "mpmath" column was measured with gmpy2 installed, so mpmath ran on
  its own gmpy backend (`mpmath.libmp.BACKEND == "gmpy"`). Unit step M = 128 at
  t = 2.0: 0.0625 s that way, 0.135 s with gmpy2 absent (`/code-fix` 2026-10-21,
  one run, loaded laptop). A plain `pip install gwr-inversion` is the slower case.
- With gmpy2 absent, `test_f0_hard.py` took 135 s (M = 1024: 94 s); the README
  sin(t) table gives 653.559 s for mpmath at M = 1024 (same run as above).
- Error columns match the README for unit step (M = 8-128), Ln(t) and sin(t)
  with gmpy2 (M = 32-1024), re-run 2026-10-21.

## What is still unverified

- The sin(t) mpmath rows at M = 512-1024 were not re-run with the README's
  inputs.

The README rewrite of 2026-10-19 kept the unit step, Ln(t) and sin(t) tables
under "Benchmarks against other inversion methods", with a note that their
timings predate 2.0. It dropped the radial-with-besselk table, which the
current Performance table supersedes.

## What would unblock it

Run `python test/compare_methods.py` on an idle machine, once with gmpy2 and
once without, replace the three tables (label which mpmath backend each column
used), and remove the "measured before the 2.0 GWR core" note in the README.
Then correct the f0_hard timing in `CLAUDE.md` (Testing).

## Files

- `README.md`: "Benchmarks" section
- `test/compare_methods.py`
