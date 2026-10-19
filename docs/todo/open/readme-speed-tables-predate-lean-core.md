# README timing tables for unit step, Ln(t) and sin(t) predate the 2.0 GWR core

Severity: low
Opened: 2026-10-19

## What happens now

The README's "Benchmarks" section (unit step, Ln(t), radial flow with WBS,
sin(t) at t = 1000) was measured before commit 1013d3a ("Lean GWR core"), which
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

- The "~10x" and "~10-100x" gmpy2 claims contradict the measured 2.8x at M = 32
  in the current Performance table.
- `test_f0_hard.py` (M = 256-1024) took 10.3 s with gmpy2.

## What is still unverified

- The "~11 minutes to 6 seconds (~110x)" claim for M = 1024 on mpmath was not
  re-run; the mpmath number in particular may now be much lower.
- Error columns (accuracy) are expected unchanged, since the algorithm is the
  same; not re-checked.

The README rewrite of 2026-10-19 kept the unit step, Ln(t) and sin(t) tables
under "Benchmarks against other inversion methods", with a note that their
timings predate 2.0. It dropped the radial-with-besselk table, which the
current Performance table supersedes.

## What would unblock it

Run `python test/compare_methods.py` on an idle machine, replace the three
tables, and remove the "measured before the 2.0 GWR core" note in the README.

## Files

- `README.md`: "Benchmarks" section
- `test/compare_methods.py`
