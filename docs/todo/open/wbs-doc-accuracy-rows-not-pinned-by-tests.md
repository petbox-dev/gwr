# The wellbore-storage doc's accuracy table cites cases its tests do not pin

Severity: low
Opened: 2026-10-31

## What happens now

`docs/wellbore-storage-real-domain.md` section 6 lists, per closed form, the
parameter and time ranges over which it matches `gwr()`. Some rows go beyond
what `test/test_wbs_real_domain.py` asserts:

- Row (7), radial infinite: "C_D = 0..10⁴, S = 0..10, t to 1e6". The test covers
  (C_D, S) = (0, 0), (100, 0), (100, 5), (1e4, 10) at t <= 1e4.
- Row (6), bounded no-flow eigen-series: cites t = 0.1 and C_D = 50, S = 1,
  R = 50; neither is a test case.
- Section 2 gives the Volterra scheme 1e-5 against `gwr()` on linear flow with
  skin, and section 5 gives 1.5e-3 at 100 steps on the linear kernel, without
  saying they are different step counts. Neither is pinned.

## Evidence

The /code-fix review of 2026-10-31 re-ran the (6) numbers (2.7e-11 to 2.9e-15;
7.3e-5 / 6.3e-9 / 1.4e-11) and they reproduce, so nothing is known to be wrong.
They come from an unpinned run, and the next edit to the closed forms can change
them silently.

## What would unblock it

Either add the cited cases to the test parametrizations (the (7) cases need the
60-panel integral, 2-4 ms each) or trim each row to what the tests pin.

## Files

- `docs/wellbore-storage-real-domain.md` sections 2, 5, 6
- `test/test_wbs_real_domain.py`
