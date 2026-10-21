# Project TODOs

Work that is not done lives here: this index, plus one document per item in
`open/` holding that item's evidence. Closed items move to `closed/` and are
kept, because the reasoning that settled them is what gets re-litigated.

Rules, in brief:

- Never leave a `# TODO:` in the code. Fix it now, or open a document here.
- Name the file for the problem, not the fix. The name does not change when the
  item closes; only the directory does.
- A document starts with exactly this header, then prose sections (what happens
  now, evidence, what was ruled out, blast radius, what is still unverified,
  what would unblock it, files). Target 60 lines; 120 is the ceiling.

  ```markdown
  # <the problem in one line>

  Severity: high | medium | low
  Opened: YYYY-MM-DD
  ```

- Adding an item means the document AND its row in the open table below. A
  document with no row is invisible.
- Close an item when no person has an action left: move it to `closed/`, move
  its row to `## Done`, add `Closed: <date> -- settled by: <one line>` and
  `If this was wrong, you will see: <symptom>` to the top of the document.

## Open

Sorted by severity, then oldest first.

| Item | Opened | Severity | Summary |
|---|---|---|---|
| [variable-wellbore-storage-sign-domain](open/variable-wellbore-storage-sign-domain.md) | 2026-10-19 | medium | Negative cDd or aD gives negative early pressure or a pole on the GWR contour; definition unverified against Hegeman 1993 |
| [reference-tests-require-gmpy2](open/reference-tests-require-gmpy2.md) | 2026-10-19 | medium | Regression tests pin gmpy2 rounding at ATOL 1e-12 and fail on an install without gmpy2 |
| [fn-second-parameter-always-receives-precin](open/fn-second-parameter-always-receives-precin.md) | 2026-10-21 | medium | A transform's second parameter receives precin even with a default: fn(s, k=1.0) silently returns a wrong value |
| [readme-speed-tables-predate-lean-core](open/readme-speed-tables-predate-lean-core.md) | 2026-10-19 | low | Unit step, Ln(t) and sin(t) timing tables in the README were measured before the 2.0 GWR core |
| [fracture-boundary-figures-show-no-boundary](open/fracture-boundary-figures-show-no-boundary.md) | 2026-10-19 | low | frac_pressure.png and frac_rate.png plot three identical curves: the time range never reaches the boundary |
| [talbot-output-type-differs-from-gwr](open/talbot-output-type-differs-from-gwr.md) | 2026-10-19 | low | talbot() returns mpf and fails on 0-d arrays; gwr() returns float and accepts them |
| [f0-hard-pins-nonconverged-values](open/f0-hard-pins-nonconverged-values.md) | 2026-10-19 | low | test_f0_hard pins wrong values at M=768 and vacuous ones at M=256/512, and never reads sin_expected |
| [gmpy2-backend-converts-time-through-str](open/gmpy2-backend-converts-time-through-str.md) | 2026-10-19 | low | gm_mpf builds mpfr(str(t)): the gmpy2 backend inverts the decimal string, mpmath the binary double |
| [gwr-single-error-handling-differs-by-backend](open/gwr-single-error-handling-differs-by-backend.md) | 2026-10-19 | low | A constant Gaver sequence gives nan on gmpy2 and ZeroDivisionError on mpmath; a bare except hides the class |
| [lap-transform-validation-cost-on-cheap-geometries](open/lap-transform-validation-cost-on-cheap-geometries.md) | 2026-10-19 | low | Direct lap_transform on linear_bounded_nf is 1.0-1.7 us/call (10-15%) slower than v2.0 HEAD c14c218 |
| [parallel-pool-keeps-running-after-a-failed-call](open/parallel-pool-keeps-running-after-a-failed-call.md) | 2026-10-19 | low | A raising chunk leaves the other chunks running in the persistent pool; untested on fork platforms |
| [model-equality-depends-on-effect-order](open/model-equality-depends-on-effect-order.md) | 2026-10-19 | low | Model(g, Skin, WBS) != Model(g, WBS, Skin) though composition order is fixed by level |
| [linear-constant-pressure-accepts-the-boundary](open/linear-constant-pressure-accepts-the-boundary.md) | 2026-10-21 | low | LinearConstantPressure(distance=1) is accepted; pressure is nan and rate raises |

## Done

| Item | Closed | Settled by |
|---|---|---|
