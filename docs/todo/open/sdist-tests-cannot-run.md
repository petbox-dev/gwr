# The sdist ships the test modules but not the reference data they load

Severity: low
Opened: 2026-10-21

## What happens now

`MANIFEST.in` includes only `README.md` and `LICENSE`. The setuptools default
still adds `test/test_*.py` to the sdist, but not `test/reference_data.npz`.
Five test modules (`test_ln`, `test_unit_step`, `test_radial_flow`,
`test_input_types`, `test_f0_hard`) call `np.load(.../'reference_data.npz')` at
import, so pytest run from an unpacked sdist fails at collection.

## Evidence

- `/code-fix` run of 2026-10-21: a scratch `python -m build` sdist listed
  `test/test_*.py` and no `.npz`, PNG or notebook files.
- CI never tests from the sdist: `python -m build` uses the sdist only to build
  the wheel, and the smoke test imports the wheel.

## What would unblock it

Decide whether the sdist is meant to carry runnable tests (as distribution
packagers expect). If yes: `include test/reference_data.npz` in `MANIFEST.in`
and add a CI step that unpacks the sdist and runs pytest. If no: `prune test`.

## Files

- `MANIFEST.in`
- `.github/workflows/ci.yml`: `build` job
