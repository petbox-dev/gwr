# LinearConstantPressure accepts distance = 1, where every output is nan or raises

Severity: low
Opened: 2026-10-21

## What happens now

`LinearConstantPressure.__post_init__` (`gwr_inversion/models.py`) accepts
`distance` in `[0, 1]`. At `x_D = 1` the transform `linear_bounded_cp`
(`gwr_inversion/transforms.py`) is `sinh(sqrt(s_eff) (1 - x_D)) / ...`, which is
identically 0: the point is the constant-pressure boundary itself.

## Evidence

Executed in the `/code-fix` run of 2026-10-21, M = 16:

- `Model(LinearConstantPressure(distance=1.0)).pressure([1.0, 10.0])` returns
  `[nan, nan]`.
- `.rate([1.0])` raises `ZeroDivisionError`.

`docs/api.md` now states that 1 gives `nan`. The `models.py` and
`transforms.py` docstrings still give `0 <= x_D <= 1` without the caveat.

## What would unblock it

Decide whether `distance = 1` is an input error (reject it: domain `[0, 1)`) or
a valid query whose answer is the boundary pressure (return 0 for pressure and
document rate as undefined). Then align `models.py`, `transforms.py`,
`docs/api.md` and add a test.

## Files

- `gwr_inversion/models.py`: `LinearConstantPressure`
- `gwr_inversion/transforms.py`: `linear_bounded_cp`
- `docs/api.md`: Geometries table
