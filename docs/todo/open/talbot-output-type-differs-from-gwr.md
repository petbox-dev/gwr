# talbot() returns mpf and fails on 0-d arrays; gwr() returns float and accepts them

Severity: low
Opened: 2026-10-19

## What happens now

Version 2.0 changed `gwr()` to return `float` / float64 `ndarray` by default
(`dtype="float"`), with `dtype="mpf"` for multiprecision output, and to
normalise its time input once (`_materialize_time`: Iterable to list, length-1
axes squeezed). `talbot()` was not changed: it returns `mp.mpf` (an object
ndarray for array input), has no `dtype` parameter, and iterates `time`
directly.

## Evidence

- Phase 4 review of the `/code-fix` run: `talbot(fn, np.array(2.0))` raises
  `TypeError` (iteration over a 0-d array); `gwr()` returns shape (1,).
  The 0-d failure was already at BASE v1.2.0.
- `talbot(model.laplace("rate"), [1, 10], degree=24)` agrees with `gwr` to
  about 2e-13, so the numerics are fine; only the container and type differ.

## Blast radius

Code that switches between the two inversion methods must convert types; a
0-d array input to `talbot` fails.

## What would unblock it

Decide the `talbot` output policy. The consistent choice is a `dtype` parameter
with the same default as `gwr` (a breaking change for `talbot` callers, so it
belongs in a version note), plus routing `time` through `_materialize_time`.

## Files

- `gwr_inversion/_gwr.py`: `talbot`, `_materialize_time`
- `gwr_inversion/_talbot.py`: `talbot_single`
