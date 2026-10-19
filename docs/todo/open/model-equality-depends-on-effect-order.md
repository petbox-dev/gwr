# Model(g, Skin, WBS) != Model(g, WBS, Skin) though composition order is fixed by level

Severity: low
Opened: 2026-10-19

## What happens now

`Model` is a frozen dataclass whose compared fields are `geometry` and the
`effects` tuple in the order given. The composition order does not depend on
that tuple order (region effects, then geometry, then skin, then storage), so
two models with the same physics can compare unequal and hash differently.

## Evidence

Executed in the `/code-fix` Phase 3 and Phase 4 reviews:

- `Model(g, Skin(1), WellboreStorage(2)) != Model(g, WellboreStorage(2), Skin(1))`.
- `repr(model)` prints `Model(geometry=..., effects=(...))`, which cannot be
  passed back to `__init__` (`*effects` is positional).
- `dataclasses.replace(model, ...)` raises `TypeError` for the same reason.
- Pickle, copy, deepcopy, `==` and `hash` after unpickle are consistent.
- The `level` ClassVar on each effect is not read by any code; the order is
  decided by `isinstance` in `_validate` and `_wellbore_effects`.

## Blast radius

Using models as dict keys or de-duplicating them in a set treats equivalent
models as distinct. No numeric result is affected.

## What would unblock it

Sort `effects` by level and type name in `Model.__init__` (or compare a
canonical form), and either use `level` for that sort or delete it. Optionally
add a `__repr__` that round-trips.

## Files

- `gwr_inversion/models.py`: `Model.__init__`, effect `level` ClassVars
