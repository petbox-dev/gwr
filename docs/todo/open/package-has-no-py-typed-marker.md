# The package is fully annotated but ships no py.typed marker

Severity: low
Opened: 2026-10-21

## What happens now

`gwr_inversion` passes `mypy --strict` (the `lint` CI job), but the wheel has
no `gwr_inversion/py.typed`. Under PEP 561 a type checker treats the installed
package as untyped, so a user's `mypy` sees every `gwr_inversion` name as
`Any` (or reports a missing-stubs error).

## Evidence

- `ls gwr_inversion/py.typed`: no such file (2026-10-21).
- The wheel built in the `/code-fix` run held only the 13 modules and
  `dist-info`.

## What would unblock it

Add an empty `gwr_inversion/py.typed` and list it as package data
(`[tool.setuptools.package-data] gwr_inversion = ["py.typed"]`). Check the
built wheel contains it. Decide first whether the private-module names the
public API re-exports are stable enough to expose as typed.

## Files

- `gwr_inversion/` (new `py.typed`)
- `pyproject.toml`: `[tool.setuptools]`
