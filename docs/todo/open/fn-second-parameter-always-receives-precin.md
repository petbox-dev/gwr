# A transform's second parameter always receives precin, even when it has a default

Severity: medium
Opened: 2026-10-21

## What happens now

`gwr_single` (`gwr_inversion/_gwr.py`) counts `signature(fn).parameters`. One
parameter: `fn(s)`. Two: `fn(s, precin)`. Three or more: `TypeError`. A default
value does not change the count, so a second parameter meant for something
else silently receives the working precision in digits.

## Evidence

Executed in the `/code-fix` run of 2026-10-21 (gmpy2 backend, M = 32, t = 1):

| fn | result | expected |
|---|---|---|
| `def f(s, k=1.0): return 1/(s+k)` | 2.33e-18 | 0.36788 (e^-1) |
| `functools.partial(g, k=1.0)` with `def g(s, k)` | `TypeError: got multiple values for argument 'k'` | 0.36788 |

The first is a silent wrong answer: `k` receives `precin = 67`.

The rule is now documented in `docs/api.md` (`gwr`, `fn` row); the code is
unchanged and the behaviour is the same at f9590a2.

## What would unblock it

Decide the contract. Options: pass `precin` only when the second parameter is
named `precin` (or is keyword-only `precin`); or count only parameters without
defaults; or keep the positional rule and raise when the second parameter has a
default. Any change needs a test per case in the table.

## Files

- `gwr_inversion/_gwr.py`: `gwr_single` (signature dispatch)
- `docs/api.md`: `gwr` parameter table
