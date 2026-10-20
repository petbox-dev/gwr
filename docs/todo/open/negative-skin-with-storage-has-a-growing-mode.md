# Skin < 0 with wellbore storage puts a pole on the positive real axis; gwr returns the wrong inverse silently

Severity: medium
Opened: 2026-10-20

## What happens now

`WellboreEffects.apply` composes `p_wD = p_S / (1 + C_D s^2 p_S)` with
`p_S = p_D + S/s`. Neither `Skin` nor `lap_transform` restricts the sign of S.
For every S < 0 tried, `1 + C_D s^2 p_S(s)` has a simple real zero s* > 0 for
`radial_inf`, so the exact inverse contains `a (e^{s* t} - 1)` and grows without
bound. `gwr()` neither detects this nor returns it consistently.

## Evidence

Measured in `test/test_wbs_real_domain.py` and
`docs/wellbore-storage-real-domain.md` section 7 (`Model(RadialInfinite(), Skin(S),
WellboreStorage(C_D)).pressure(t, M=32)`, gmpy2 backend). The branch-cut
integral plus the pole term reproduces `gwr()` to 1e-16 where `gwr()` has not
failed:

| C_D | S | s* | a | exact p_wD(100) | gwr p_wD(100) |
|---|---|---|---|---|---|
| 10 | -2 | 0.1453 | 0.4456 | 9.07e5 | 9.07e5 (3e-10 rel) |
| 100 | -1 | 0.3841 | 1.97e-3 | 9.4e13 | 0.43 |
| 1e4 | -3 | 3.45e-3 | 1.63e-3 | 1.005e-2 | 1.005e-2 |
| 100 | -0.01 | 9902 | 2.0e-10 | 2e33 at t = 0.01 | 9.99e-5 at t = 0.01 |

- For small |S| the pole sits at s* ~ 1/S^2 with a tiny residue; `gwr()` never
  samples that far out and returns the branch-cut part alone (it matches the
  C_D = 100, S = -0.01 integral to 1e-14). That curve looks physical. It is not
  the inverse of the function passed in.
- For larger |S| `gwr()` follows the growing solution until the Gaver sequence
  breaks, then returns O(1) values (0.43 at t = 100 for C_D = 100, S = -1).
- No zeros off the real axis were found (grid of complex starting points for
  `findroot`, C_D = 10, S = -2); the residue sum closes with the real pole alone.

Mechanism: `p_wD - p_sf = S q` with S < 0 is a negative resistance in series with
the storage, so `C_D dp_wD/dt = 1 + (p_wD - p_sf)/|S|` has positive feedback until
the reservoir pressure drop exceeds |S|. The model, not the inversion, is at
fault. Standard practice for stimulated wells is the effective wellbore radius
`r_w e^{-S}` with S = 0 (the `C_D e^{2S}` type-curve parameter), which has no such
pole.

## What was ruled out

- Quadrature error in the integral: it matches `gwr()` to 1e-14 for every S >= 0
  case and for S < 0 after adding the pole term.
- A missing complex pole pair: searched; the real pole alone closes the sum.

## Blast radius

Any `Model` or `lap_transform` call with `skin < 0` and `cD > 0`, and any user
transform of the same shape. Bounded geometries have the same mechanism (the
eigen-equation `C a Phi = (1 - C S a^2) Psi` acquires a root at s > 0); not
measured. `test/generate_figures.py` and the README examples use S >= 0.

## What is still unverified

- Whether Agarwal, Al-Hussainy and Ramey (1970) discuss this instability; the
  paper was not available.
- Whether every S < 0 gives a pole for the bounded and fractured geometries, or
  only below a threshold depending on C_D and r_eD.
- The behaviour under `VariableWellboreStorage` (see
  `variable-wellbore-storage-sign-domain.md`, same family of problem).

## What would unblock it

Decide between (a) rejecting `skin < 0` whenever storage is applied, (b) warning
and documenting that the result is the branch-cut part only, or (c) adding an
effective-radius form of skin (`r_wD = e^{-S}` scaling of `t_D` and `C_D`) as
the supported way to model stimulated wells with storage. (c) is what the
literature does; (a) is a two-line check in `WellboreEffects.apply` and
`Model.__init__`.

## Files

- `gwr_inversion/_modifiers.py`: `apply_skin`, `apply_wbs`, `WellboreEffects.apply`
- `gwr_inversion/models.py`: `Skin`, `Model.__init__`
- `test/test_wbs_real_domain.py`: `test_negative_skin_with_storage_has_a_real_positive_pole`
- `docs/wellbore-storage-real-domain.md`: section 7
