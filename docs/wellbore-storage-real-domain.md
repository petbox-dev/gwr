# Wellbore storage in the time domain

Status: verified numerically against `gwr()` (table in section 6). Script:
`test/test_wbs_real_domain.py`.

## Summary

The package composes wellbore storage in the Laplace domain as

    p̄_wD(s) = p̄_S(s) / (1 + C_D s² p̄_S(s)),      p̄_S = p̄_D + S / s          (1)

(`_modifiers.apply_wbs`, after `apply_skin`). In the time domain (1) is exactly
the linear Volterra equation

    p_wD(t) + C_D [ S p_wD'(t) + ∫₀ᵗ p_D'(t − τ) p_wD'(τ) dτ ] = p_D(t) + S,   (2)

with p_wD(0) = 0. Equation (2) holds for every geometry; (1) is its resolvent
written symbolically. So "a real-domain solution" means solving (2), and whether a
closed form exists depends only on the kernel p_D'. The answer, geometry by
geometry:

| Kernel p_D(t) | Real-domain p_wD(t) | Form |
|---|---|---|
| Power law a tᵛ (linear ν = 1/2, bilinear ν = 1/4, ...) | (t / C_D) E_{1−ν, 2}(−t^{1−ν} / (a Γ(ν+1) C_D)) | closed (Mittag-Leffler; erfc for linear) |
| Any p̄_S rational in s^{1/m} (power law with skin, variable storage) | finite sum of Mittag-Leffler terms | closed |
| Bounded (no-flow or constant-pressure outer boundary) | PSS line + Σ aₙ exp(−αₙ² t), αₙ from a transcendental eigen-equation | exact series |
| Infinite radial, finite wellbore, S ≥ 0 | (4/π²) ∫₀^∞ (1 − e^{−u²t}) du / (u³ [A(u)² + B(u)²]) | exact real integral |
| Infinite radial, S < 0 | the same integral + a (e^{s* t} − 1), s* > 0 | exact; grows without bound |

What does not exist is an elementary closed form for the infinite radial case. The
reason is structural, not a gap in the literature: the C_D = 0, S = 0 case of the
integral above is the van Everdingen and Hurst (1949) constant-rate cylindrical
source, which has had no elementary form for 77 years, and any closed form with
storage would contain it as the C_D → 0 limit. This is a strong heuristic, not a
theorem.

## 1. Derivation of (2)

Clear the denominator of (1): p̄_wD + C_D (s p̄_S)(s p̄_wD) = p̄_S. The initial-value
theorem on (1) gives p_wD(0⁺) = lim s p̄_wD = lim 1/(C_D s) = 0, and p_D(0⁺) = 0 for
every geometry here, so s p̄_wD = L{p_wD'} and s p̄_S = L{p_D'} + S. The product
inverts to a convolution, which is (2). With the sandface rate q = 1 − C_D p_wD'
(the material balance of a wellbore of storage C_D fed at unit surface rate), (2)
is Duhamel's principle, p_wD = S q + ∫₀ᵗ q(τ) p_D'(t − τ) dτ. Nothing in the
derivation uses the geometry.

## 2. Power-law kernels: Mittag-Leffler closed form

If p_D(t) = a tᵛ (0 < ν < 1) then p̄_D = a Γ(ν+1) s^{−ν−1} and, with S = 0,

    p̄_wD = (1/C_D) s^{−ν−1} / (s^{1−ν} + λ),     λ = 1 / (a Γ(ν+1) C_D).

The standard pair L⁻¹{s^{α−β} / (s^α + λ)} = t^{β−1} E_{α,β}(−λ t^α) (Podlubny
1999, eq. 1.80; Haubold, Mathai and Saxena 2011) with α = 1 − ν, β = 2 gives

    p_wD(t) = (t / C_D) E_{1−ν, 2}(−t^{1−ν} / (a Γ(ν+1) C_D)).                    (3)

E_{α,β} is entire, so (3) is a convergent power series in t^{1−ν} for all t: the
early-time limit t/C_D is its first term, and the large-argument expansion
E_{α,β}(−x) ~ −Σₖ (−x)^{−k} / Γ(β − αk) gives the late-time behaviour term by term:
p_wD = a tᵛ − C_D a² Γ(ν+1)² t^{2ν−1} / Γ(2ν) + O(t^{3ν−2}), the usual storage shift
(−C_D for linear flow).

Linear flow (`linear_inf`, x_D = 0): a = 2/√π, ν = 1/2, and E_{1/2,2} reduces to
the error function:

    p_wD(t) = C_D [ exp(t / C_D²) erfc(√t / C_D) − 1 ] + 2 √(t / π).              (4)

Check: expanding at small t gives t/C_D + O(t^{3/2}); at large t, 2√(t/π) − C_D +
C_D² / √(π t) + ....

Skin, variable storage, or a bounded linear geometry in its early-time regime make
p̄_S rational in √s (or in s^{1/4} for bilinear flow). A rational function of
z = s^{1/m} splits by partial fractions into terms 1/zᵏ and 1/(z − zⱼ), each of
which is a Mittag-Leffler term; for m = 2, L⁻¹{1/(√s + b)} = 1/√(πt) −
b exp(b²t) erfc(b√t), valid for complex b. So the whole power-law family with skin
and (variable) storage is closed-form. The erfc expression for linear flow with
skin is not written out here; `gwr()` and the Volterra scheme of section 5 agree
to 1e-5 or better on that case, which is the scheme's discretization error.

## 3. Bounded geometries: exponential eigen-series

For a bounded reservoir p̄_wD is meromorphic in s (the log branch of K₀ cancels
between the K and I terms), so the Bromwich integral is a sum of residues.

Linear, no-flow (`linear_bounded_nf`, x_D = 0): p̄_D = coth(√s) / (s √s), and
p̄_wD = 1 / (s [√s tanh √s + C_D s]). Poles at s = −αₙ² where, with √s = iα,

    tan αₙ = −C_D αₙ,     αₙ ∈ ((n − ½)π, nπ),   n = 1, 2, ...

and a double pole at s = 0. The residues evaluate in closed form:

    p_wD(t) = t / (1 + C_D) + 1 / (3 (1 + C_D)²)
              − 2 Σₙ exp(−αₙ² t) / (αₙ² (1 + C_D + C_D² αₙ²)).                      (5)

C_D = 0 gives αₙ = nπ and the Carslaw and Jaeger closed-slab result
t + 1/3 − (2/π²) Σ exp(−n²π²t)/n². The pseudo-steady slope 1/(1 + C_D) is the
reservoir plus wellbore storage in series. At t = 0 the series sums to the
constant, so p_wD(0) = 0 (checked: 2 Σ 1/(n²π²) = 1/3).

Radial, no-flow at r_eD = R (`radial_bounded_nf`) with skin: continuing to
√s = iα with K_ν(iα) = (π/2)(−i)^{ν+1} H_ν⁽²⁾(α), I_ν(iα) = iᵛ J_ν(α) gives

    p̄_S(−α²) = −Φ(α) / (α³ Ψ(α)) − S / α²,
    Φ(α) = J₀(α) Y₁(αR) − Y₀(α) J₁(αR),    Ψ(α) = J₁(α) Y₁(αR) − Y₁(α) J₁(αR),

and the eigen-equation 1 + C_D s² p̄_S = 0 becomes

    C_D α Φ(α) = (1 − C_D S α²) Ψ(α).                                             (6)

C_D = 0 recovers Ψ(αₙ) = 0, the van Everdingen and Hurst eigenvalues. The s = 0
double pole contributes t/H(0) − H'(0)/H(0)² with H(0) = (R² − 1)/2 + C_D (the
storage again adds to the pore volume), and each αₙ contributes
exp(−αₙ² t) / G'(−αₙ²) with G = 1/p̄_S + C_D s². The test evaluates H'(0) and G'
numerically; both have closed forms in Bessel functions that are not worth the
page. Constant-pressure boundaries replace Y₁(αR), J₁(αR) by Y₀(αR), J₀(αR) and
drop the s = 0 pole.

## 4. Infinite radial: the branch-cut integral

`radial_inf` with skin and storage,

    p̄_wD = [K₀(√s) + S √s K₁(√s)] / ( s [ √s K₁(√s) + C_D s (K₀(√s) + S √s K₁(√s)) ] ),

has a branch cut on the negative real axis and no other singularity when S ≥ 0.
Collapsing the Bromwich contour onto the cut with s = u² e^{±iπ} and subtracting
p_wD(0) = 0 gives

    p_wD(t) = (4/π²) ∫₀^∞ (1 − e^{−u²t}) du / ( u³ [A(u)² + B(u)²] ),
    A(u) = u C_D J₀(u) − (1 − C_D S u²) J₁(u),
    B(u) = u C_D Y₀(u) − (1 − C_D S u²) Y₁(u).                                     (7)

This is the form given by Agarwal, Al-Hussainy and Ramey (1970) for the problem
they posed; I re-derived it rather than copying it, since the paper was not at
hand, and the test pins it to `gwr()` at 1e-14. C_D = S = 0 is the van Everdingen
and Hurst cylindrical-source integral. The integrand is smooth and positive; the
difficulty is numerical only. The mass sits at u ~ 1/(C_D ...) for large C_D while
the (1 − e^{−u²t}) factor turns on at u ~ 1/√t, so a single adaptive quadrature on
[0, ∞) misses the peak (observed: `scipy.integrate.quad` returned 0 ± 2e-12 for
C_D = 10⁴, S = 10, t = 0.01 where the value is 1e-6). Sixty log-spaced panels from
1e-8 to 10/√t plus the tail recover 1e-13 relative error at about 0.3 s per time
point, against about 10 ms for `gwr()` at M = 32. The integral is the real-domain
solution; it is not the cheaper way to compute it.

## 5. Any geometry: the Volterra equation solved directly

Equation (2) can be marched in time without a transform. With a piecewise-constant
sandface rate per step, product integration of the exact kernel gives, at step k,

    p_k = (S + Δ_kk) q_k + Σ_{j<k} q_j Δ_kj,   Δ_kj = p_D(t_k − t_{j−1}) − p_D(t_k − t_j),
    q_k = 1 − C_D (p_k − p_{k−1}) / (t_k − t_{k−1}),

one linear equation per step. On the linear kernel it converges to (4) at first
order (relative error 1.5e-3 at 100 log-spaced steps over seven decades, 4e-4
limited by the first step). Two caveats:

- It needs p_D(t) itself in the time domain, which for the cylindrical source is
  the integral (7) at C_D = S = 0, so for radial geometries this route is not
  transform-free in any useful sense.
- It is ill-posed when the kernel vanishes at the origin. The line source
  p_D = ½ E₁(1/(4t)) has p_D'(0⁺) = 0, so for steps shorter than the time the line
  source takes to respond, q_k drops out of its own equation and the scheme
  differentiates the history with gain C_D/Δt: it blew up at Δt = 0.01, C_D = 10³.
  The finite-wellbore kernel p_D' ~ 1/√(πt) is an Abel kernel and the scheme is
  well posed.

## 6. Numerical evidence

All against `Model(...).pressure(t, M=32)` (gmpy2 backend, 30-digit mpmath for the
closed forms). Relative error:

| Case | Parameters | t | rel. error |
|---|---|---|---|
| (4) linear, erfc | C_D = 0.1, 10, 10³ | 1e-3 ... 1e5 | 0 to 2e-16 |
| (3) power law, Mittag-Leffler | ν = 1/2, 1/4, 3/4; C_D = 5 | 1e-2 ... 1e3 | 0 |
| (5) linear no-flow series | C_D = 0, 0.5, 20 | 1e-3 ... 10 | 0 |
| (6) radial no-flow series, 40 modes | C_D = 5, S = 2, R = 10 | 0.1 ... 1e4 | 3e-11 to 3e-15 |
| (6) radial no-flow series, 40 modes | C_D = 50, S = 1, R = 50 | 0.1 / 1 / 100 | 7e-5 / 6e-9 / 1e-11 (truncation) |
| (7) radial infinite integral | C_D = 0..10⁴, S = 0..10 | 1e-2 ... 1e6 | ≤ 1e-14 |
| (7) + pole, S < 0 | C_D = 10, S = −2 | 0.01 / 1 / 10 | 0 / 1e-16 / 0 |
| (7) + pole, S < 0 | C_D = 10, S = −2 | 100 | 3e-10 (`gwr` at 9.07e5) |

## 7. Negative skin: the composition has a growing mode

For every S < 0 tested, 1 + C_D s² p̄_S(s) has a simple real zero s* > 0 (none off
the axis were found; the residue sum closes exactly with the real pole alone). The
exact inverse of (1) is then (7) plus a (e^{s* t} − 1), a = p̄_S(s*) / F'(s*), and it
grows without bound:

| C_D | S | s* | a | exact at t = 100 |
|---|---|---|---|---|
| 10 | −2 | 0.1453 | 0.4456 | 9.07e5 |
| 100 | −1 | 0.3841 | 1.97e-3 | 9.4e13 |
| 10⁴ | −3 | 3.45e-3 | 1.63e-3 | 1.005e-2 (growth starts near t ~ 300) |
| 100 | −0.01 | 9902 | 2.0e-10 | 2e33 at t = 0.01 |

Mechanism: the infinitesimal skin p_wD − p_sf = S q with S < 0 is a negative
resistance between wellbore and sandface. The storage balance then has positive
feedback, C_D dp_wD/dt = 1 + (p_wD − p_sf)/|S|, which the reservoir only overcomes
once p_D(t) exceeds |S|; the growth rate s* follows K₀/(√s K₁) = |S| − 1/(C_D s),
about 1/S² for small |S|. So the model is unstable for every S < 0, and for small
|S| the instability is at times no test reaches but the amplitude blows up
immediately. What `gwr()` returns in that regime is the branch-cut part (it
matches (7) to 1e-14 for C_D = 100, S = −0.01), that is, the stable curve the
model does not actually produce; at larger |S| it returns the growing solution up
to the point where the Gaver sequence fails, then garbage (0.43 at t = 100 for
C_D = 100, S = −1, where the exact value is 9.4e13). The well-test practice of
using the effective wellbore radius r_w e^{−S} and C_D e^{2S} for stimulated wells
rather than the S/s skin term avoids this. Whether the 1970 paper discusses the
instability I could not check. Tracked in
`docs/todo/open/negative-skin-with-storage-has-a-growing-mode.md`.

## 8. What this means for the package

- Keep GWR as the general evaluator. The exact forms above are references, and
  the test file uses them that way: `gwr()` plus the compositor now has an
  independent analytic check on four geometries, which it did not have before
  (the reference data in `test/reference_data.npz` was generated by `gwr()`).
- For linear and other power-law regimes (4) and (3) are cheaper than inversion
  and exact; a user who only needs those can call them directly.
- The negative-skin instability is a modelling limit of (1), not of the
  inversion. Rejecting or warning on `Skin(S < 0)` combined with storage is the
  todo item.

## References

- Agarwal, R. G., Al-Hussainy, R. and Ramey, H. J. Jr. 1970. An Investigation of
  Wellbore Storage and Skin Effect in Unsteady Liquid Flow: I. Analytical
  Treatment. SPE Journal 10 (3): 279–290. SPE-2466-PA.
  https://doi.org/10.2118/2466-PA
- van Everdingen, A. F. and Hurst, W. 1949. The Application of the Laplace
  Transformation to Flow Problems in Reservoirs. JPT 1 (12): 305–324.
  https://doi.org/10.2118/949305-G
- Carslaw, H. S. and Jaeger, J. C. 1959. Conduction of Heat in Solids, 2nd ed.
  Oxford. Section 3.8 (slab with prescribed flux).
- Podlubny, I. 1999. Fractional Differential Equations. Academic Press. Section
  1.2.2, eq. (1.80).
- Haubold, H. J., Mathai, A. M. and Saxena, R. K. 2011. Mittag-Leffler Functions
  and Their Applications. J. Appl. Math. 2011: 298628.
  https://doi.org/10.1155/2011/298628
- Valkó, P. P. and Abate, J. 2002. Comparison of Sequence Accelerators for the
  Gaver Method of Numerical Laplace Transform Inversion. Comput. Math. Appl.
  48: 629–636. https://doi.org/10.1016/j.camwa.2002.10.017
