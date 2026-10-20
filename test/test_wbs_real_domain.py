"""
Wellbore storage: exact time-domain solutions against the Laplace-domain compositor.

Derivations and evidence: docs/wellbore-storage-real-domain.md. Each test pins
``gwr()`` on a composed model to a closed form, an eigen-series or a real integral
that does not go through a Laplace inversion.
"""
from __future__ import annotations

import importlib.util
from typing import Any, Callable, Iterator, List, Tuple

import numpy as np
import pytest
from mpmath import mp

from gwr_inversion import (
    LinearInfinite, LinearNoFlow, Model, RadialInfinite, RadialNoFlow, Skin,
    WellboreStorage, besselk, gwr,
)

HAS_SCIPY = importlib.util.find_spec("scipy") is not None
requires_scipy = pytest.mark.skipif(not HAS_SCIPY, reason="scipy not installed")


@pytest.fixture(autouse=True)
def precision() -> Iterator[None]:
    with mp.workdps(30):
        yield


def rel(a: Any, b: Any) -> float:
    return abs(float(a) - float(b)) / abs(float(b))


def model(geometry: Any, skin: float, storage: float) -> Model:
    effects: List[Any] = []
    if skin:
        effects.append(Skin(skin))
    if storage:
        effects.append(WellboreStorage(storage))
    return Model(geometry, *effects)


# ---------------------------------------------------------------- power-law kernels

def mittag_leffler(alpha: Any, beta: Any, z: Any) -> Any:
    """E_{alpha,beta}(z) by its power series; extra digits absorb the cancellation."""
    z = mp.mpf(z)
    with mp.workdps(max(mp.dps, int(2 * abs(z)) + 40)):
        total = mp.mpf(0)
        power = mp.mpf(1)
        for k in range(20000):
            term = power / mp.gamma(alpha * k + beta)
            total += term
            power *= z
            if k > 10 and abs(term) < mp.eps * abs(total) and abs(power) < abs(term):
                break
        return total


def power_law_storage(t: float, storage: float, nu: Any, a: Any) -> Any:
    """p_D = a t^nu with storage:  (t/C) E_{1-nu,2}(-t^{1-nu} / (a Gamma(nu+1) C))."""
    tt, c = mp.mpf(t), mp.mpf(storage)
    lam = 1 / (a * mp.gamma(nu + 1) * c)
    return tt / c * mittag_leffler(1 - nu, 2, -lam * tt ** (1 - nu))


def linear_storage(t: float, storage: float) -> Any:
    """Linear flow with storage:  C [e^{t/C^2} erfc(sqrt t / C) - 1] + 2 sqrt(t/pi)."""
    tt, c = mp.mpf(t), mp.mpf(storage)
    x = mp.sqrt(tt) / c
    return c * (mp.exp(x * x) * mp.erfc(x) - 1) + 2 * mp.sqrt(tt / mp.pi)


@pytest.mark.parametrize("storage", [0.1, 10.0, 1000.0])
def test_linear_infinite_storage_is_erfc(storage: float) -> None:
    times = [1e-3, 1e-1, 10.0, 1e3, 1e5]
    computed = Model(LinearInfinite(), WellboreStorage(storage)).pressure(times, M=32)
    for t, value in zip(times, computed):
        assert rel(linear_storage(t, storage), value) < 1e-13


def test_erfc_form_is_the_mittag_leffler_form() -> None:
    # linear flow: a = 2/sqrt(pi), nu = 1/2, so E_{1/2,2} must reduce to erfc
    t, storage = 3.0, 7.0
    assert rel(power_law_storage(t, storage, mp.mpf(1) / 2, 2 / mp.sqrt(mp.pi)),
               linear_storage(t, storage)) < 1e-25


@pytest.mark.parametrize("nu", [mp.mpf(1) / 4, mp.mpf(3) / 4])
def test_power_law_kernel_storage_is_mittag_leffler(nu: Any) -> None:
    storage, a = 5.0, mp.mpf(1)

    def fn(s: Any) -> Any:
        p_d = a * mp.gamma(nu + 1) * s ** (-nu - 1)   # L{a t^nu}
        return 1 / (1 / p_d + storage * s ** 2)

    times = [1e-2, 1.0, 1e2, 1e3]
    computed = gwr(fn, times, M=32)
    for t, value in zip(times, computed):
        assert rel(power_law_storage(t, storage, nu, a), value) < 1e-13


# ---------------------------------------------------------------- bounded: eigen-series

def linear_noflow_storage(t: float, storage: float, nterms: int = 200) -> Any:
    """t/(1+C) + 1/(3(1+C)^2) - 2 sum e^{-a_n^2 t} / (a_n^2 (1 + C + C^2 a_n^2)),
    tan a_n = -C a_n."""
    tt, c = mp.mpf(t), mp.mpf(storage)
    total = tt / (1 + c) + 1 / (3 * (1 + c) ** 2)
    for n in range(1, nterms + 1):
        lo, hi = (n - mp.mpf(1) / 2) * mp.pi, n * mp.pi
        if c > 0:
            bracket = (lo + mp.mpf("1e-9"), hi - mp.mpf("1e-9"))
            root = mp.findroot(lambda x: mp.tan(x) + c * x, bracket, solver="bisect", tol=1e-28)
        else:
            root = hi
        total -= 2 * mp.exp(-root * root * tt) / (root * root * (1 + c + c * c * root * root))
    return total


@pytest.mark.parametrize("storage", [0.0, 0.5, 20.0])
def test_linear_noflow_storage_is_exponential_series(storage: float) -> None:
    times = [1e-3, 1e-2, 1e-1, 1.0, 10.0]
    computed = model(LinearNoFlow(), 0.0, storage).pressure(times, M=32)
    for t, value in zip(times, computed):
        assert rel(linear_noflow_storage(t, storage), value) < 1e-13


def phi_psi(alpha: Any, r: Any) -> Tuple[Any, Any]:
    """Phi = J0(a) Y1(aR) - Y0(a) J1(aR),  Psi = J1(a) Y1(aR) - Y1(a) J1(aR)."""
    j0, j1, y0, y1 = (mp.besselj(0, alpha), mp.besselj(1, alpha),
                      mp.bessely(0, alpha), mp.bessely(1, alpha))
    j1r, y1r = mp.besselj(1, alpha * r), mp.bessely(1, alpha * r)
    return j0 * y1r - y0 * j1r, j1 * y1r - y1 * j1r


def radial_noflow_eigen_series(t: float, storage: float, skin: float, radius: float,
                               alphas: List[Any]) -> Any:
    """PSS line + sum of residues at s = -alpha_n^2, with G = 1/p_S + C s^2."""
    tt, c, sk, r = mp.mpf(t), mp.mpf(storage), mp.mpf(skin), mp.mpf(radius)

    def p_s(alpha: Any) -> Any:   # p_S continued to s = -alpha^2
        phi, psi = phi_psi(alpha, r)
        return -phi / (alpha ** 3 * psi) - sk / alpha ** 2

    def g_of_s(s: Any) -> Any:
        z = mp.sqrt(s)
        k0, k1, i0, i1 = mp.besselk(0, z), mp.besselk(1, z), mp.besseli(0, z), mp.besseli(1, z)
        k1r, i1r = mp.besselk(1, z * r), mp.besseli(1, z * r)
        p_d = (k0 * i1r + i0 * k1r) / (s * z * (k1 * i1r - i1 * k1r))
        return 1 / (p_d + sk / s) + c * s * s

    with mp.workdps(60):
        h0 = (r * r - 1) / 2 + c
        eps = mp.mpf("1e-15")
        h1 = (g_of_s(eps) / eps ** 2 - h0) / eps
    total = tt / h0 - h1 / h0 ** 2
    for alpha in alphas:
        g_alpha = lambda x: 1 / p_s(x) + c * x ** 4   # noqa: E731
        total += mp.exp(-alpha * alpha * tt) / (mp.diff(g_alpha, alpha) / (-2 * alpha))
    return total


def radial_noflow_eigenvalues(storage: float, skin: float, radius: float, nterms: int) -> List[Any]:
    """Roots of  C a Phi(a) = (1 - C S a^2) Psi(a)  (Psi = 0 when C = 0)."""
    from scipy import special

    def char(a: Any) -> Any:
        phi = special.j0(a) * special.y1(a * radius) - special.y0(a) * special.j1(a * radius)
        psi = special.j1(a) * special.y1(a * radius) - special.y1(a) * special.j1(a * radius)
        return storage * a * phi - (1 - storage * skin * a * a) * psi

    def char_mp(a: Any) -> Any:
        phi, psi = phi_psi(a, mp.mpf(radius))
        return storage * a * phi - (1 - storage * skin * a * a) * psi

    grid = np.linspace(1e-4, (nterms + 2) * np.pi / (radius - 1), 400 * (nterms + 2))
    vals = char(grid)
    roots: List[Any] = []
    for lo, hi, flo, fhi in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
        if flo * fhi < 0:
            bracket = (mp.mpf(lo), mp.mpf(hi))
            roots.append(mp.findroot(char_mp, bracket, solver="anderson", tol=1e-26))
        if len(roots) >= nterms:
            break
    return roots


@requires_scipy
@pytest.mark.parametrize("storage, skin", [(0.0, 0.0), (5.0, 2.0)])
def test_radial_noflow_storage_is_exponential_series(storage: float, skin: float) -> None:
    radius = 10.0
    alphas = radial_noflow_eigenvalues(storage, skin, radius, nterms=40)
    times = [1.0, 10.0, 100.0, 1e4]
    computed = model(RadialNoFlow(radius), skin, storage).pressure(times, M=32)
    for t, value in zip(times, computed):
        assert rel(radial_noflow_eigen_series(t, storage, skin, radius, alphas), value) < 1e-10


# ---------------------------------------------------------------- infinite radial: branch cut

def radial_branch_cut_integral(t: float, storage: float, skin: float) -> Tuple[float, float]:
    """(4/pi^2) Int_0^inf (1 - e^{-u^2 t}) / (u^3 [A^2 + B^2]) du; value and quadrature error."""
    from scipy import integrate, special

    def integrand(u: float) -> float:
        j0, j1, y0, y1 = special.j0(u), special.j1(u), special.y0(u), special.y1(u)
        a = u * storage * j0 - (1 - storage * skin * u * u) * j1
        b = u * storage * y0 - (1 - storage * skin * u * u) * y1
        return float(-np.expm1(-u * u * t) / (u ** 3 * (a * a + b * b)))

    # log-spaced panels: the mass can sit at u << 1 while (1 - e^{-u^2 t}) turns on at 1/sqrt(t)
    edges = [0.0] + list(np.logspace(-8, np.log10(10.0 / np.sqrt(t)), 60))
    total = err = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        v, e = integrate.quad(integrand, lo, hi, limit=400, epsabs=0, epsrel=1e-13)
        total += v
        err += e
    v, e = integrate.quad(integrand, edges[-1], np.inf, limit=4000, epsabs=0, epsrel=1e-13)
    return 4 / np.pi ** 2 * (total + v), 4 / np.pi ** 2 * (err + e)


@requires_scipy
@pytest.mark.parametrize("storage, skin", [(0.0, 0.0), (100.0, 0.0), (100.0, 5.0), (1e4, 10.0)])
def test_radial_infinite_storage_is_branch_cut_integral(storage: float, skin: float) -> None:
    times = [1e-2, 1.0, 1e2, 1e4]
    computed = model(RadialInfinite(), skin, storage).pressure(times, M=32)
    for t, value in zip(times, computed):
        integral, err = radial_branch_cut_integral(t, storage, skin)
        assert err < 1e-11 * abs(integral)
        assert rel(integral, value) < 1e-11


def storage_characteristic(storage: float, skin: float) -> Callable[[Any], Any]:
    """F(s) = 1 + C s^2 p_S(s) for radial_inf; its zeros are the poles of p_wD."""
    def p_s(s: Any) -> Any:
        z = mp.sqrt(s)
        return besselk(0, z) / (z * besselk(1, z)) / s + skin / s

    def f(s: Any) -> Any:
        return 1 + storage * s * s * p_s(s)

    f.p_s = p_s   # type: ignore[attr-defined]
    return f


@requires_scipy
def test_negative_skin_with_storage_has_a_real_positive_pole() -> None:
    storage, skin = 10.0, -2.0
    f = storage_characteristic(storage, skin)
    pole = mp.findroot(f, (mp.mpf("0.1"), mp.mpf("0.2")), solver="bisect", tol=1e-25)
    amplitude = f.p_s(pole) / mp.diff(f, pole)   # type: ignore[attr-defined]
    assert rel(pole, 0.145259977) < 1e-8
    assert rel(amplitude, 0.4456344762) < 1e-8
    times = [1e-2, 1.0, 10.0]
    computed = model(RadialInfinite(), skin, storage).pressure(times, M=32)
    for t, value in zip(times, computed):
        integral, _ = radial_branch_cut_integral(t, storage, skin)
        # the integral alone fixes p_wD(0) = 0, so the pole enters as a (e^{s* t} - 1)
        exact = integral + float(amplitude * (mp.exp(pole * t) - 1))
        assert rel(exact, value) < 1e-11
    # the exact inverse grows without bound: 9.07e5 at t = 100
    assert float(amplitude * (mp.exp(pole * 100) - 1)) > 9e5
