"""
Agreement tests for the duplicated evaluation paths: the ARB and mpf
compositors, the three Bessel backends, and complex ``s`` (Talbot) through
every geometry and every backend combination.
"""
from __future__ import annotations

import itertools
from typing import Any, Callable, Dict, Iterator, Optional, Tuple

import numpy as np
import pytest
from mpmath import mp

import gwr_inversion._flint as _flint
import gwr_inversion._gmpy2_bessel as _gmpy2_bessel
from gwr_inversion import gwr, talbot, transforms

LAPLACE_POINTS = ["1e-6", "0.37", "5.0", "880.0"]

try:
    import scipy  # noqa: F401
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
# hybrid_rate needs scipy (the [hybrid] extra).
requires_scipy = pytest.mark.skipif(not HAS_SCIPY, reason="scipy not installed")
requires_flint = pytest.mark.skipif(not _flint.HAS_FLINT, reason="python-flint not installed")

ARB_GEOMETRIES = [
    (transforms.radial_inf, {}),
    (transforms.linear_inf, {"xD": 0.0}),
    (transforms.linear_inf, {"xD": 0.3}),
    (transforms.linear_bounded_nf, {"xD": 0.2}),
    (transforms.linear_bounded_cp, {"xD": 0.2}),
]
DUAL_POROSITY = [{}, {"lam": 1e-6, "omega": 0.1}]
SKIN = [{}, {"skin": 3.0}, {"skin": -2.0}]
STORAGE = [{}, {"cD": 100.0}, {"cD": 100.0, "cDd": 10.0, "aD": 0.01}]

ALL_GEOMETRIES = [
    (transforms.radial_inf, {}),
    (transforms.radial_bounded_nf, {"reD": 50.0}),
    (transforms.radial_bounded_cp, {"reD": 50.0}),
    (transforms.linear_inf, {"xD": 0.0}),
    (transforms.linear_bounded_nf, {"xD": 0.0}),
    (transforms.linear_bounded_cp, {"xD": 0.0}),
    (transforms.frac_inf, {}),
    (transforms.frac_bounded_nf, {"reD": 10.0}),
    (transforms.frac_bounded_cp, {"reD": 10.0}),
]

# (python-flint on, gmpy2 Bessel on) for each install the dispatch must handle;
# a combination that needs a backend this machine lacks is skipped.
BACKEND_COMBINATIONS = [(True, True), (True, False), (False, True), (False, False)]


@pytest.fixture(autouse=True)
def restore_precision() -> Iterator[None]:
    """Tests here set mp.dps; restore it so later test modules run at their own."""
    with mp.workdps(mp.dps):
        yield


@pytest.fixture
def set_backends(monkeypatch: pytest.MonkeyPatch) -> Callable[[bool, bool], None]:
    """Switch the ARB and gmpy2 Bessel backends off for one test (restored after)."""
    def set_flags(flint_enabled: bool, gmpy2_enabled: bool) -> None:
        if flint_enabled and not _flint.HAS_FLINT:
            pytest.skip("python-flint not installed")
        if gmpy2_enabled and not _gmpy2_bessel.HAS_GMPY2_BESSEL:
            pytest.skip("gmpy2 not installed")
        if not flint_enabled:
            monkeypatch.setattr(_flint, "HAS_FLINT", False)
        if not gmpy2_enabled:
            monkeypatch.setattr(_gmpy2_bessel, "HAS_GMPY2_BESSEL", False)
    return set_flags


@pytest.mark.skipif(not _flint.HAS_FLINT, reason="python-flint not installed")
@pytest.mark.parametrize(
    "geometry, parameters, dual_porosity, skin, storage",
    [(geometry, parameters, dual_porosity, skin, storage)
     for (geometry, parameters), dual_porosity, skin, storage
     in itertools.product(ARB_GEOMETRIES, DUAL_POROSITY, SKIN, STORAGE)])
def test_arb_compositor_matches_mpf(
    monkeypatch: pytest.MonkeyPatch, geometry: Any, parameters: Dict[str, float],
    dual_porosity: Dict[str, float], skin: Dict[str, float], storage: Dict[str, float],
) -> None:
    mp.dps = 34
    keyword_arguments = {**parameters, **dual_porosity, **skin, **storage}
    ball_values = [transforms.lap_transform(mp.mpf(point), geometry, **keyword_arguments)
                   for point in LAPLACE_POINTS]
    monkeypatch.setattr(_flint, "HAS_FLINT", False)
    for point, ball_value in zip(LAPLACE_POINTS, ball_values):
        mpf_value = transforms.lap_transform(mp.mpf(point), geometry, **keyword_arguments)
        assert abs(ball_value / mpf_value - 1) < mp.mpf(10) ** -30, point


BESSEL_ARGUMENTS = ["1e-6", "0.3", "4.7", "40.1", "900.0"]


def _backend_values(order: int, x: Any) -> Tuple[Optional[Any], ...]:
    """K and I of ``order`` at ``x`` from ARB and gmpy2 (None where not installed)."""
    return (_flint.besselk(order, x), _flint.besseli(order, x),
            _gmpy2_bessel.besselk_gmpy2(order, x), _gmpy2_bessel.besseli_gmpy2(order, x))


@pytest.mark.parametrize("digits", [17, 34, 67])
@pytest.mark.parametrize("order", [0, 1])
def test_bessel_backends_agree_with_mpmath(digits: int, order: int) -> None:
    mp.dps = digits
    tolerance = mp.mpf(10) ** -(digits - 2)
    for argument in BESSEL_ARGUMENTS:
        x = mp.mpf(argument)
        k_reference = mp.besselk(order, x)
        i_reference = mp.besseli(order, x)
        arb_k, arb_i, gmpy2_k, gmpy2_i = _backend_values(order, x)
        expected_backends = ([arb_k, arb_i] if _flint.HAS_FLINT else []) + (
            [gmpy2_k, gmpy2_i] if _gmpy2_bessel.HAS_GMPY2_BESSEL else [])
        # An installed backend must evaluate here; None would make the check vacuous.
        assert all(value is not None for value in expected_backends), argument
        for value, reference in ((arb_k, k_reference), (arb_i, i_reference),
                                 (gmpy2_k, k_reference), (gmpy2_i, i_reference)):
            if value is not None:
                assert abs(value / reference - 1) < tolerance, argument


@pytest.mark.parametrize("argument", [mp.mpc(1, 2), mp.mpf(-1), mp.mpf(0), mp.inf])
def test_fast_backends_decline_outside_positive_reals(argument: Any) -> None:
    assert _backend_values(0, argument) == (None, None, None, None)
    assert _flint.besselk0_integral(argument) is None
    assert _flint.besseli0_integral(argument) is None


@pytest.mark.parametrize("flint_enabled, gmpy2_enabled", BACKEND_COMBINATIONS)
@pytest.mark.parametrize("geometry, parameters", ALL_GEOMETRIES)
def test_talbot_complex_s_matches_gwr(
    set_backends: Callable[[bool, bool], None], flint_enabled: bool, gmpy2_enabled: bool,
    geometry: Any, parameters: Dict[str, float],
) -> None:
    # Talbot evaluates the transform on a complex contour; every backend must
    # decline complex s and leave it to mpmath instead of raising.
    set_backends(flint_enabled, gmpy2_enabled)
    function = lambda s: transforms.lap_transform(s, geometry, cD=10.0, **parameters)  # noqa: E731
    by_talbot = float(talbot(function, 2.0, degree=32))
    by_gwr = gwr(function, 2.0, M=32)
    assert by_talbot == pytest.approx(by_gwr, rel=1e-8)


def _radial_bounded_reference(s: Any, outer_radius: float, no_flow: bool) -> Any:
    """Ozkan-Raghavan bounded radial solution from mpmath Bessel functions, 30 digits above."""
    with mp.workdps(mp.dps + 30):
        root = mp.sqrt(s)
        outer = root * outer_radius
        k0, k1 = mp.besselk(0, root), mp.besselk(1, root)
        i0, i1 = mp.besseli(0, root), mp.besseli(1, root)
        if no_flow:
            numerator = k0 * mp.besseli(1, outer) + i0 * mp.besselk(1, outer)
            denominator = k1 * mp.besseli(1, outer) - i1 * mp.besselk(1, outer)
        else:
            numerator = k0 * mp.besseli(0, outer) - i0 * mp.besselk(0, outer)
            denominator = k1 * mp.besseli(0, outer) + i1 * mp.besselk(0, outer)
        return numerator / (s * root * denominator)


@pytest.mark.parametrize("flint_enabled, gmpy2_enabled", BACKEND_COMBINATIONS)
@pytest.mark.parametrize("geometry, no_flow", [
    (transforms.radial_bounded_nf, True), (transforms.radial_bounded_cp, False)])
def test_bounded_radial_boundary_term_matches_direct_evaluation(
    set_backends: Callable[[bool, bool], None], flint_enabled: bool, gmpy2_enabled: bool,
    geometry: Any, no_flow: bool,
) -> None:
    # Pins the "boundary term is negligible -> mp.zero" gate in _bessel.bessel_k_over_i:
    # z = sqrt(s) r_eD spans 0.3 to 77, through the region where the K/I ratio is
    # just at, above and below the bits the result needs.
    set_backends(flint_enabled, gmpy2_enabled)
    mp.dps = 34
    for outer_radius in (3.0, 10.0):
        for point in ("1e-2", "0.05", "0.5", "5.0", "60.0"):
            s = mp.mpf(point)
            value = transforms.lap_transform(s, geometry, reD=outer_radius)
            reference = _radial_bounded_reference(s, outer_radius, no_flow)
            assert abs(value / reference - 1) < mp.mpf(10) ** -31, (outer_radius, point)


def test_variable_storage_has_no_cancellation_at_large_s(monkeypatch: pytest.MonkeyPatch) -> None:
    # W(s) = 1 + cD cDd s / (1 + aD s); the form s^2 cD cDd (1/s - 1/(s + 1/aD))
    # lost ~4e-9 relative at s = 1e8 and 17 digits on the mpf path.
    monkeypatch.setattr(_flint, "HAS_FLINT", False)
    mp.dps = 17
    s = mp.mpf("1e8")
    value = transforms.apply_variable_wbs(mp.mpf(1), s, 100.0, 1.0, 1e7)
    with mp.workdps(60):
        exact = (1 + mp.mpf(100) * s / (1 + mp.mpf(10) ** 7 * s)) / (1 + s ** 2 * 100)
    assert abs(value / exact - 1) < mp.mpf(10) ** -15


@pytest.mark.parametrize("keyword_arguments, message", [
    ({"lam": 1e-6}, "both lam and omega"),
    ({"omega": 0.1}, "both lam and omega"),
    ({"cD": 10.0, "cDd": 1.0}, "both cDd and aD"),
    ({"cD": 10.0, "aD": 1.0}, "both cDd and aD"),
    ({"cD": -5.0}, "cD must be >= 0"),
])
def test_lap_transform_rejects_partial_effects(
    keyword_arguments: Dict[str, float], message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        transforms.lap_transform(mp.mpf(2), transforms.radial_inf, **keyword_arguments)


@requires_scipy
def test_hybrid_rate_integral_average_includes_first_interval() -> None:
    # qDdi(t) = (1/t) int_0^t qD_hybrid; exact from gwr of rate / s plus the closed-form
    # Arps and exponential integrals. The first interval [0, t0] was missing (5% at t = 10).
    times = 10 ** np.linspace(-1, 3, 81)
    rate_laplace = lambda s: transforms.output_rate(  # noqa: E731
        transforms.lap_transform(s, transforms.radial_inf), s)
    rate = np.asarray(gwr(rate_laplace, times, 32))
    _, _, _, integral_average, _ = transforms.hybrid_rate(rate, times, b=0.5)
    cumulative = np.asarray(gwr(lambda s: rate_laplace(s) / s, times, 32))
    b = 0.5
    arps_integral = ((1 + b * times) ** (1 - 1 / b) - 1) / (b * (1 - 1 / b))
    exact = (cumulative + arps_integral - (1 - np.exp(-times))) / times
    later = times >= 10.0
    np.testing.assert_allclose(integral_average[later], exact[later], rtol=5e-3)


@requires_scipy
def test_hybrid_rate_rejects_fewer_than_two_points() -> None:
    with pytest.raises(ValueError, match="length >= 2"):
        transforms.hybrid_rate(np.array([1.0]), np.array([1.0]))


def test_bessel_k_over_i_rejects_other_orders() -> None:
    from gwr_inversion._bessel import bessel_k_over_i
    with pytest.raises(ValueError, match="order 0 or 1"):
        bessel_k_over_i(2, mp.mpf("0.7"), lambda ratio: 0.0)


# ---------------------------------------------------------------------------
# Independent oracles: the formulas written out here from the references, with
# mpmath Bessel functions and quadrature, 30 digits above the working precision.
# The agreement tests above compare paths that share code; these do not.
# ---------------------------------------------------------------------------


def _radial_infinite_reference(s: Any, s_eff: Any) -> Any:
    """Finite-wellbore radial solution K_0(sqrt(s_eff)) / (s sqrt(s_eff) K_1(sqrt(s_eff)))."""
    root = mp.sqrt(s_eff)
    return mp.besselk(0, root) / (s * root * mp.besselk(1, root))


def _with_modifiers_reference(
    s: Any, lam: float, omega: float, skin: float, storage: float,
    storage_change_ratio: Optional[float], change_rate: Optional[float],
) -> Any:
    """Warren-Root s_eff, skin S/s, then van Everdingen-Hurst or Hegeman storage,
    in the published forms (Hegeman: 1 + s^2 C_D C_Dd (1/s - 1/(s + 1/a_D)))."""
    # The library forms omega (1 - omega) and 1 - omega from Python-float parameters
    # (a 1e-16 parameter rounding, not a formula difference); do the same here.
    product_coefficient, complement = float(omega) * (1 - float(omega)), 1 - float(omega)
    s_eff = s * (lam + product_coefficient * s) / (lam + complement * s)
    sandface = _radial_infinite_reference(s, s_eff) + skin / s
    storage_factor = mp.one
    if storage_change_ratio is not None and change_rate is not None:
        storage_factor = 1 + s ** 2 * storage * storage_change_ratio * (
            1 / s - 1 / (s + 1 / mp.mpf(change_rate)))
    return storage_factor / (1 / sandface + s ** 2 * storage)


@pytest.mark.parametrize("flint_enabled, gmpy2_enabled", BACKEND_COMBINATIONS)
@pytest.mark.parametrize("storage_change_ratio, change_rate", [(None, None), (10.0, 0.01)])
def test_modifiers_match_published_formulas(
    set_backends: Callable[[bool, bool], None], flint_enabled: bool, gmpy2_enabled: bool,
    storage_change_ratio: Optional[float], change_rate: Optional[float],
) -> None:
    set_backends(flint_enabled, gmpy2_enabled)
    mp.dps = 34
    for point in LAPLACE_POINTS:
        s = mp.mpf(point)
        value = transforms.lap_transform(
            s, transforms.radial_inf, lam=1e-3, omega=0.2, skin=3.0, cD=100.0,
            cDd=storage_change_ratio, aD=change_rate)
        with mp.workdps(mp.dps + 30):
            reference = _with_modifiers_reference(
                s, mp.mpf(1e-3), mp.mpf(0.2), mp.mpf(3), mp.mpf(100),
                storage_change_ratio, change_rate)
        assert abs(value / reference - 1) < mp.mpf(10) ** -30, point


def _fracture_reference(s: Any, fracture_position: float, outer_radius: Optional[float],
                        no_flow: bool) -> Any:
    """Uniform-flux fracture (Gringarten et al. 1974): (1/2s) int_{-1}^{1} of the
    point-source kernel at |x_D - x'|, by quadrature. The bounded kernel adds
    (K_1/I_1 or -K_0/I_0 at sqrt(s) r_eD) I_0, as in Ozkan and Raghavan (1991)."""
    root = mp.sqrt(s)
    if outer_radius is None:
        boundary = mp.zero
    elif no_flow:
        boundary = mp.besselk(1, root * outer_radius) / mp.besseli(1, root * outer_radius)
    else:
        boundary = -mp.besselk(0, root * outer_radius) / mp.besseli(0, root * outer_radius)

    def kernel(position: Any) -> Any:
        distance = root * abs(fracture_position - position)
        return mp.besselk(0, distance) + boundary * mp.besseli(0, distance)

    return mp.quad(kernel, [-1, fracture_position, 1]) / (2 * s)


@pytest.mark.parametrize("flint_enabled, gmpy2_enabled", BACKEND_COMBINATIONS)
@pytest.mark.parametrize("geometry, outer_radius, no_flow", [
    (transforms.frac_inf, None, True),
    (transforms.frac_bounded_nf, 3.0, True),
    (transforms.frac_bounded_cp, 3.0, False),
])
def test_fracture_geometries_match_quadrature(
    set_backends: Callable[[bool, bool], None], flint_enabled: bool, gmpy2_enabled: bool,
    geometry: Any, outer_radius: Optional[float], no_flow: bool,
) -> None:
    set_backends(flint_enabled, gmpy2_enabled)
    mp.dps = 20
    parameters = {"xD": 0.4} if outer_radius is None else {"xD": 0.4, "reD": outer_radius}
    for point in ("0.01", "0.5", "8.0"):
        s = mp.mpf(point)
        value = transforms.lap_transform(s, geometry, **parameters)
        with mp.workdps(40):
            reference = _fracture_reference(s, 0.4, outer_radius, no_flow)
        # 1e-15: the library forms 1 - x_D and 1 + x_D from the Python float x_D.
        assert abs(value / reference - 1) < mp.mpf(10) ** -15, point


@pytest.mark.skipif(not _flint.HAS_FLINT, reason="python-flint not installed")
def test_arb_compositor_evaluates_every_arb_case() -> None:
    # The agreement test above would compare mpf with mpf if the ARB path declined.
    import gwr_inversion._arb_transforms as arb_transforms
    mp.dps = 34
    for (geometry, parameters), dual_porosity, skin, storage in itertools.product(
            ARB_GEOMETRIES, DUAL_POROSITY, SKIN, STORAGE):
        wellbore = transforms.WellboreEffects(
            skin.get("skin", 0.0), storage.get("cD", 0.0), storage.get("cDd"),
            storage.get("aD"))
        for point in LAPLACE_POINTS:
            value = arb_transforms.lap_transform(
                mp.mpf(point), transforms._ARB_GEOMETRIES[geometry], wellbore,
                dual_porosity.get("lam"), dual_porosity.get("omega"), parameters)
            assert value is not None, (geometry.__name__, point)


def test_log2_of_nan_keeps_every_bit() -> None:
    # A NaN sensitivity must not read as "negligible" (-inf) and drop a boundary term.
    assert transforms._log2(mp.nan) == float("inf")


def _raise_zero_division(ratio: Any) -> float:
    """A sensitivity formula whose denominator is exactly zero."""
    raise ZeroDivisionError


@pytest.mark.parametrize("log2_sensitivity", [
    _raise_zero_division, lambda ratio: float("inf"), lambda ratio: transforms._log2(mp.nan)])
def test_infinite_sensitivity_uses_full_precision_instead_of_raising(
    log2_sensitivity: Callable[[Any], float],
) -> None:
    from gwr_inversion._bessel import bessel_k_over_i

    mp.dps = 20
    for argument in ("0.5", "5.0", "50.0"):
        z = mp.mpf(argument)
        ratio = bessel_k_over_i(1, z, log2_sensitivity)
        assert abs(ratio / (mp.besselk(1, z) / mp.besseli(1, z)) - 1) < mp.mpf(10) ** -18


@requires_flint
@pytest.mark.parametrize("order", [float("nan"), float("inf"), mp.mpc(1, 1), "one"])
def test_arb_declines_non_numeric_order(order: Any) -> None:
    assert _flint.besselk(order, mp.mpf(1)) is None


@pytest.mark.skipif(not _flint.HAS_FLINT, reason="python-flint not installed")
def test_negative_mpf_converts_to_an_exact_ball() -> None:
    mp.dps = 120
    value = -mp.mpf(2) ** -3 / 3
    with _flint.flint_context.workprec(30):
        ball = _flint.mpf_to_arb(value)
    assert ball.rad() == 0
    assert _flint.mpf_from_raw(ball._mpf_) == value


@pytest.mark.skipif(not _gmpy2_bessel.HAS_GMPY2_BESSEL, reason="gmpy2 not installed")
def test_gmpy2_declines_where_exp_underflows() -> None:
    mp.dps = 17
    assert _gmpy2_bessel.besselk_gmpy2(0, mp.mpf("1e9")) is None


@requires_scipy
@pytest.mark.parametrize("rate, time, message", [
    (np.array([1.0, 0.5, 0.25]), np.array([1.0, 0.5, 2.0]), "strictly increasing"),
    (np.array([-1.0, 0.5, 0.25]), np.array([1.0, 2.0, 3.0]), "positive and finite"),
    # 1/t^2: the [0, t0] integral of a power law with exponent <= -1 diverges.
    (np.array([1.0, 0.25, 1 / 9]), np.array([1.0, 2.0, 3.0]), "faster than 1/t"),
])
def test_hybrid_rate_rejects_grids_it_cannot_integrate(
    rate: Any, time: Any, message: str,
) -> None:
    # Subtract the Arps correction so the hybrid rate equals ``rate`` closely enough.
    arps = 1 / (1 + time * 0.5) ** 2 - np.exp(-time)
    with pytest.raises(ValueError, match=message):
        transforms.hybrid_rate(rate - arps, time, b=0.5)


@pytest.mark.parametrize("keyword_arguments, message", [
    ({"cD": float("nan")}, "cD must be >= 0"),
    ({"lam": -1e-3, "omega": 0.5}, "lam > 0 and 0 < omega < 1"),
    ({"lam": 1e-3, "omega": 2.0}, "lam > 0 and 0 < omega < 1"),
])
def test_lap_transform_rejects_out_of_domain_effects(
    keyword_arguments: Dict[str, float], message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        transforms.lap_transform(mp.mpf(2), transforms.radial_inf, **keyword_arguments)


@requires_flint
@pytest.mark.parametrize("order", [mp.mpc(1, 1), "one", float("nan"), float("inf")])
def test_arb_besseli_declines_non_numeric_order(order: Any) -> None:
    assert _flint.besseli(order, mp.mpf("2.5")) is None


@pytest.mark.skipif(not _gmpy2_bessel.HAS_GMPY2_BESSEL, reason="gmpy2 not installed")
def test_gmpy2_besseli_declines_where_exp_overflows() -> None:
    mp.dps = 17
    assert _gmpy2_bessel.besseli_gmpy2(0, mp.mpf("8e8")) is None


def test_variable_storage_without_constant_storage_raises() -> None:
    with pytest.raises(ValueError, match="needs cD > 0"):
        transforms.lap_transform(mp.mpf(2), transforms.radial_inf, cD=0.0, cDd=10.0, aD=0.01)


@requires_scipy
def test_hybrid_rate_rejects_non_positive_arps_exponent() -> None:
    with pytest.raises(ValueError, match="b > 0"):
        transforms.hybrid_rate(np.array([1.0, 0.9, 0.8]), np.array([1.0, 2.0, 3.0]), b=0.0)


def test_mpf_parameters_give_the_float_parameter_result() -> None:
    # mpf parameters must work on every path (the ARB path declines them). Dyadic
    # values, so float and mpf parameters form 1 - x_D and omega (1 - omega) exactly.
    mp.dps = 34
    s = mp.mpf("0.37")
    dyadic_geometries = [(transforms.radial_inf, {}), (transforms.linear_inf, {"xD": 0.25}),
                         (transforms.linear_bounded_nf, {"xD": 0.25}),
                         (transforms.linear_bounded_cp, {"xD": 0.5})]
    for geometry, parameters in dyadic_geometries:
        with_floats = transforms.lap_transform(
            s, geometry, lam=2.0 ** -10, omega=0.25, cD=10.0, cDd=2.0, aD=0.5, **parameters)
        with_mpfs = transforms.lap_transform(
            s, geometry, lam=mp.mpf(2) ** -10, omega=mp.mpf(0.25), cD=mp.mpf(10),
            cDd=mp.mpf(2), aD=mp.mpf(0.5),
            **{name: mp.mpf(value) for name, value in parameters.items()})
        assert abs(with_mpfs / with_floats - 1) < mp.mpf(10) ** -30, geometry.__name__


def test_sensitivity_beyond_float_range_uses_full_precision() -> None:
    # K_0/I_0 at z = 1e400 is near 2^(-2.9e400): float(mp.mag(ratio)) in _log2 raises
    # OverflowError, which the guard must turn into "full precision", not propagate.
    from gwr_inversion._bessel import bessel_k_over_i
    mp.dps = 20
    ratio = bessel_k_over_i(0, mp.mpf("1e400"), transforms._log2)
    assert isinstance(ratio, mp.mpf) and ratio >= 0
