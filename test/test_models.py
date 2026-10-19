"""Tests for gwr_inversion.models (composable model objects)."""
from __future__ import annotations

import importlib.util
import inspect
import pickle
from dataclasses import dataclass
from typing import Any, Dict, Iterator

import numpy as np
import pytest
from mpmath import mp

from gwr_inversion import gwr, shutdown_pool, transforms
from gwr_inversion.models import (
    OUTPUTS, RESERVOIR, DualPorosity, FractureConstantPressure, FractureInfinite,
    FractureNoFlow, Geometry, LaplaceFunction, LinearConstantPressure, LinearInfinite,
    LinearNoFlow, Model, RadialConstantPressure, RadialInfinite, RadialNoFlow, Skin,
    VariableWellboreStorage, WellboreStorage,
)

LAPLACE_POINTS = ["0.001", "0.37", "5.0", "880.0"]
TIMES = [0.1, 1.0, 10.0, 100.0]


@pytest.fixture(autouse=True)
def restore_precision() -> Iterator[None]:
    """Tests here set mp.dps = 34; restore it so later test modules run at their own."""
    with mp.workdps(mp.dps):
        yield

BUILTIN_CASES = [
    (RadialInfinite(), transforms.radial_inf, {}),
    (RadialNoFlow(outer_radius=1000.0), transforms.radial_bounded_nf, {"reD": 1000.0}),
    (RadialConstantPressure(outer_radius=50.0), transforms.radial_bounded_cp, {"reD": 50.0}),
    (LinearInfinite(distance=0.3), transforms.linear_inf, {"xD": 0.3}),
    (LinearNoFlow(distance=0.2), transforms.linear_bounded_nf, {"xD": 0.2}),
    (LinearConstantPressure(distance=0.3), transforms.linear_bounded_cp, {"xD": 0.3}),
    (FractureInfinite(), transforms.frac_inf, {"xD": 0.732}),
    (FractureInfinite(fracture_position=0.4), transforms.frac_inf, {"xD": 0.4}),
    (FractureNoFlow(outer_radius=50.0, fracture_position=0.4), transforms.frac_bounded_nf,
     {"reD": 50.0, "xD": 0.4}),
    (FractureConstantPressure(outer_radius=50.0, fracture_position=0.4),
     transforms.frac_bounded_cp, {"reD": 50.0, "xD": 0.4}),
]


@pytest.mark.parametrize("geometry, geometry_function, parameters", BUILTIN_CASES)
def test_builtin_model_equals_lap_transform(
    geometry: Geometry, geometry_function: Any, parameters: Dict[str, float],
) -> None:
    model = Model(geometry, DualPorosity(1e-6, 0.1), Skin(3.0), WellboreStorage(100.0))
    mp.dps = 34
    for point in LAPLACE_POINTS:
        s = mp.mpf(point)
        expected = transforms.lap_transform(
            s, geometry_function, cD=100.0, lam=1e-6, omega=0.1, skin=3.0, **parameters)
        assert model.laplace_pressure(s) == expected


def test_variable_storage_translates() -> None:
    model = Model(RadialInfinite(), VariableWellboreStorage(1e3, 10.0, 0.01))
    mp.dps = 34
    s = mp.mpf("0.25")
    expected = transforms.lap_transform(
        s, transforms.radial_inf, cD=1e3, cDd=10.0, aD=0.01)
    assert model.laplace_pressure(s) == expected


@dataclass(frozen=True)
class LinearInfiniteByHand(Geometry):
    """User geometry: psd = 1 / (s sqrt(s_eff)), the same as transforms.linear_inf at xD=0."""

    def laplace_pressure(self, s: Any, s_eff_by_region: Dict[str, Any]) -> Any:
        return 1 / (s * mp.sqrt(s_eff_by_region[RESERVOIR]))


@pytest.mark.parametrize("effects", [
    (),
    (DualPorosity(1e-3, 0.2),),
    (Skin(2.0), WellboreStorage(10.0)),
    (DualPorosity(1e-3, 0.2), Skin(2.0), VariableWellboreStorage(10.0, 5.0, 0.1)),
])
def test_generic_path_matches_builtin(effects: Any) -> None:
    mp.dps = 34
    for point in LAPLACE_POINTS:
        s = mp.mpf(point)
        generic = Model(LinearInfiniteByHand(), *effects).laplace_pressure(s)
        builtin = Model(LinearInfinite(), *effects).laplace_pressure(s)
        assert abs(generic / builtin - 1) < mp.mpf(10) ** -30


def test_model_pickles() -> None:
    model = Model(RadialNoFlow(1000.0), DualPorosity(1e-6, 0.1), Skin(3.0),
                  WellboreStorage(100.0))
    assert pickle.loads(pickle.dumps(model)) == model


@pytest.mark.parametrize("effects, error, message", [
    ((Skin(1.0), Skin(2.0)), ValueError, "at most one Skin"),
    ((WellboreStorage(10.0), VariableWellboreStorage(10.0, 1.0, 1.0)), ValueError,
     "at most one wellbore storage"),
    ((DualPorosity(1e-6, 0.1), DualPorosity(1e-5, 0.2)), ValueError,
     "at most one DualPorosity per region"),
    ((DualPorosity(1e-6, 0.1, region="stimulated"),), ValueError, "has no region"),
    (("not an effect",), TypeError, "unknown effect"),
])
def test_validation(effects: Any, error: type, message: str) -> None:
    with pytest.raises(error, match=message):
        Model(RadialInfinite(), *effects)


def test_geometry_type_is_checked() -> None:
    with pytest.raises(TypeError, match="geometry must be a Geometry"):
        Model(transforms.radial_inf)  # type: ignore[arg-type]


def test_linear_infinite_outputs_match_exact_solutions() -> None:
    # Exact inverses for psd = 1 / (s sqrt(s)) (linear flow, x_D = 0):
    # p = 2 sqrt(t/pi), q = 1/sqrt(pi t), t dp/dt = sqrt(t/pi), cumulative = (4/3) t^1.5 / sqrt(pi)
    model = Model(LinearInfinite())
    times = np.array(TIMES)
    np.testing.assert_allclose(model.pressure(times, M=32), 2 * np.sqrt(times / np.pi),
                               rtol=1e-14)
    np.testing.assert_allclose(model.rate(times, M=32), 1 / np.sqrt(np.pi * times), rtol=1e-14)
    np.testing.assert_allclose(model.derivative(times, M=32), np.sqrt(times / np.pi),
                               rtol=1e-14)
    np.testing.assert_allclose(model.cumulative(times, M=32),
                               4 / 3 * times ** 1.5 / np.sqrt(np.pi), rtol=1e-14)


def test_scalar_and_list_inputs() -> None:
    model = Model(RadialInfinite(), WellboreStorage(1e3))
    scalar = model.pressure(100.0, M=16)
    assert isinstance(scalar, float)
    assert model.derivative(100.0, M=16) == pytest.approx(
        100.0 * gwr(model.laplace("derivative"), 100.0, 16), rel=0, abs=0)
    as_list = model.pressure([10.0, 100.0], M=16)
    assert isinstance(as_list, list) and as_list[1] == scalar


def test_laplace_function_is_one_argument_and_pickles() -> None:
    function = Model(RadialNoFlow(1000.0), Skin(3.0)).laplace("rate")
    assert isinstance(function, LaplaceFunction)
    # gwr_single passes precin to a second parameter; LaplaceFunction must take one.
    assert len(inspect.signature(function).parameters) == 1
    assert pickle.loads(pickle.dumps(function)) == function
    mp.dps = 34
    s = mp.mpf("0.5")
    assert function(s) == transforms.output_rate(
        transforms.lap_transform(s, transforms.radial_bounded_nf, reD=1000.0, skin=3.0), s)


def test_unknown_output_raises() -> None:
    with pytest.raises(ValueError, match="output must be one of"):
        Model(RadialInfinite()).laplace("velocity")
    assert OUTPUTS == ("pressure", "rate", "derivative", "cumulative")


def test_derivative_rejects_mpf_dtype() -> None:
    with pytest.raises(ValueError, match="dtype"):
        Model(RadialInfinite()).derivative([1.0], M=8, dtype="mpf")


def test_parallel_matches_sequential_bitwise() -> None:
    model = Model(RadialNoFlow(1000.0), DualPorosity(1e-6, 0.1), WellboreStorage(100.0))
    times = 10 ** np.linspace(1, 6, 7)
    try:
        sequential = model.rate(times, M=8)
        parallel = model.rate(times, M=8, workers=2)
    finally:
        shutdown_pool()
    assert np.array_equal(sequential, parallel)


@pytest.mark.skipif(importlib.util.find_spec("scipy") is None, reason="scipy not installed")
def test_hybrid_rate_matches_transforms() -> None:
    model = Model(RadialNoFlow(24.0))
    times = 10 ** np.linspace(-1, 4, 51)
    rate = np.asarray(model.rate(times, M=8))
    expected = transforms.hybrid_rate(rate, times, b=0.5)
    for value, reference in zip(model.hybrid_rate(times, b=0.5, M=8), expected):
        assert np.array_equal(value, reference)


def test_public_exports() -> None:
    import gwr_inversion
    for name in ["Model", "Geometry", "LaplaceFunction", "OUTPUTS", "RadialInfinite",
                 "RadialNoFlow", "RadialConstantPressure", "LinearInfinite", "LinearNoFlow",
                 "LinearConstantPressure", "FractureInfinite", "FractureNoFlow",
                 "FractureConstantPressure", "DualPorosity", "Skin", "WellboreStorage",
                 "VariableWellboreStorage"]:
        assert hasattr(gwr_inversion, name), name


@pytest.mark.parametrize("time", [
    np.array([[1.0], [10.0], [100.0]]), np.array([[1.0, 10.0, 100.0]]),
    range(1, 4), (value for value in [1.0, 10.0, 100.0]),
])
def test_derivative_matches_gwr_for_every_time_shape(time: Any) -> None:
    model = Model(RadialInfinite(), WellboreStorage(10.0))
    times = [1.0, 10.0, 100.0] if not isinstance(time, range) else [1.0, 2.0, 3.0]
    expected = [time_value * value for time_value, value
                in zip(times, gwr(model.laplace("derivative"), times, 12))]
    derivative = model.derivative(time, M=12)
    assert np.asarray(derivative).shape == (3,)
    np.testing.assert_array_equal(np.asarray(derivative), np.asarray(expected))


@dataclass(frozen=True)
class DoubledRadial(RadialInfinite):
    """A built-in subclass that overrides laplace_pressure: Model must call it."""

    def laplace_pressure(self, s: Any, s_eff_by_region: Dict[str, Any]) -> Any:
        return 2 * super().laplace_pressure(s, s_eff_by_region)


def test_builtin_subclass_override_is_honoured() -> None:
    mp.dps = 34
    s = mp.mpf("0.5")
    plain = Model(RadialInfinite(), Skin(1.0)).laplace_pressure(s)
    doubled = Model(DoubledRadial(), Skin(1.0)).laplace_pressure(s)
    assert abs(doubled - (plain - 1 / s) * 2 - 1 / s) < mp.mpf(10) ** -30


@pytest.mark.parametrize("build, message", [
    (lambda: RadialNoFlow(outer_radius=0.5), "outer_radius"),
    (lambda: RadialConstantPressure(outer_radius=float("nan")), "outer_radius"),
    (lambda: LinearNoFlow(distance=2.0), "distance"),
    (lambda: LinearInfinite(distance=-1.0), "distance"),
    (lambda: FractureNoFlow(outer_radius=100.0, fracture_position=1.5), "fracture_position"),
    (lambda: DualPorosity(1e-6, 1.5), "storativity_ratio"),
    (lambda: DualPorosity(-1e-6, 0.1), "interporosity_flow"),
    (lambda: Skin(float("nan")), "skin"),
    (lambda: WellboreStorage(0.0), "storage"),
    (lambda: WellboreStorage(float("inf")), "storage"),
    (lambda: VariableWellboreStorage(10.0, float("nan"), 1.0), "storage_change_ratio"),
])
def test_parameters_outside_documented_domain_raise(build: Any, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        build()


@pytest.mark.skipif(importlib.util.find_spec("scipy") is None, reason="scipy not installed")
def test_hybrid_rate_takes_order_positionally_like_the_other_outputs() -> None:
    model = Model(RadialNoFlow(24.0))
    times = 10 ** np.linspace(-1, 2, 21)
    by_position = model.hybrid_rate(times, 8)
    by_keyword = model.hybrid_rate(times, M=8, b=0.5)
    for value, reference in zip(by_position, by_keyword):
        assert np.array_equal(value, reference)


def test_laplace_function_rejects_unknown_output_at_construction() -> None:
    with pytest.raises(ValueError, match="output must be one of"):
        LaplaceFunction(Model(RadialInfinite()), "velocity")


@pytest.mark.skipif(importlib.util.find_spec("scipy") is None, reason="scipy not installed")
def test_hybrid_rate_accepts_generator_and_row_vector_time() -> None:
    model = Model(RadialNoFlow(24.0))
    times = 10 ** np.linspace(-1, 2, 11)
    reference = model.hybrid_rate(times, 8)
    for time in ((value for value in times), times.reshape(1, -1)):
        for value, expected in zip(model.hybrid_rate(time, 8), reference):
            assert np.array_equal(value, expected)


def test_cache_clear_clears_the_gaver_tables() -> None:
    import gwr_inversion
    from gwr_inversion import _gwr
    _gwr._gaver_integers(8)
    gwr_inversion.cache_clear()
    assert _gwr._gaver_integers.cache_info().currsize == 0
