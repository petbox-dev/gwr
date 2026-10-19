"""Tests for gwr_inversion.models (composable model objects)."""
from __future__ import annotations

import pickle
from dataclasses import dataclass
from typing import Any, Dict

import pytest
from mpmath import mp

from gwr_inversion import transforms
from gwr_inversion.models import (
    RESERVOIR, DualPorosity, FractureConstantPressure, FractureInfinite, FractureNoFlow,
    Geometry, LinearConstantPressure, LinearInfinite, LinearNoFlow, Model,
    RadialConstantPressure, RadialInfinite, RadialNoFlow, Skin, VariableWellboreStorage,
    WellboreStorage,
)

LAPLACE_POINTS = ["0.001", "0.37", "5.0", "880.0"]

BUILTIN_CASES = [
    (RadialInfinite(), transforms.radial_inf, {}),
    (RadialNoFlow(outer_radius=1000.0), transforms.radial_bounded_nf, {"reD": 1000.0}),
    (RadialConstantPressure(outer_radius=50.0), transforms.radial_bounded_cp, {"reD": 50.0}),
    (LinearInfinite(distance=0.0), transforms.linear_inf, {"xD": 0.0}),
    (LinearNoFlow(distance=0.2), transforms.linear_bounded_nf, {"xD": 0.2}),
    (LinearConstantPressure(distance=0.0), transforms.linear_bounded_cp, {"xD": 0.0}),
    (FractureInfinite(), transforms.frac_inf, {"xD": 0.732}),
    (FractureNoFlow(outer_radius=100.0), transforms.frac_bounded_nf,
     {"reD": 100.0, "xD": 0.732}),
    (FractureConstantPressure(outer_radius=100.0), transforms.frac_bounded_cp,
     {"reD": 100.0, "xD": 0.732}),
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
    ((WellboreStorage(0.0),), ValueError, "storage must be > 0"),
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
