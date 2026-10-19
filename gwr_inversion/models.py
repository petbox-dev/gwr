"""
Composable reservoir models.

A model is a geometry plus effects::

    from gwr_inversion import Model, RadialNoFlow, DualPorosity, Skin, WellboreStorage

    model = Model(RadialNoFlow(outer_radius=1000.0),
                  DualPorosity(interporosity_flow=1e-6, storativity_ratio=0.1),
                  Skin(3.0), WellboreStorage(100.0))

Every effect has a level. Region effects (``DualPorosity``) change the
Laplace variable of one region of the geometry, ``s -> s_eff``; wellbore
effects (``Skin``, ``WellboreStorage``, ``VariableWellboreStorage``) act on
the sandface pressure. The composition order is fixed by level: region
effects, then the geometry, then skin, then storage. Every built-in geometry
has one region, ``"reservoir"``; multi-region geometries (trilinear, EFR)
will name more.

Models are frozen dataclasses of module-level classes, so they pickle and
work with ``gwr(..., workers > 1)``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable, ClassVar, Dict, Iterable, List, Tuple, Union

import numpy as np

from . import transforms
from ._gwr import gwr

RESERVOIR = "reservoir"


def _require_in_range(
    name: str, value: float, lower: float, upper: float, *,
    lower_inclusive: bool = True, upper_inclusive: bool = True,
) -> None:
    """Raise ValueError unless ``value`` is finite and inside the documented range."""
    above_lower = value >= lower if lower_inclusive else value > lower
    below_upper = value <= upper if upper_inclusive else value < upper
    if not (math.isfinite(value) and above_lower and below_upper):
        lower_bracket = "[" if lower_inclusive else "("
        upper_bracket = "]" if upper_inclusive and math.isfinite(upper) else ")"
        raise ValueError(f"{name} must be finite and in {lower_bracket}{lower}, {upper}"
                         f"{upper_bracket}, got {value!r}")


# ---------------------------------------------------------------------------
# Geometries
# ---------------------------------------------------------------------------


class Geometry:
    """
    Base class for a geometry.

    Subclasses set ``regions`` (names of the regions whose Laplace variable a
    region effect may change) and implement ``laplace_pressure``.
    """

    regions: ClassVar[Tuple[str, ...]] = (RESERVOIR,)

    def laplace_pressure(self, s: Any, s_eff_by_region: Dict[str, Any]) -> Any:
        """
        Dimensionless sandface pressure in Laplace space, before wellbore effects.

        Parameters
        ----------
            s : mpf
                Laplace variable.

            s_eff_by_region : dict
                Effective Laplace variable per region name (``s`` itself for a
                region without a region effect).

        Returns
        -------
            mpf
                Dimensionless pressure in Laplace space.
        """
        raise NotImplementedError


class BuiltinGeometry(Geometry):
    """A geometry backed by one of the functions in ``transforms.py``."""

    def geometry_function(self) -> Callable[..., Any]:
        """The ``transforms`` function with signature ``(s, s_eff, **parameters)``."""
        raise NotImplementedError

    def geometry_parameters(self) -> Dict[str, float]:
        """Keyword arguments for ``geometry_function`` (``reD``, ``xD``)."""
        return {}

    def laplace_pressure(self, s: Any, s_eff_by_region: Dict[str, Any]) -> Any:
        """``geometry_function`` at the reservoir's effective Laplace variable.
        ``Model`` calls ``transforms.compose`` instead (same value, ARB fast path)."""
        return self.geometry_function()(
            s, s_eff_by_region[RESERVOIR], **self.geometry_parameters())


@dataclass(frozen=True)
class RadialInfinite(BuiltinGeometry):
    """Infinite-acting radial flow (``transforms.radial_inf``)."""

    def geometry_function(self) -> Callable[..., Any]:
        return transforms.radial_inf


@dataclass(frozen=True)
class RadialNoFlow(BuiltinGeometry):
    """Bounded radial flow, no-flow outer boundary (``transforms.radial_bounded_nf``).

    ``outer_radius`` is r_eD = r_e / r_w (> 1).
    """

    outer_radius: float = 1000.0

    def __post_init__(self) -> None:
        """Reject parameters outside the documented domain."""
        _require_in_range("outer_radius", self.outer_radius, 1.0, math.inf,
                          lower_inclusive=False)

    def geometry_function(self) -> Callable[..., Any]:
        return transforms.radial_bounded_nf

    def geometry_parameters(self) -> Dict[str, float]:
        return {"reD": self.outer_radius}


@dataclass(frozen=True)
class RadialConstantPressure(BuiltinGeometry):
    """Bounded radial flow, constant-pressure outer boundary
    (``transforms.radial_bounded_cp``). ``outer_radius`` is r_eD = r_e / r_w (> 1).
    """

    outer_radius: float = 1000.0

    def __post_init__(self) -> None:
        """Reject parameters outside the documented domain."""
        _require_in_range("outer_radius", self.outer_radius, 1.0, math.inf,
                          lower_inclusive=False)

    def geometry_function(self) -> Callable[..., Any]:
        return transforms.radial_bounded_cp

    def geometry_parameters(self) -> Dict[str, float]:
        return {"reD": self.outer_radius}


@dataclass(frozen=True)
class LinearInfinite(BuiltinGeometry):
    """Infinite-acting linear flow (``transforms.linear_inf``).

    ``distance`` is x_D, the dimensionless distance from the well (0 = at the well).
    """

    distance: float = 0.0

    def __post_init__(self) -> None:
        """Reject parameters outside the documented domain."""
        _require_in_range("distance", self.distance, 0.0, math.inf)

    def geometry_function(self) -> Callable[..., Any]:
        return transforms.linear_inf

    def geometry_parameters(self) -> Dict[str, float]:
        return {"xD": self.distance}


@dataclass(frozen=True)
class LinearNoFlow(BuiltinGeometry):
    """Bounded linear flow, no-flow boundary (``transforms.linear_bounded_nf``).

    ``distance`` is x_D, 0 <= x_D <= 1 (0 = at the well, 1 = at the boundary).
    """

    distance: float = 0.0

    def __post_init__(self) -> None:
        """Reject parameters outside the documented domain."""
        _require_in_range("distance", self.distance, 0.0, 1.0)

    def geometry_function(self) -> Callable[..., Any]:
        return transforms.linear_bounded_nf

    def geometry_parameters(self) -> Dict[str, float]:
        return {"xD": self.distance}


@dataclass(frozen=True)
class LinearConstantPressure(BuiltinGeometry):
    """Bounded linear flow, constant-pressure boundary (``transforms.linear_bounded_cp``).

    ``distance`` is x_D, 0 <= x_D <= 1.
    """

    distance: float = 0.0

    def __post_init__(self) -> None:
        """Reject parameters outside the documented domain."""
        _require_in_range("distance", self.distance, 0.0, 1.0)

    def geometry_function(self) -> Callable[..., Any]:
        return transforms.linear_bounded_cp

    def geometry_parameters(self) -> Dict[str, float]:
        return {"xD": self.distance}


@dataclass(frozen=True)
class FractureInfinite(BuiltinGeometry):
    """Fractured well, infinite-acting (``transforms.frac_inf``).

    ``fracture_position`` is x_D along the fracture half-length (0.732 = the
    infinite-conductivity equivalent point; 0 = uniform flux).
    """

    fracture_position: float = 0.732

    def __post_init__(self) -> None:
        """Reject parameters outside the documented domain."""
        _require_in_range("fracture_position", self.fracture_position, 0.0, 1.0)

    def geometry_function(self) -> Callable[..., Any]:
        return transforms.frac_inf

    def geometry_parameters(self) -> Dict[str, float]:
        return {"xD": self.fracture_position}


@dataclass(frozen=True)
class FractureNoFlow(BuiltinGeometry):
    """Fractured well, no-flow circular outer boundary (``transforms.frac_bounded_nf``).

    ``outer_radius`` is r_eD = r_e / x_f (> 1); ``fracture_position`` as in
    ``FractureInfinite``.
    """

    outer_radius: float = 100.0
    fracture_position: float = 0.732

    def __post_init__(self) -> None:
        """Reject parameters outside the documented domain."""
        _require_in_range("outer_radius", self.outer_radius, 1.0, math.inf,
                          lower_inclusive=False)
        _require_in_range("fracture_position", self.fracture_position, 0.0, 1.0)

    def geometry_function(self) -> Callable[..., Any]:
        return transforms.frac_bounded_nf

    def geometry_parameters(self) -> Dict[str, float]:
        return {"reD": self.outer_radius, "xD": self.fracture_position}


@dataclass(frozen=True)
class FractureConstantPressure(BuiltinGeometry):
    """Fractured well, constant-pressure circular outer boundary
    (``transforms.frac_bounded_cp``). Parameters as in ``FractureNoFlow``.
    """

    outer_radius: float = 100.0
    fracture_position: float = 0.732

    def __post_init__(self) -> None:
        """Reject parameters outside the documented domain."""
        _require_in_range("outer_radius", self.outer_radius, 1.0, math.inf,
                          lower_inclusive=False)
        _require_in_range("fracture_position", self.fracture_position, 0.0, 1.0)

    def geometry_function(self) -> Callable[..., Any]:
        return transforms.frac_bounded_cp

    def geometry_parameters(self) -> Dict[str, float]:
        return {"reD": self.outer_radius, "xD": self.fracture_position}


# ---------------------------------------------------------------------------
# Effects
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DualPorosity:
    """Warren–Root pseudo-steady-state interporosity flow in one region
    (``transforms.apply_dual_porosity``): ``s_eff = f(s) s``.

    ``interporosity_flow`` is lambda; ``storativity_ratio`` is omega (0 < omega < 1).
    """

    interporosity_flow: float
    storativity_ratio: float
    region: str = RESERVOIR
    level: ClassVar[str] = "region"

    def __post_init__(self) -> None:
        """Reject parameters outside the documented domain."""
        _require_in_range("interporosity_flow", self.interporosity_flow, 0.0, math.inf,
                          lower_inclusive=False)
        _require_in_range("storativity_ratio", self.storativity_ratio, 0.0, 1.0,
                          lower_inclusive=False, upper_inclusive=False)

    def effective_s(self, s: Any) -> Any:
        """The region's effective Laplace variable ``f(s) s``."""
        return transforms.apply_dual_porosity(s, self.interporosity_flow,
                                              self.storativity_ratio)


@dataclass(frozen=True)
class Skin:
    """Wellbore skin factor S: ``psd + S / s`` (``transforms.apply_skin``)."""

    skin: float
    level: ClassVar[str] = "wellbore"

    def __post_init__(self) -> None:
        """Reject a non-finite skin."""
        _require_in_range("skin", self.skin, -math.inf, math.inf)


@dataclass(frozen=True)
class WellboreStorage:
    """Constant wellbore storage, ``storage`` = C_D > 0 (``transforms.apply_wbs``)."""

    storage: float
    level: ClassVar[str] = "wellbore"

    def __post_init__(self) -> None:
        """Reject a storage coefficient that is not finite and > 0."""
        _require_in_range("storage", self.storage, 0.0, math.inf, lower_inclusive=False)


@dataclass(frozen=True)
class VariableWellboreStorage:
    """Changing wellbore storage (Hegeman, Hallford and Joseph 1993;
    ``transforms.apply_variable_wbs``).

    ``storage`` is C_D > 0; ``storage_change_ratio`` is C_Dd = (C_Ds - C_DL) / C_DL;
    ``change_rate`` is a_D, the dimensionless rate of the change.
    """

    storage: float
    storage_change_ratio: float
    change_rate: float
    level: ClassVar[str] = "wellbore"

    def __post_init__(self) -> None:
        """Reject a storage coefficient that is not finite and > 0, and non-finite
        change parameters (their signs are not restricted here)."""
        _require_in_range("storage", self.storage, 0.0, math.inf, lower_inclusive=False)
        _require_in_range("storage_change_ratio", self.storage_change_ratio,
                          -math.inf, math.inf)
        _require_in_range("change_rate", self.change_rate, -math.inf, math.inf)


Effect = Union[DualPorosity, Skin, WellboreStorage, VariableWellboreStorage]


def _validate(geometry: Any, effects: Tuple[Any, ...]) -> None:
    """Raise TypeError / ValueError for an invalid geometry-effect combination."""
    if not isinstance(geometry, Geometry):
        raise TypeError(f"geometry must be a Geometry, got {geometry!r}")
    porosity_regions: List[str] = []
    skin_count = 0
    storage_count = 0
    for effect in effects:
        if isinstance(effect, DualPorosity):
            if effect.region not in geometry.regions:
                raise ValueError(
                    f"{type(geometry).__name__} has no region {effect.region!r}; "
                    f"its regions are {geometry.regions}")
            if effect.region in porosity_regions:
                raise ValueError(
                    f"at most one DualPorosity per region; region {effect.region!r} has two")
            porosity_regions.append(effect.region)
        elif isinstance(effect, Skin):
            skin_count += 1
        elif isinstance(effect, (WellboreStorage, VariableWellboreStorage)):
            storage_count += 1
        else:
            raise TypeError(f"unknown effect {effect!r}")
    if skin_count > 1:
        raise ValueError("at most one Skin")
    if storage_count > 1:
        raise ValueError(
            "at most one wellbore storage effect (WellboreStorage or VariableWellboreStorage)")


def _wellbore_effects(effects: Tuple[Effect, ...]) -> transforms.WellboreEffects:
    """The skin and storage effects of a validated model as one bundle."""
    skin = next((effect.skin for effect in effects if isinstance(effect, Skin)), 0.0)
    storage = next((effect for effect in effects
                    if isinstance(effect, (WellboreStorage, VariableWellboreStorage))), None)
    if isinstance(storage, VariableWellboreStorage):
        return transforms.WellboreEffects(
            skin, storage.storage, storage.storage_change_ratio, storage.change_rate)
    return transforms.WellboreEffects(skin, 0.0 if storage is None else storage.storage)


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------


@dataclass(frozen=True, init=False)
class Model:
    """A geometry composed with region and wellbore effects."""

    geometry: Geometry
    effects: Tuple[Effect, ...]
    # Derived once from ``effects`` (the model is frozen); laplace_pressure runs
    # millions of times per inversion.
    dual_porosity_by_region: Dict[str, DualPorosity] = field(
        init=False, repr=False, compare=False)
    wellbore: transforms.WellboreEffects = field(init=False, repr=False, compare=False)
    # True for an unmodified built-in geometry: one transforms.compose call gives
    # the same value as the generic path, plus the ARB fast path. A subclass that
    # adds regions or overrides laplace_pressure takes the generic path.
    uses_compose: bool = field(init=False, repr=False, compare=False)

    def __init__(self, geometry: Geometry, *effects: Effect) -> None:
        _validate(geometry, effects)
        object.__setattr__(self, "geometry", geometry)
        object.__setattr__(self, "effects", tuple(effects))
        object.__setattr__(self, "dual_porosity_by_region", {
            effect.region: effect for effect in effects if isinstance(effect, DualPorosity)})
        object.__setattr__(self, "wellbore", _wellbore_effects(effects))
        object.__setattr__(self, "uses_compose", (
            isinstance(geometry, BuiltinGeometry)
            and type(geometry).laplace_pressure is BuiltinGeometry.laplace_pressure
            and geometry.regions == (RESERVOIR,)))

    def laplace_pressure(self, s: Any) -> Any:
        """
        Dimensionless wellbore pressure ``pwd(s)`` in Laplace space.

        Unmodified built-in geometries go through one ``transforms.compose``
        call (same result, and the ARB fast path for radial_inf and the linear
        geometries). Other geometries use the generic path: region effects,
        ``geometry.laplace_pressure``, skin, then storage.
        """
        if self.uses_compose:
            assert isinstance(self.geometry, BuiltinGeometry)
            dual_porosity = self.dual_porosity_by_region.get(RESERVOIR)
            return transforms.compose(
                s,
                self.geometry.geometry_function(),
                self.wellbore,
                None if dual_porosity is None else dual_porosity.interporosity_flow,
                None if dual_porosity is None else dual_porosity.storativity_ratio,
                self.geometry.geometry_parameters(),
            )

        s_eff_by_region = {
            region: (self.dual_porosity_by_region[region].effective_s(s)
                     if region in self.dual_porosity_by_region else s)
            for region in self.geometry.regions
        }
        pressure = self.geometry.laplace_pressure(s, s_eff_by_region)
        return self.wellbore.apply(pressure, s)

    def laplace(self, output: str = "pressure") -> LaplaceFunction:
        """The Laplace-domain function for ``output`` (one of ``OUTPUTS``),
        ready for ``gwr()`` or ``talbot()``; it pickles, so ``workers > 1`` works."""
        if output not in OUTPUTS:
            raise ValueError(f"output must be one of {OUTPUTS}, got {output!r}")
        return LaplaceFunction(self, output)

    def pressure(self, time: Any, M: int = 32, **gwr_options: Any) -> Any:
        """Dimensionless wellbore pressure p_wD(t)."""
        return gwr(self.laplace("pressure"), time, M, **gwr_options)

    def rate(self, time: Any, M: int = 32, **gwr_options: Any) -> Any:
        """Dimensionless rate q_D(t) at constant wellbore pressure."""
        return gwr(self.laplace("rate"), time, M, **gwr_options)

    def cumulative(self, time: Any, M: int = 32, **gwr_options: Any) -> Any:
        """Time integral of p_wD."""
        return gwr(self.laplace("cumulative"), time, M, **gwr_options)

    def derivative(self, time: Any, M: int = 32, **gwr_options: Any) -> Any:
        """Bourdet derivative ``t dp_wD/dt`` (float output only)."""
        if gwr_options.get("dtype", "float") != "float":
            raise ValueError(
                "derivative() supports dtype='float' only; for mpf output use "
                "gwr(model.laplace('derivative'), ...) and multiply by t")
        if isinstance(time, Iterable) and not isinstance(time, np.ndarray):
            time = list(time)  # read twice: by gwr and by the multiplication below
        values = gwr(self.laplace("derivative"), time, M, **gwr_options)
        if isinstance(values, np.ndarray):
            # gwr drops length-1 axes and returns 0-d input as shape (1,); so does ravel
            return np.ravel(np.asarray(time, dtype=np.float64)) * values
        if isinstance(values, list):
            return [float(time_value) * value for time_value, value in zip(time, values)]
        return float(time) * values

    def hybrid_rate(
        self, time: Any, M: int = 32, *, b: float = 0.5, **gwr_options: Any,
    ) -> Tuple["np.ndarray[Any, Any]", "np.ndarray[Any, Any]", "np.ndarray[Any, Any]",
               "np.ndarray[Any, Any]", "np.ndarray[Any, Any]"]:
        """``transforms.hybrid_rate`` applied to this model's rate on ``time``."""
        if isinstance(time, Iterable) and not isinstance(time, np.ndarray):
            time = list(time)
        times = np.ravel(np.asarray(time, dtype=np.float64))  # as gwr drops length-1 axes
        rate = np.asarray(self.rate(times, M, **gwr_options), dtype=np.float64)
        return transforms.hybrid_rate(rate, times, b=b)


_OUTPUT_TRANSFORMS: Dict[str, Callable[[Any, Any], Any]] = {
    "pressure": transforms.output_pressure,
    "rate": transforms.output_rate,
    "derivative": transforms.output_derivative,
    "cumulative": transforms.output_cumulative,
}

OUTPUTS: Tuple[str, ...] = tuple(_OUTPUT_TRANSFORMS)


@dataclass(frozen=True)
class LaplaceFunction:
    """``s -> output(model.laplace_pressure(s), s)``: one parameter for ``gwr()``,
    and picklable for ``workers > 1``."""

    model: Model
    output: str

    def __post_init__(self) -> None:
        """Reject an unknown output here, not at the first evaluation in a worker."""
        if self.output not in _OUTPUT_TRANSFORMS:
            raise ValueError(f"output must be one of {OUTPUTS}, got {self.output!r}")

    def __call__(self, s: Any) -> Any:
        """The output's Laplace-domain value at ``s``."""
        return _OUTPUT_TRANSFORMS[self.output](self.model.laplace_pressure(s), s)
