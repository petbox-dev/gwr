"""
Composable Laplace-Domain Transforms for Petroleum Reservoir Engineering.

Provides geometry functions, physics modifiers, and a compositor for building
Laplace-domain wellbore pressure/rate solutions from interchangeable parts.

References
----------
    Ozkan, E. and Raghavan, R. 1991. New Solutions for Well-Test-Analysis
    Problems: Part 1 — Analytical Considerations. SPE Formation Evaluation,
    6 (3): 359–368. SPE-18615-PA. https://doi.org/10.2118/18615-PA.

    Ozkan, E. and Raghavan, R. 1991. New Solutions for Well-Test-Analysis
    Problems: Part 2 — Computational Considerations and Applications.
    SPE Formation Evaluation, 6 (3): 369–378. SPE-18616-PA.
    https://doi.org/10.2118/18616-PA.
"""
from __future__ import annotations

from typing import Any, Callable, Dict, Optional, Tuple

import numpy as np
from mpmath import mp

from . import _arb_transforms, _flint
from ._bessel import (
    bessel_k_over_i as _bessel_k_over_i,
    besseli as _besseli,
    besseli0_integral as _besseli0_integral,
    besselk as _besselk,
    besselk0_integral as _besselk0_integral,
)
from ._modifiers import (  # noqa: F401  (physics modifiers, re-exported)
    WellboreEffects,
    apply_dual_porosity,
    apply_skin,
    apply_variable_wbs,
    apply_wbs,
    effective_s,
)


def _log2(value: Any) -> float:
    """Upper bound on log2 of a non-negative mpf (mpmath's ``mag``); -inf for zero.

    NaN maps to +inf ("cannot bound it, keep every bit"), never to -inf, which
    would mark the boundary term negligible.
    """
    if mp.isnan(value):
        return float("inf")
    return float(mp.mag(value)) if value > 0 else float("-inf")


# ============================================================================
# Geometry Functions
# ============================================================================
#
# Each geometry function has the signature:
#     geometry_fn(s, s_eff, **kwargs) -> psd
#
# where ``s`` is the Laplace variable, ``s_eff`` is the effective Laplace
# variable (modified by interporosity flow when applicable), and ``psd`` is
# the dimensionless pressure in Laplace space.


def radial_inf(s: Any, s_eff: Any) -> Any:
    """
    Infinite-acting radial flow.

    The fundamental radial flow solution for a line-source well in an
    infinite homogeneous reservoir.

    Ozkan, E. and Raghavan, R. 1991. New Solutions for Well-Test-Analysis
    Problems: Part 1 — Analytical Considerations. SPE Formation Evaluation,
    6 (3): 359–368. SPE-18615-PA. https://doi.org/10.2118/18615-PA.

    .. math::

        \\bar{p}_D = \\frac{K_0(\\sqrt{s_\\mathrm{eff}})}{s \\, \\sqrt{s_\\mathrm{eff}} \\, K_1(\\sqrt{s_\\mathrm{eff}})}

    Parameters
    ----------
        s : mpf
            Laplace variable.

        s_eff : mpf
            Effective Laplace variable. Equal to ``s`` for homogeneous
            reservoirs, or ``f(s) * s`` for naturally fractured reservoirs.

    Returns
    -------
        mpf
            Dimensionless pressure in Laplace space.
    """
    sqrt_s = mp.sqrt(s_eff)
    k0 = _besselk(0, sqrt_s)
    k1 = _besselk(1, sqrt_s)
    return k0 / (sqrt_s * k1) / s


def radial_bounded_nf(s: Any, s_eff: Any, reD: float = 1000.0) -> Any:
    """
    Bounded radial flow with no-flow outer boundary.

    Radial flow solution for a line-source well in a finite reservoir with
    a sealed (no-flow) external boundary.

    Ozkan, E. and Raghavan, R. 1991. New Solutions for Well-Test-Analysis
    Problems: Part 1 — Analytical Considerations. SPE Formation Evaluation,
    6 (3): 359–368. SPE-18615-PA. https://doi.org/10.2118/18615-PA.

    .. math::

        \\bar{p}_D = \\frac{K_0(\\sqrt{s_\\mathrm{eff}}) \\, I_1(\\sqrt{s_\\mathrm{eff}} \\, r_{eD})
        + I_0(\\sqrt{s_\\mathrm{eff}}) \\, K_1(\\sqrt{s_\\mathrm{eff}} \\, r_{eD})}
        {s \\, \\sqrt{s_\\mathrm{eff}} \\left[ K_1(\\sqrt{s_\\mathrm{eff}}) \\, I_1(\\sqrt{s_\\mathrm{eff}} \\, r_{eD})
        - I_1(\\sqrt{s_\\mathrm{eff}}) \\, K_1(\\sqrt{s_\\mathrm{eff}} \\, r_{eD}) \\right]}

    Parameters
    ----------
        s : mpf
            Laplace variable.

        s_eff : mpf
            Effective Laplace variable.

        reD : float
            Dimensionless external radius, ``r_e / r_w``. Must be > 1.

    Returns
    -------
        mpf
            Dimensionless pressure in Laplace space.
    """
    sqrt_s = mp.sqrt(s_eff)

    k0 = _besselk(0, sqrt_s)
    i0 = _besseli(0, sqrt_s)
    k1 = _besselk(1, sqrt_s)
    i1 = _besseli(1, sqrt_s)

    # Divided through by I_1(sqrt_s reD): the boundary enters only through
    # ratio = K_1/I_1 at sqrt_s reD, evaluated to the bits it contributes.
    def log2_sensitivity(ratio: Any) -> float:
        return _log2(i0 * ratio / (k0 + i0 * ratio) + i1 * ratio / abs(k1 - i1 * ratio))

    ratio = _bessel_k_over_i(1, sqrt_s * reD, log2_sensitivity)
    return (k0 + i0 * ratio) / (sqrt_s * (k1 - i1 * ratio)) / s


def radial_bounded_cp(s: Any, s_eff: Any, reD: float = 1000.0) -> Any:
    """
    Bounded radial flow with constant-pressure outer boundary.

    Radial flow solution for a line-source well in a finite reservoir with
    a constant-pressure external boundary.

    Ozkan, E. and Raghavan, R. 1991. New Solutions for Well-Test-Analysis
    Problems: Part 1 — Analytical Considerations. SPE Formation Evaluation,
    6 (3): 359–368. SPE-18615-PA. https://doi.org/10.2118/18615-PA.

    .. math::

        \\bar{p}_D = \\frac{K_0(\\sqrt{s_\\mathrm{eff}}) \\, I_0(\\sqrt{s_\\mathrm{eff}} \\, r_{eD})
        - I_0(\\sqrt{s_\\mathrm{eff}}) \\, K_0(\\sqrt{s_\\mathrm{eff}} \\, r_{eD})}
        {s \\, \\sqrt{s_\\mathrm{eff}} \\left[ K_1(\\sqrt{s_\\mathrm{eff}}) \\, I_0(\\sqrt{s_\\mathrm{eff}} \\, r_{eD})
        + I_1(\\sqrt{s_\\mathrm{eff}}) \\, K_0(\\sqrt{s_\\mathrm{eff}} \\, r_{eD}) \\right]}

    Parameters
    ----------
        s : mpf
            Laplace variable.

        s_eff : mpf
            Effective Laplace variable.

        reD : float
            Dimensionless external radius, ``r_e / r_w``. Must be > 1.

    Returns
    -------
        mpf
            Dimensionless pressure in Laplace space.
    """
    sqrt_s = mp.sqrt(s_eff)

    k0 = _besselk(0, sqrt_s)
    i0 = _besseli(0, sqrt_s)
    k1 = _besselk(1, sqrt_s)
    i1 = _besseli(1, sqrt_s)

    # Divided through by I_0(sqrt_s reD): the boundary enters only through
    # ratio = K_0/I_0 at sqrt_s reD, evaluated to the bits it contributes.
    def log2_sensitivity(ratio: Any) -> float:
        return _log2(i0 * ratio / abs(k0 - i0 * ratio) + i1 * ratio / (k1 + i1 * ratio))

    ratio = _bessel_k_over_i(0, sqrt_s * reD, log2_sensitivity)
    return (k0 - i0 * ratio) / (sqrt_s * (k1 + i1 * ratio)) / s


# ---------------------------------------------------------------------------
# Linear Flow
# ---------------------------------------------------------------------------


def linear_inf(s: Any, s_eff: Any, xD: float = 0.0) -> Any:
    """
    Infinite-acting linear flow.

    Linear flow solution for an infinite-acting slab reservoir (e.g.,
    infinite-conductivity fracture in an infinite reservoir).

    Ozkan, E. and Raghavan, R. 1991. New Solutions for Well-Test-Analysis
    Problems: Part 1 — Analytical Considerations. SPE Formation Evaluation,
    6 (3): 359–368. SPE-18615-PA. https://doi.org/10.2118/18615-PA.

    .. math::

        \\bar{p}_D = \\frac{\\exp(-\\sqrt{s_\\mathrm{eff}} \\, x_D)}{s \\, \\sqrt{s_\\mathrm{eff}}}

    Parameters
    ----------
        s : mpf
            Laplace variable.

        s_eff : mpf
            Effective Laplace variable.

        xD : float
            Dimensionless distance from the wellbore, ``0 <= xD``.
            Default ``0.0`` (at the well).

    Returns
    -------
        mpf
            Dimensionless pressure in Laplace space.
    """
    sqrt_s = mp.sqrt(s_eff)
    return mp.exp(-sqrt_s * xD) / (s * sqrt_s)


def linear_bounded_nf(s: Any, s_eff: Any, xD: float = 0.0) -> Any:
    """
    Linear flow with no-flow outer boundary.

    Linear flow solution for a finite slab reservoir with a sealed
    (no-flow) external boundary.

    Ozkan, E. and Raghavan, R. 1991. New Solutions for Well-Test-Analysis
    Problems: Part 1 — Analytical Considerations. SPE Formation Evaluation,
    6 (3): 359–368. SPE-18615-PA. https://doi.org/10.2118/18615-PA.

    .. math::

        \\bar{p}_D = \\frac{\\cosh(\\sqrt{s_\\mathrm{eff}} \\, (1 - x_D))}
        {s \\, \\sqrt{s_\\mathrm{eff}} \\, \\sinh(\\sqrt{s_\\mathrm{eff}})}

    Parameters
    ----------
        s : mpf
            Laplace variable.

        s_eff : mpf
            Effective Laplace variable.

        xD : float
            Dimensionless distance from the wellbore, ``0 <= xD <= 1``.
            Default ``0.0`` (at the well).

    Returns
    -------
        mpf
            Dimensionless pressure in Laplace space.
    """
    sqrt_s = mp.sqrt(s_eff)
    return mp.cosh(sqrt_s * (1 - xD)) / (s * sqrt_s * mp.sinh(sqrt_s))


def linear_bounded_cp(s: Any, s_eff: Any, xD: float = 0.0) -> Any:
    """
    Linear flow with constant-pressure outer boundary.

    Linear flow solution for a finite slab reservoir with a constant-pressure
    external boundary.

    Ozkan, E. and Raghavan, R. 1991. New Solutions for Well-Test-Analysis
    Problems: Part 1 — Analytical Considerations. SPE Formation Evaluation,
    6 (3): 359–368. SPE-18615-PA. https://doi.org/10.2118/18615-PA.

    .. math::

        \\bar{p}_D = \\frac{\\sinh(\\sqrt{s_\\mathrm{eff}} \\, (1 - x_D))}
        {s \\, \\sqrt{s_\\mathrm{eff}} \\, \\cosh(\\sqrt{s_\\mathrm{eff}})}

    Parameters
    ----------
        s : mpf
            Laplace variable.

        s_eff : mpf
            Effective Laplace variable.

        xD : float
            Dimensionless distance from the wellbore, ``0 <= xD <= 1``.
            Default ``0.0`` (at the well).

    Returns
    -------
        mpf
            Dimensionless pressure in Laplace space.
    """
    sqrt_s = mp.sqrt(s_eff)
    return mp.sinh(sqrt_s * (1 - xD)) / (s * sqrt_s * mp.cosh(sqrt_s))


# ---------------------------------------------------------------------------
# Fractured Well (line-source integration)
# ---------------------------------------------------------------------------


def frac_inf(s: Any, s_eff: Any, xD: float = 0.732) -> Any:
    """
    Infinite-acting fractured well (uniform-flux vertical fracture).

    Line-source integration of K_0 over the fracture half-length for a
    vertically fractured well in an infinite reservoir.

    Gringarten, A. C., Ramey, H. J., Jr., and Raghavan, R. 1974.
    Unsteady-State Pressure Distributions Created by a Well With a Single
    Infinite-Conductivity Vertical Fracture. SPE Journal, 14 (4): 347–360.
    SPE-4051-PA. https://doi.org/10.2118/4051-PA.

    .. math::

        \\bar{p}_D = \\frac{1}{2 s \\sqrt{s_\\mathrm{eff}}} \\left[
        \\int_0^{\\sqrt{s_\\mathrm{eff}}(1 - x_D)} K_0(u) \\, du
        + \\int_0^{\\sqrt{s_\\mathrm{eff}}(1 + x_D)} K_0(u) \\, du \\right]

    Parameters
    ----------
        s : mpf
            Laplace variable.

        s_eff : mpf
            Effective Laplace variable.

        xD : float
            Dimensionless position along the fracture half-length.
            Default ``0.732`` (uniform-flux equivalent point).

    Returns
    -------
        mpf
            Dimensionless pressure in Laplace space.
    """
    sqrt_s = mp.sqrt(s_eff)
    line_source_prefactor = 1 / (2 * s * sqrt_s)

    k0_int_1 = _besselk0_integral(sqrt_s * (1 - xD))
    k0_int_2 = _besselk0_integral(sqrt_s * (1 + xD))

    return line_source_prefactor * (k0_int_1 + k0_int_2)


def frac_bounded_nf(s: Any, s_eff: Any, reD: float = 100.0,
                    xD: float = 0.732) -> Any:
    """
    Fractured well with no-flow outer boundary.

    Vertically fractured well in a bounded cylindrical reservoir with a
    sealed external boundary. Combines the infinite-acting K_0 integral
    with a bounded correction using I_0 integrals.

    Gringarten, A. C., Ramey, H. J., Jr., and Raghavan, R. 1974.
    Unsteady-State Pressure Distributions Created by a Well With a Single
    Infinite-Conductivity Vertical Fracture. SPE Journal, 14 (4): 347–360.
    SPE-4051-PA. https://doi.org/10.2118/4051-PA.

    .. math::

        \\bar{p}_D = \\frac{1}{2 s \\sqrt{s_\\mathrm{eff}}} \\left[
        \\left( \\int_0^{\\sqrt{s_\\mathrm{eff}}(1-x_D)} K_0 \\, du
        + \\int_0^{\\sqrt{s_\\mathrm{eff}}(1+x_D)} K_0 \\, du \\right)
        + \\frac{K_1(\\sqrt{s_\\mathrm{eff}} \\, r_{eD})}{I_1(\\sqrt{s_\\mathrm{eff}} \\, r_{eD})}
        \\left( \\int_0^{\\sqrt{s_\\mathrm{eff}}(1-x_D)} I_0 \\, du
        + \\int_0^{\\sqrt{s_\\mathrm{eff}}(1+x_D)} I_0 \\, du \\right)
        \\right]

    The K_0 and I_0 integrals use closed forms (modified Struve and
    1F2 hypergeometric functions), not numerical quadrature.

    Parameters
    ----------
        s : mpf
            Laplace variable.

        s_eff : mpf
            Effective Laplace variable.

        reD : float
            Dimensionless external radius, ``r_e / x_f``. Must be > 1.

        xD : float
            Dimensionless position along the fracture half-length.
            Default ``0.732`` (uniform-flux equivalent point).

    Returns
    -------
        mpf
            Dimensionless pressure in Laplace space.
    """
    sqrt_s = mp.sqrt(s_eff)
    line_source_prefactor = 1 / (2 * s * sqrt_s)

    # Infinite-acting component
    k0_int_1 = _besselk0_integral(sqrt_s * (1 - xD))
    k0_int_2 = _besselk0_integral(sqrt_s * (1 + xD))
    inf_part = k0_int_1 + k0_int_2

    # Bounded correction: ratio = K_1/I_1 at sqrt_s reD, to the bits it contributes
    i0_int_sum = (_besseli0_integral(sqrt_s * (1 - xD))
                  + _besseli0_integral(sqrt_s * (1 + xD)))

    def log2_sensitivity(ratio: Any) -> float:
        return _log2(ratio * i0_int_sum / abs(inf_part + ratio * i0_int_sum))

    ratio = _bessel_k_over_i(1, sqrt_s * reD, log2_sensitivity)
    return line_source_prefactor * (inf_part + ratio * i0_int_sum)


def frac_bounded_cp(s: Any, s_eff: Any, reD: float = 100.0,
                    xD: float = 0.732) -> Any:
    """
    Fractured well with constant-pressure outer boundary.

    Vertically fractured well in a bounded cylindrical reservoir with a
    constant-pressure external boundary.

    Gringarten, A. C., Ramey, H. J., Jr., and Raghavan, R. 1974.
    Unsteady-State Pressure Distributions Created by a Well With a Single
    Infinite-Conductivity Vertical Fracture. SPE Journal, 14 (4): 347–360.
    SPE-4051-PA. https://doi.org/10.2118/4051-PA.

    .. math::

        \\bar{p}_D = \\frac{1}{2 s \\sqrt{s_\\mathrm{eff}}} \\left[
        \\left( \\int_0^{\\sqrt{s_\\mathrm{eff}}(1-x_D)} K_0 \\, du
        + \\int_0^{\\sqrt{s_\\mathrm{eff}}(1+x_D)} K_0 \\, du \\right)
        - \\frac{K_0(\\sqrt{s_\\mathrm{eff}} \\, r_{eD})}{I_0(\\sqrt{s_\\mathrm{eff}} \\, r_{eD})}
        \\left( \\int_0^{\\sqrt{s_\\mathrm{eff}}(1-x_D)} I_0 \\, du
        + \\int_0^{\\sqrt{s_\\mathrm{eff}}(1+x_D)} I_0 \\, du \\right)
        \\right]

    Parameters
    ----------
        s : mpf
            Laplace variable.

        s_eff : mpf
            Effective Laplace variable.

        reD : float
            Dimensionless external radius, ``r_e / x_f``. Must be > 1.

        xD : float
            Dimensionless position along the fracture half-length.
            Default ``0.732`` (uniform-flux equivalent point).

    Returns
    -------
        mpf
            Dimensionless pressure in Laplace space.
    """
    sqrt_s = mp.sqrt(s_eff)
    line_source_prefactor = 1 / (2 * s * sqrt_s)

    # Infinite-acting component
    k0_int_1 = _besselk0_integral(sqrt_s * (1 - xD))
    k0_int_2 = _besselk0_integral(sqrt_s * (1 + xD))
    inf_part = k0_int_1 + k0_int_2

    # Constant-pressure correction: ratio = K_0/I_0 at sqrt_s reD, to the bits it contributes
    i0_int_sum = (_besseli0_integral(sqrt_s * (1 - xD))
                  + _besseli0_integral(sqrt_s * (1 + xD)))

    def log2_sensitivity(ratio: Any) -> float:
        return _log2(ratio * i0_int_sum / abs(inf_part - ratio * i0_int_sum))

    ratio = _bessel_k_over_i(0, sqrt_s * reD, log2_sensitivity)
    return line_source_prefactor * (inf_part - ratio * i0_int_sum)


# ============================================================================
# Output Transforms
# ============================================================================


def output_pressure(pwd: Any, s: Any) -> Any:
    """
    Identity output — dimensionless pressure in Laplace space.

    .. math::

        \\bar{p}_{wD}

    Parameters
    ----------
        pwd : mpf
            Dimensionless wellbore pressure in Laplace space.

        s : mpf
            Laplace variable.

    Returns
    -------
        mpf
            Dimensionless pressure in Laplace space.
    """
    return pwd


def output_rate(pwd: Any, s: Any) -> Any:
    """
    Dimensionless rate in Laplace space via the constant-rate/constant-pressure
    duality.

    van Everdingen, A. F. and Hurst, W. 1949. The Application of the
    Laplace Transformation to Flow Problems in Reservoirs. Journal of
    Petroleum Technology, 1 (12): 305–324. SPE-949305-G.
    https://doi.org/10.2118/949305-G.

    .. math::

        \\bar{q}_D = \\frac{1}{s^2 \\, \\bar{p}_{wD}}

    Parameters
    ----------
        pwd : mpf
            Dimensionless wellbore pressure in Laplace space.

        s : mpf
            Laplace variable.

    Returns
    -------
        mpf
            Dimensionless rate in Laplace space.
    """
    return 1 / (s ** 2 * pwd)


def output_cumulative(pwd: Any, s: Any) -> Any:
    """
    Time integral of the dimensionless wellbore pressure, in Laplace space
    (for constant-pressure cumulative production, apply the same division by
    ``s`` to the rate).

    .. math::

        \\mathcal{L}\\left[\\int_0^t p_{wD} \\, d\\tau\\right] = \\frac{\\bar{p}_{wD}}{s}

    Parameters
    ----------
        pwd : mpf
            Dimensionless wellbore pressure in Laplace space.

        s : mpf
            Laplace variable.

    Returns
    -------
        mpf
            Laplace transform of the time integral of ``p_wD``.
    """
    return pwd / s


def output_derivative(pwd: Any, s: Any) -> Any:
    """
    Laplace transform of ``dp_wD/dt``, from which the Bourdet derivative follows.

    Multiplication by ``s`` in the Laplace domain is differentiation in
    time (``p_wD(0) = 0``), so the inverse is :math:`dp_{wD}/dt`; multiplying
    it by ``t`` after inversion gives the Bourdet derivative
    :math:`t \\, dp_{wD}/dt`.

    Bourdet, D., Whittle, T. M., Douglas, A. A., and Pirard, Y. M. 1983.
    A New Set of Type Curves Simplifies Well Test Analysis. World Oil,
    196 (6): 95–106.

    .. math::

        s \\, \\bar{p}_{wD}

    Parameters
    ----------
        pwd : mpf
            Dimensionless wellbore pressure in Laplace space.

        s : mpf
            Laplace variable.

    Returns
    -------
        mpf
            ``s p_wD``, the transform of ``dp_wD/dt``. Multiply by ``t``
            after inversion to obtain the time-domain derivative.
    """
    return s * pwd


# ============================================================================
# Compositor
# ============================================================================


# Ball-arithmetic versions of the geometries where they measured faster; used
# only while ``_flint.HAS_FLINT``. See _arb_transforms for why the others are not here.
_ARB_GEOMETRIES: Dict[Callable[..., Any], Callable[..., Any]] = {
    radial_inf: _arb_transforms.radial_inf,
    linear_inf: _arb_transforms.linear_inf,
    linear_bounded_nf: _arb_transforms.linear_bounded_nf,
    linear_bounded_cp: _arb_transforms.linear_bounded_cp,
}


def lap_transform(
    s: Any,
    geometry_fn: Callable[..., Any],
    cD: float = 0.0,
    cDd: Optional[float] = None,
    aD: Optional[float] = None,
    lam: Optional[float] = None,
    omega: Optional[float] = None,
    skin: float = 0.0,
    **geom_kwargs: Any,
) -> Any:
    """
    Compose a Laplace-domain wellbore pressure solution from building blocks.

    Applies physics modifiers in standard order:

    1. Interporosity flow (naturally fractured): ``s`` → ``s_eff``
    2. Geometry function: ``psd(s, s_eff, **geom_kwargs)``
    3. Skin: ``psd + skin / s``
    4. Wellbore storage (constant or variable): ``pwd``

    Ozkan, E. and Raghavan, R. 1991. New Solutions for Well-Test-Analysis
    Problems: Part 2 — Computational Considerations and Applications.
    SPE Formation Evaluation, 6 (3): 369–378. SPE-18616-PA.
    https://doi.org/10.2118/18616-PA.

    Parameters
    ----------
        s : mpf
            Laplace variable.

        geometry_fn : callable
            A geometry function with signature ``(s, s_eff, **kwargs) -> psd``.
            One of :func:`radial_inf`, :func:`radial_bounded_nf`,
            :func:`radial_bounded_cp`, :func:`linear_inf`,
            :func:`linear_bounded_nf`, :func:`linear_bounded_cp`,
            :func:`frac_inf`, :func:`frac_bounded_nf`, :func:`frac_bounded_cp`.

        cD : float
            Dimensionless wellbore storage coefficient ``C_D``.
            Default ``0.0`` (no wellbore storage).

        cDd : float or None
            Dimensionless wellbore storage change ratio for variable WBS.
            Requires ``aD`` to also be set.

        aD : float or None
            Dimensionless wellbore storage change rate for variable WBS.
            Requires ``cDd`` to also be set.

        lam : float or None
            Interporosity flow coefficient :math:`\\lambda` for naturally
            fractured reservoirs. Requires ``omega`` to also be set.

        omega : float or None
            Storativity ratio :math:`\\omega` for naturally fractured
            reservoirs. Requires ``lam`` to also be set.

        skin : float
            Skin factor ``S``. Default ``0.0``.

        **geom_kwargs
            Additional keyword arguments passed to ``geometry_fn``
            (e.g., ``reD``, ``xD``).

    Returns
    -------
        mpf
            Dimensionless wellbore pressure in Laplace space.

    Examples
    --------
    Infinite-acting radial flow with wellbore storage::

        >>> pwd = lambda s: lap_transform(s, radial_inf, cD=1e3)

    Naturally fractured, bounded radial flow::

        >>> pwd = lambda s: lap_transform(
        ...     s, radial_bounded_nf, reD=1000, lam=1e-6, omega=0.1, cD=100)

    Rate from bounded linear flow::

        >>> from gwr_inversion import gwr
        >>> pwd_fn = lambda s: lap_transform(s, linear_bounded_nf, xD=0.0)
        >>> qD = gwr(lambda s: 1 / (s**2 * pwd_fn(s)), time, M=8)
    """
    if (lam is None) != (omega is None):
        raise ValueError("dual porosity needs both lam and omega")
    if lam is not None and omega is not None and not (lam > 0 and 0 < omega < 1):
        raise ValueError(f"dual porosity needs lam > 0 and 0 < omega < 1, got {lam!r}, {omega!r}")
    if (cDd is None) != (aD is None):
        raise ValueError("variable wellbore storage needs both cDd and aD")
    if cDd is not None and not cD > 0:
        raise ValueError(f"variable wellbore storage needs cD > 0, got {cD!r}")
    if not cD >= 0:  # also rejects NaN, which would otherwise read as "no storage"
        raise ValueError(f"cD must be >= 0, got {cD!r}")
    return compose(s, geometry_fn, WellboreEffects(skin, cD, cDd, aD), lam, omega, geom_kwargs)


def compose(
    s: Any,
    geometry_fn: Callable[..., Any],
    wellbore: WellboreEffects,
    lam: Optional[float],
    omega: Optional[float],
    geometry_parameters: Dict[str, Any],
) -> Any:
    """
    ``lap_transform`` with the wellbore effects bundled and the geometry
    parameters as a dict; ``models.Model`` builds the bundle once per model.
    """
    # radial_inf and the linear geometries run in ARB ball arithmetic when
    # python-flint is installed (same quantity); others use the mpf code.
    if _flint.HAS_FLINT:
        ball_geometry = _ARB_GEOMETRIES.get(geometry_fn)
        if ball_geometry is not None:
            result = _arb_transforms.lap_transform(
                s, ball_geometry, wellbore, lam, omega, geometry_parameters)
            if result is not None:
                return result

    # 1. Naturally fractured: s -> s_eff;  2. geometry;  3-4. skin, then storage
    pressure = geometry_fn(s, effective_s(s, lam, omega), **geometry_parameters)
    return wellbore.apply(pressure, s)


# ============================================================================
# Hybrid Rate Analysis (real-domain)
# ============================================================================


def hybrid_rate(
    qD: "np.ndarray[Any, Any]",
    time: "np.ndarray[Any, Any]",
    b: float = 0.5,
) -> Tuple[
    "np.ndarray[Any, Any]",
    "np.ndarray[Any, Any]",
    "np.ndarray[Any, Any]",
    "np.ndarray[Any, Any]",
    "np.ndarray[Any, Any]",
]:
    """
    Hybrid rate-decline analysis from GWR-inverted dimensionless rate.

    Combines the numerically-inverted rate with an Arps hyperbolic
    correction to stabilize the solution at early time, then computes
    the standard suite of rate-transient analysis diagnostic functions.

    The correction is written in ``t_D`` with a unit time scale, whatever the
    geometry's reference length. With ``b = 0.5`` it is 5% of ``q_D`` at
    ``t_D = 10`` and 0.1% at ``t_D = 100`` for radial flow in ``r_w`` units
    (``reD = 1000``), so ``q_{Ddr}`` there differs from the uncorrected rate
    (1.91 vs 2.40 at ``t_D = 10``) until about ``t_D = 100``.

    Fulford, D. S., and Blasingame, T. A. 2013. Evaluation of Time-Rate
    Performance of Shale Wells using the Transient Hyperbolic Relation.
    Presented at SPE Unconventional Resources Conference – Canada in Calgary,
    Alberta, Canada, 5–7 November. SPE-167242-MS.
    https://doi.org/10.2118/167242-MS.

    The Arps correction is:

    .. math::

        q_{D,\\mathrm{hyp}} = \\frac{1}{(1 + t_D \\, b)^{1/b}} - e^{-t_D}

    The hybrid rate is:

    .. math::

        q_{D,\\mathrm{hybrid}} = q_D + q_{D,\\mathrm{hyp}}

    Diagnostic functions (all in the real domain):

    .. math::

        q_{Dd} = t_D \\left| \\frac{dq_{D,\\mathrm{hybrid}}}{dt_D} \\right|

    .. math::

        q_{Ddr} = \\frac{1}{t_D \\, \\frac{d}{dt_D}\\left(\\frac{1}{q_{D,\\mathrm{hybrid}}}\\right)}

    .. math::

        q_{Ddi} = \\frac{1}{t_D} \\int_0^{t_D} q_{D,\\mathrm{hybrid}} \\, d\\tau

    .. math::

        q_{Ddid} = q_{Ddi} - q_{D,\\mathrm{hybrid}}

    Parameters
    ----------
        qD : ndarray
            Dimensionless rate from GWR inversion, shape ``(n,)``.
            Typically obtained via
            ``gwr(lambda s: 1/(s**2 * lap_transform(s, ...)), time, M)``.

        time : ndarray
            Dimensionless time array, shape ``(n,)``. Must be positive and
            monotonically increasing.

        b : float
            Arps hyperbolic exponent for the early-time correction.
            Default ``0.5``.

    Returns
    -------
        qD_hybrid : ndarray
            Corrected hybrid dimensionless rate.

        qDd : ndarray
            Rate derivative, :math:`t_D |dq_D/dt_D|`.

        qDdr : ndarray
            Reciprocal rate derivative, :math:`1 / (t_D \\, d(1/q_D)/dt_D)`.

        qDdi : ndarray
            Cumulative-normalized rate (integral-average rate). The integral
            over ``[0, t_0]`` extrapolates ``q_D ~ t^n`` through the first two
            samples, so the first samples carry that extrapolation error
            (measured 1-16% at ``t_0``, falling off as ``t_0 / t``).

        qDdid : ndarray
            Integral derivative, :math:`q_{Ddi} - q_{D,\\mathrm{hybrid}}`.

    Raises
    ------
        ValueError
            For arrays of different shapes or fewer than two points, a time
            grid that is not positive and strictly increasing, or first two
            rates that are not positive or fall faster than ``1/t`` (the
            ``[0, t_0]`` extrapolation needs ``n > -1``).

        ImportError
            When scipy is not installed (``pip install gwr-inversion[hybrid]``).
    """
    try:
        from scipy.integrate import simpson
    except ImportError as error:
        raise ImportError(
            "hybrid_rate needs scipy; install with: pip install gwr-inversion[hybrid]") from error

    if not b > 0:
        raise ValueError(f"hybrid_rate needs the Arps exponent b > 0, got {b!r}")
    qD = np.asarray(qD, dtype=np.float64)
    time = np.asarray(time, dtype=np.float64)
    if qD.ndim != 1 or qD.shape != time.shape or len(time) < 2:
        raise ValueError("hybrid_rate needs qD and time as 1-D arrays of the same length >= 2")
    if not (np.all(time > 0) and np.all(np.diff(time) > 0)):
        raise ValueError("hybrid_rate needs time positive and strictly increasing")

    # Arps correction
    qD_hyp = 1 / (1 + time * b) ** (1 / b) - np.exp(-time)
    qD_hybrid = qD + qD_hyp

    # Rate derivative: tD * |dqD/dtD|
    qDd = time * np.abs(np.gradient(qD_hybrid, time))

    # Reciprocal rate derivative: 1 / (tD * d(1/qD)/dtD)
    qDdr = 1 / (time * np.gradient(1 / qD_hybrid, time))

    # Cumulative-normalized rate: (∫₀ᵗ qD dτ) / t. Before the first sample, qD is
    # extrapolated as a power law qD ~ t^exponent through the first two samples, so
    # ∫₀^t0 qD dτ = qD(t0) t0 / (exponent + 1); Simpson's rule covers t0 to t.
    first_rates = qD_hybrid[:2]
    if not np.all(np.isfinite(first_rates) & (first_rates > 0)):
        raise ValueError(f"hybrid_rate needs the first two hybrid rates positive and finite, "
                         f"got {first_rates.tolist()}")
    power_law_exponent = float(np.log(first_rates[1] / first_rates[0]) / np.log(time[1] / time[0]))
    if not power_law_exponent > -1:
        raise ValueError(
            f"the rate falls faster than 1/t between the first two samples (power-law exponent "
            f"{power_law_exponent:.3g} <= -1), so the integral from t = 0 cannot be "
            f"extrapolated; start the time grid earlier, in transient flow")
    first_interval_integral = qD_hybrid[0] * time[0] / (power_law_exponent + 1)
    sampled_integral = np.array([0.0] + [
        simpson(qD_hybrid[:i + 1], x=time[:i + 1]) for i in range(1, len(time))
    ])
    qDdi = (first_interval_integral + sampled_integral) / time

    # Integral derivative
    qDdid = qDdi - qD_hybrid

    return qD_hybrid, qDd, qDdr, qDdi, qDdid
