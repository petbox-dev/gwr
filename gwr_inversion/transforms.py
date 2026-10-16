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

from typing import Any, Callable, Optional, Tuple

import numpy as np
from mpmath import mp

from ._bessel import (
    besselk as _besselk, besseli as _besseli,
    besselk0_integral as _besselk0_integral, besseli0_integral as _besseli0_integral,
)


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
    i1_reD = _besseli(1, sqrt_s * reD)
    k1_reD = _besselk(1, sqrt_s * reD)
    k1 = _besselk(1, sqrt_s)
    i1 = _besseli(1, sqrt_s)

    return (
        (k0 * i1_reD + i0 * k1_reD)
        / (sqrt_s * k1 * i1_reD - sqrt_s * i1 * k1_reD)
    ) / s


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
    k0_reD = _besselk(0, sqrt_s * reD)
    i0_reD = _besseli(0, sqrt_s * reD)
    k1 = _besselk(1, sqrt_s)
    i1 = _besseli(1, sqrt_s)

    return (
        (k0 * i0_reD - i0 * k0_reD)
        / (sqrt_s * k1 * i0_reD + sqrt_s * i1 * k0_reD)
    ) / s


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
    c_ = 1 / (2 * s * sqrt_s)

    k0_int_1 = _besselk0_integral(sqrt_s * (1 - xD))
    k0_int_2 = _besselk0_integral(sqrt_s * (1 + xD))

    return c_ * (k0_int_1 + k0_int_2)


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
    c_ = 1 / (2 * s * sqrt_s)

    # Infinite-acting component
    k0_int_1 = _besselk0_integral(sqrt_s * (1 - xD))
    k0_int_2 = _besselk0_integral(sqrt_s * (1 + xD))
    inf_part = k0_int_1 + k0_int_2

    # Bounded correction
    k1_reD = _besselk(1, sqrt_s * reD)
    i1_reD = _besseli(1, sqrt_s * reD)
    i0_int_1 = _besseli0_integral(sqrt_s * (1 - xD))
    i0_int_2 = _besseli0_integral(sqrt_s * (1 + xD))
    nf_part = k1_reD / i1_reD * (i0_int_1 + i0_int_2)

    return c_ * (inf_part + nf_part)


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
    c_ = 1 / (2 * s * sqrt_s)

    # Infinite-acting component
    k0_int_1 = _besselk0_integral(sqrt_s * (1 - xD))
    k0_int_2 = _besselk0_integral(sqrt_s * (1 + xD))
    inf_part = k0_int_1 + k0_int_2

    # Constant-pressure correction
    k0_reD = _besselk(0, sqrt_s * reD)
    i0_reD = _besseli(0, sqrt_s * reD)
    i0_int_1 = _besseli0_integral(sqrt_s * (1 - xD))
    i0_int_2 = _besseli0_integral(sqrt_s * (1 + xD))
    cp_part = k0_reD / i0_reD * (i0_int_1 + i0_int_2)

    return c_ * (inf_part - cp_part)


# ============================================================================
# Physics Modifiers
# ============================================================================


def apply_dual_porosity(s: Any, lam: float, omega: float) -> Any:
    """
    Naturally fractured (dual-porosity) interporosity flow modifier.

    Computes the effective Laplace variable ``s_eff = f(s) * s`` for a
    naturally fractured reservoir using the Warren and Root (1963)
    pseudo-steady-state interporosity flow model.

    Warren, J. E. and Root, P. J. 1963. The Behavior of Naturally Fractured
    Reservoirs. SPE Journal, 3 (3): 245–255. SPE-426-PA.
    https://doi.org/10.2118/426-PA.

    .. math::

        f(s) = \\frac{\\lambda + \\omega (1 - \\omega) s}{\\lambda + (1 - \\omega) s}

    .. math::

        s_\\mathrm{eff} = f(s) \\cdot s

    Parameters
    ----------
        s : mpf
            Laplace variable.

        lam : float
            Interporosity flow coefficient :math:`\\lambda`, characterizing
            matrix-to-fracture flow. Dimensionless, typically ``1e-9`` to
            ``1e-3``.

        omega : float
            Storativity ratio :math:`\\omega`, the fraction of total storage
            in the fracture system. Dimensionless, ``0 < omega < 1``.

    Returns
    -------
        mpf
            Effective Laplace variable ``s_eff``.
    """
    fs = (lam + omega * (1 - omega) * s) / (lam + (1 - omega) * s)
    return fs * s


def apply_skin(psd: Any, s: Any, skin: float) -> Any:
    """
    Apply skin effect to a dimensionless pressure solution.

    .. math::

        \\bar{p}_{wD} = \\bar{p}_D + \\frac{S}{s}

    Parameters
    ----------
        psd : mpf
            Dimensionless pressure in Laplace space (before skin).

        s : mpf
            Laplace variable.

        skin : float
            Skin factor ``S``. Dimensionless. Positive for damage,
            negative for stimulation.

    Returns
    -------
        mpf
            Dimensionless wellbore pressure in Laplace space.
    """
    return psd + skin / s


def apply_wbs(psd: Any, s: Any, cD: float) -> Any:
    """
    Apply wellbore storage to a dimensionless pressure solution.

    van Everdingen, A. F. and Hurst, W. 1949. The Application of the
    Laplace Transformation to Flow Problems in Reservoirs. Journal of
    Petroleum Technology, 1 (12): 305–324. SPE-949305-G.
    https://doi.org/10.2118/949305-G.

    .. math::

        \\bar{p}_{wD} = \\frac{1}{\\frac{1}{\\bar{p}_D} + s^2 C_D}

    Parameters
    ----------
        psd : mpf
            Dimensionless pressure in Laplace space (before WBS).

        s : mpf
            Laplace variable.

        cD : float
            Dimensionless wellbore storage coefficient ``C_D``.
            Must be >= 0.

    Returns
    -------
        mpf
            Dimensionless wellbore pressure in Laplace space.
    """
    return 1 / (1 / psd + s ** 2 * cD)


def apply_variable_wbs(psd: Any, s: Any, cD: float, cDd: float,
                       aD: float) -> Any:
    """
    Apply variable wellbore storage to a dimensionless pressure solution.

    Models a changing wellbore storage coefficient (e.g., phase segregation,
    changing compressibility).

    Hegeman, P. S., Hallford, D. L., and Joseph, J. A. 1993.
    Well-Test Analysis With Changing Wellbore Storage. SPE Formation
    Evaluation, 8 (3): 201–207. SPE-21829-PA.
    https://doi.org/10.2118/21829-PA.

    .. math::

        W(s) = 1 + s^2 C_D C_{Dd} \\left( \\frac{1}{s} - \\frac{1}{s + 1/a_D} \\right)

    .. math::

        \\bar{p}_{wD} = \\frac{W(s)}{\\frac{1}{\\bar{p}_D} + s^2 C_D}

    Parameters
    ----------
        psd : mpf
            Dimensionless pressure in Laplace space (before WBS).

        s : mpf
            Laplace variable.

        cD : float
            Dimensionless wellbore storage coefficient ``C_D``.

        cDd : float
            Dimensionless wellbore storage change ratio ``C_{Dd}``,
            defined as ``(C_Ds - C_DL) / C_DL``.

        aD : float
            Dimensionless wellbore storage change rate ``a_D``,
            controlling the transition speed.

    Returns
    -------
        mpf
            Dimensionless wellbore pressure in Laplace space.
    """
    wbs = 1 + s ** 2 * cD * cDd * (1 / s - 1 / (s + 1 / aD))
    return wbs / (1 / psd + s ** 2 * cD)


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
    Dimensionless cumulative production in Laplace space.

    .. math::

        \\bar{N}_{pD} = \\frac{\\bar{p}_{wD}}{s}

    Parameters
    ----------
        pwd : mpf
            Dimensionless wellbore pressure in Laplace space.

        s : mpf
            Laplace variable.

    Returns
    -------
        mpf
            Dimensionless cumulative production in Laplace space.
    """
    return pwd / s


def output_derivative(pwd: Any, s: Any) -> Any:
    """
    Bourdet pressure derivative in Laplace space.

    The Laplace-domain multiplication by ``s`` corresponds to
    multiplication by ``t`` and differentiation in the time domain,
    yielding the Bourdet derivative :math:`t \\, dp_{wD}/dt` after
    inversion.

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
            Bourdet derivative in Laplace space. Multiply by ``t``
            after inversion to obtain the time-domain derivative.
    """
    return s * pwd


# ============================================================================
# Compositor
# ============================================================================


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

        >>> from functools import partial
        >>> pwd = lambda s: lap_transform(s, radial_inf, cD=1e3)

    Naturally fractured, bounded radial flow::

        >>> pwd = lambda s: lap_transform(
        ...     s, radial_bounded_nf, reD=1000, lam=1e-6, omega=0.1, cD=100)

    Rate from bounded linear flow::

        >>> from gwr_inversion import gwr
        >>> pwd_fn = lambda s: lap_transform(s, linear_bounded_nf, xD=0.0)
        >>> qD = gwr(lambda s: 1 / (s**2 * pwd_fn(s)), time, M=8)
    """
    # 1. Naturally fractured: modify s → s_eff
    if lam is not None and omega is not None:
        s_eff = apply_dual_porosity(s, lam, omega)
    else:
        s_eff = s

    # 2. Geometry
    psd = geometry_fn(s, s_eff, **geom_kwargs)

    # 3. Skin
    if skin != 0:
        psd = apply_skin(psd, s, skin)

    # 4. Wellbore storage
    if cD > 0:
        if cDd is not None and aD is not None:
            pwd = apply_variable_wbs(psd, s, cD, cDd, aD)
        else:
            pwd = apply_wbs(psd, s, cD)
    else:
        pwd = psd

    return pwd


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
            Cumulative-normalized rate (integral-average rate).

        qDdid : ndarray
            Integral derivative, :math:`q_{Ddi} - q_{D,\\mathrm{hybrid}}`.
    """
    from scipy.integrate import simpson

    qD = np.asarray(qD, dtype=np.float64)
    time = np.asarray(time, dtype=np.float64)

    # Arps correction
    qD_hyp = 1 / (1 + time * b) ** (1 / b) - np.exp(-time)
    qD_hybrid = qD + qD_hyp

    # Rate derivative: tD * |dqD/dtD|
    qDd = time * np.abs(np.gradient(qD_hybrid, time))

    # Reciprocal rate derivative: 1 / (tD * d(1/qD)/dtD)
    qDdr = 1 / (time * np.gradient(1 / qD_hybrid, time))

    # Cumulative-normalized rate: (∫₀ᵗ qD dτ) / t
    qDdi = np.array([
        simpson(qD_hybrid[:i + 1], x=time[:i + 1]) for i in range(len(time))
    ]) / time

    # Power-law extrapolation for the first point (avoid 0/0)
    B = float(mp.log(qD_hybrid[1] / qD_hybrid[0]) / mp.log(time[1] / time[0]))
    A = qD_hybrid[0] / time[0] ** B
    qDdi[0] = A / (B + 1) * time[0] ** (B + 1)

    # Integral derivative
    qDdid = qDdi - qD_hybrid

    return qD_hybrid, qDd, qDdr, qDdi, qDdid
