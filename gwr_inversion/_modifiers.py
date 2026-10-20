"""
Physics modifiers shared by the mpf compositor (``transforms.lap_transform``),
the ball-arithmetic compositor (``_arb_transforms.lap_transform``) and the
generic ``Model`` path. They use only ``+ - * / **``, so the same functions
evaluate mpf and ARB balls.
"""
from __future__ import annotations

import math
import warnings
from typing import Any, NamedTuple, Optional


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
    interporosity_factor = (lam + omega * (1 - omega) * s) / (lam + (1 - omega) * s)
    return interporosity_factor * s


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
             = 1 + \\frac{C_D C_{Dd} \\, s}{1 + a_D s}

    The second form is evaluated: the first cancels at large ``s`` (measured
    relative error 4e-9 at s = 1e8 and 17 digits), the second does not.

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
    storage_change_factor = 1 + cD * cDd * s / (1 + aD * s)
    return storage_change_factor / (1 / psd + s ** 2 * cD)


def effective_s(s: Any, lam: Optional[float], omega: Optional[float]) -> Any:
    """
    Step 1 of the compositor: the effective Laplace variable. Dual porosity
    applies when both ``lam`` and ``omega`` are set; otherwise ``s`` itself.
    """
    if lam is not None and omega is not None:
        return apply_dual_porosity(s, lam, omega)
    return s


def effective_radius_scale(skin: float) -> float:
    """
    The factor ``rho^2 = e^{-2S}`` that turns a radial geometry into the same
    geometry seen from the effective wellbore radius ``r_w' = r_w e^{-S}``.

    With ``t_D`` and ``C_D`` defined on ``r_w``, the solution for a well of
    radius ``r_w'`` and no skin is ``p_D(t_D) = p_D'(t_D rho^2)`` with
    ``r_eD' = r_eD / rho``, so in Laplace space

    .. math::

        \\bar{p}_D(s) = \\rho^2 \\, \\bar{p}_D'(\\rho^2 s), \\qquad \\rho = e^{-S}.

    Dual porosity commutes with the scaling (``s_eff' = rho^2 s_eff``) and
    storage composes afterwards with ``C_D`` unchanged. Earlougher, R. C. Jr.
    1977. Advances in Well Test Analysis, SPE Monograph 5, section 2.3;
    Gringarten, A. C., Bourdet, D. P., Landel, P. A. and Kniazeff, V. J. 1979.
    A Comparison Between Different Skin and Wellbore Storage Type-Curves for
    Early-Time Transient Analysis. SPE-8205-MS. https://doi.org/10.2118/8205-MS.
    """
    return math.exp(-2.0 * skin)


def warn_if_unstable(skin: float, storage: float, effective_radius: bool = False) -> None:
    """
    Warn when the infinitesimal skin ``S / s`` with ``S < 0`` is composed with
    wellbore storage: ``1 + C_D s^2 p_S(s)`` then has a zero at some ``s* > 0``,
    the exact inverse grows as ``e^{s* t}``, and a numerical inversion returns
    neither that nor a physical curve. Use the effective-radius skin instead
    (docs/wellbore-storage-real-domain.md, section 7).
    """
    if skin < 0 and storage > 0 and not effective_radius:
        warnings.warn(
            f"skin {skin!r} < 0 composed with wellbore storage {storage!r} as S/s puts a "
            "pole of p_wD on the positive real axis; the result is not the inverse of the "
            "transform. Use the effective-radius skin (EffectiveRadiusSkin / "
            "skin_model='effective_radius') for a stimulated well with storage.",
            stacklevel=3)


class WellboreEffects(NamedTuple):
    """
    Skin and wellbore storage: steps 3 and 4 of the compositor, in that order.

    ``skin`` is S (0: no skin). ``storage`` is C_D (<= 0: no storage). Storage
    is variable when both ``storage_change_ratio`` (C_Dd) and ``change_rate``
    (a_D) are set, else constant. With ``effective_radius`` the skin is an
    effective wellbore radius applied to the geometry (step 2, through
    ``geometry_scale``) and step 3 is skipped.
    """

    skin: float = 0.0
    storage: float = 0.0
    storage_change_ratio: Optional[float] = None
    change_rate: Optional[float] = None
    effective_radius: bool = False

    def geometry_scale(self) -> float:
        """``rho^2 = e^{-2S}`` for an effective-radius skin, else 1."""
        if self.effective_radius and self.skin != 0:
            return effective_radius_scale(self.skin)
        return 1.0

    def apply(self, pressure: Any, s: Any) -> Any:
        """Wellbore pressure from the sandface ``pressure`` (mpf or ARB ball) at ``s``."""
        if self.skin != 0 and not self.effective_radius:
            pressure = apply_skin(pressure, s, self.skin)
        if self.storage > 0:
            if self.storage_change_ratio is not None and self.change_rate is not None:
                return apply_variable_wbs(pressure, s, self.storage,
                                          self.storage_change_ratio, self.change_rate)
            return apply_wbs(pressure, s, self.storage)
        return pressure
