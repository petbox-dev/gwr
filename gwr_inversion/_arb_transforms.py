"""
Ball-arithmetic (ARB) evaluation of selected built-in Laplace transforms.

For ``radial_inf`` and the three linear geometries, ``lap_transform``
converts ``s`` to an exact ARB ball once, computes dual porosity, the
geometry, skin and wellbore storage in ball arithmetic, and converts the
result to mpf once. That removes the mpmath arithmetic and the per-call
conversions that dominate these cheap transforms (measured 1.2-1.9x at M=8,
1.2-1.4x at M=32 for radial_inf with WBS / dual porosity and linear no-flow).

The bounded radial and fractured-well geometries stay on the mpf code: in
ball arithmetic they measured 0.90-1.11x, because their cost is the Bessel
and integral evaluations themselves, which both paths share.

The order-0/1 Bessel set has its own precision loop and rebuilds its
argument ``sqrt(s_eff)`` from the exact ``s`` at each precision, so a retry
shrinks the argument's radius too. A final check on the result ball raises
the working precision of the whole chain when needed. Dual porosity, skin and
wellbore storage are the same ``_modifiers`` code the mpf path calls, so
both paths compute the same quantity.

Reached only through ``transforms.compose`` when ``_flint.HAS_FLINT``.
"""
from __future__ import annotations

from typing import Any, Callable, Dict, Optional

from mpmath import mp

from . import _flint
from ._modifiers import WellboreEffects, effective_s
from ._numeric import INITIAL_GUARD_BITS, MAX_ATTEMPTS, TARGET_GUARD_BITS, positive_real

# A callable that rebuilds a ball-valued quantity (s_eff) from the exact s at
# the current context precision.
BallFactory = Callable[[], Any]
# A ball geometry: (s ball, s_eff factory, **geometry parameters) -> sandface pressure ball.
BallGeometry = Callable[..., Any]


def _order_01(argument: BallFactory) -> _flint.BesselOrder01:
    """
    (K_0, K_1, I_0, I_1) balls at ``argument()``, raising the context
    precision (and rebuilding the argument at that precision) until each has
    ``mp.prec + TARGET_GUARD_BITS`` relative bits. Starts at the caller's
    context precision.
    """
    target_bits = mp.prec + TARGET_GUARD_BITS
    working_bits = _flint.flint_context.prec
    values = _flint.bessel_order_01_arb(argument())
    for _ in range(MAX_ATTEMPTS):
        accuracy_bits = min(value.rel_accuracy_bits() for value in values)
        if accuracy_bits >= target_bits:
            break
        working_bits = _flint.next_working_bits(working_bits, target_bits, accuracy_bits)
        with _flint.flint_context.workprec(working_bits):
            values = _flint.bessel_order_01_arb(argument())
    return values


# ---------------------------------------------------------------------------
# Geometries: (s, s_eff, **kwargs) -> sandface pressure ball, the same
# formulas as the mpf functions of the same name in ``transforms``.
# ``s`` is the exact ball; ``s_eff`` rebuilds s_eff at the current precision.
# ---------------------------------------------------------------------------


def radial_inf(s: Any, s_eff: BallFactory) -> Any:
    """``transforms.radial_inf``: ``K_0(sqrt(s_eff)) / (sqrt(s_eff) K_1(sqrt(s_eff)) s)``."""
    k0, k1, _, _ = _order_01(lambda: s_eff().sqrt())
    sqrt_s = s_eff().sqrt()
    return k0 / (sqrt_s * k1) / s


def linear_inf(s: Any, s_eff: BallFactory, xD: float = 0.0) -> Any:
    """``transforms.linear_inf``: ``exp(-sqrt(s_eff) x_D) / (s sqrt(s_eff))``."""
    sqrt_s = s_eff().sqrt()
    return (-sqrt_s * xD).exp() / (s * sqrt_s)


def linear_bounded_nf(s: Any, s_eff: BallFactory, xD: float = 0.0) -> Any:
    """``transforms.linear_bounded_nf``: ``cosh(sqrt(s_eff)(1 - x_D)) /
    (s sqrt(s_eff) sinh(sqrt(s_eff)))``."""
    sqrt_s = s_eff().sqrt()
    return (sqrt_s * (1 - xD)).cosh() / (s * sqrt_s * sqrt_s.sinh())


def linear_bounded_cp(s: Any, s_eff: BallFactory, xD: float = 0.0) -> Any:
    """``transforms.linear_bounded_cp``: ``sinh(sqrt(s_eff)(1 - x_D)) /
    (s sqrt(s_eff) cosh(sqrt(s_eff)))``."""
    sqrt_s = s_eff().sqrt()
    return (sqrt_s * (1 - xD)).sinh() / (s * sqrt_s * sqrt_s.cosh())


# ---------------------------------------------------------------------------
# Compositor
# ---------------------------------------------------------------------------


def lap_transform(
    s: Any,
    ball_geometry: BallGeometry,
    wellbore: WellboreEffects,
    lam: Optional[float],
    omega: Optional[float],
    geometry_parameters: Dict[str, Any],
) -> Optional[Any]:
    """
    ``transforms.compose`` in ball arithmetic. Returns the mpf result, or None
    (``s`` outside ``positive_real``; loop did not converge) so that the
    caller falls back to the mpf code.
    """
    s_mpf = positive_real(s)
    if s_mpf is None:
        return None

    target_bits = mp.prec + TARGET_GUARD_BITS
    working_bits = target_bits + INITIAL_GUARD_BITS
    s_arb = _flint.mpf_to_arb(s_mpf)

    def s_eff() -> Any:
        return effective_s(s_arb, lam, omega)

    for _ in range(MAX_ATTEMPTS):
        try:
            with _flint.flint_context.workprec(working_bits):
                pressure = ball_geometry(s_arb, s_eff, **geometry_parameters)
                wellbore_pressure = wellbore.apply(pressure, s_arb)
            accuracy_bits = wellbore_pressure.rel_accuracy_bits()
        except (AttributeError, TypeError):
            # An mpf parameter times a ball collapses the ball to an mpf (mpmath
            # reads its _mpf_), which then lacks the ball methods; leave the call
            # to the mpf code, which accepts every numeric type. Checking here,
            # not before the loop, keeps the per-call cost off the float path.
            # (test_arb_compositor_evaluates_every_arb_case pins that float
            # parameters never take this exit.)
            return None
        if accuracy_bits >= target_bits:
            return _flint.arb_to_mpf(wellbore_pressure)
        working_bits = _flint.next_working_bits(working_bits, target_bits, accuracy_bits)
    return None
