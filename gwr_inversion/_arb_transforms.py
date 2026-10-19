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
the working precision of the whole chain when needed. The arithmetic follows
the mpf code in ``transforms.py`` operation by operation (Python-float
coefficients such as ``omega * (1 - omega)`` are formed the same way), so
both paths compute the same quantity.
"""
from __future__ import annotations

from typing import Any, Callable, Dict, Optional

from mpmath import mp

from . import _flint

if _flint.HAS_FLINT:
    from flint import ctx as _flint_ctx

# A ball-valued quantity rebuilt from the exact s at the current ctx.prec.
Argument = Callable[[], Any]
ArbGeometry = Callable[..., Any]


def _order_01(argument: Argument) -> Any:
    """
    (K_0, K_1, I_0, I_1) balls at ``argument()``, raising ctx.prec (and
    rebuilding the argument at that precision) until each has
    ``mp.prec + 10`` relative bits.
    """
    target_bits = mp.prec + _flint._TARGET_GUARD_BITS
    working_bits = _flint_ctx.prec
    values = _flint._bessel_order_01_arb(argument())
    for _ in range(_flint._MAX_ATTEMPTS):
        accuracy_bits = min(value.rel_accuracy_bits() for value in values)
        if accuracy_bits >= target_bits:
            break
        working_bits = _flint._next_working_bits(working_bits, target_bits, accuracy_bits)
        previous_bits = _flint_ctx.prec
        try:
            _flint_ctx.prec = working_bits
            values = _flint._bessel_order_01_arb(argument())
        finally:
            _flint_ctx.prec = previous_bits
    return values


# ---------------------------------------------------------------------------
# Geometries: (s, s_eff, **kwargs) -> psd ball.
# ``s`` is the exact ball; ``s_eff`` rebuilds s_eff at the current precision.
# ---------------------------------------------------------------------------


def radial_inf(s: Any, s_eff: Argument) -> Any:
    k0, k1, _, _ = _order_01(lambda: s_eff().sqrt())
    sqrt_s = s_eff().sqrt()
    return k0 / (sqrt_s * k1) / s


def linear_inf(s: Any, s_eff: Argument, xD: float = 0.0) -> Any:
    sqrt_s = s_eff().sqrt()
    return (-sqrt_s * xD).exp() / (s * sqrt_s)


def linear_bounded_nf(s: Any, s_eff: Argument, xD: float = 0.0) -> Any:
    sqrt_s = s_eff().sqrt()
    return (sqrt_s * (1 - xD)).cosh() / (s * sqrt_s * sqrt_s.sinh())


def linear_bounded_cp(s: Any, s_eff: Argument, xD: float = 0.0) -> Any:
    sqrt_s = s_eff().sqrt()
    return (sqrt_s * (1 - xD)).sinh() / (s * sqrt_s * sqrt_s.cosh())


# ---------------------------------------------------------------------------
# Compositor
# ---------------------------------------------------------------------------


def lap_transform(
    s: Any,
    arb_geometry: ArbGeometry,
    cD: float,
    cDd: Optional[float],
    aD: Optional[float],
    lam: Optional[float],
    omega: Optional[float],
    skin: float,
    geom_kwargs: Dict[str, Any],
) -> Optional[Any]:
    """
    ``transforms.lap_transform`` in ball arithmetic. Returns the mpf result,
    or None (positive real ``s`` required; loop did not converge) so that the
    caller falls back to the mpf code.
    """
    try:
        s_mpf = mp.mpf(s)
    except TypeError:
        return None
    if not s_mpf > 0:
        return None

    target_bits = mp.prec + _flint._TARGET_GUARD_BITS
    working_bits = target_bits + _flint._INITIAL_GUARD_BITS
    previous_bits = _flint_ctx.prec
    try:
        _flint_ctx.prec = mp.prec + 64
        s_arb = _flint._mpf_to_arb(s_mpf)

        def s_eff() -> Any:
            if lam is not None and omega is not None:
                fs = (lam + omega * (1 - omega) * s_arb) / (lam + (1 - omega) * s_arb)
                return fs * s_arb
            return s_arb

        for _ in range(_flint._MAX_ATTEMPTS):
            _flint_ctx.prec = working_bits
            psd = arb_geometry(s_arb, s_eff, **geom_kwargs)
            # The geometry may have changed ctx.prec (Bessel precision loop).
            _flint_ctx.prec = working_bits
            if skin != 0:
                psd = psd + skin / s_arb
            if cD > 0:
                if cDd is not None and aD is not None:
                    wbs = 1 + s_arb ** 2 * cD * cDd * (1 / s_arb - 1 / (s_arb + 1 / aD))
                    pwd = wbs / (1 / psd + s_arb ** 2 * cD)
                else:
                    pwd = 1 / (1 / psd + s_arb ** 2 * cD)
            else:
                pwd = psd
            accuracy_bits = pwd.rel_accuracy_bits()
            if accuracy_bits >= target_bits:
                return _flint._arb_to_mpf(pwd)
            working_bits = _flint._next_working_bits(working_bits, target_bits, accuracy_bits)
        return None
    finally:
        _flint_ctx.prec = previous_bits
