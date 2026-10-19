"""
Numeric policy shared by the Bessel backends (``_flint``, ``_gmpy2_bessel``,
``_bessel``) and the ball-arithmetic transforms (``_arb_transforms``): the
argument domain the fast backends accept, and the guard bits of their
precision loops.
"""
from __future__ import annotations

import math
from typing import Any, Optional

from mpmath import mp

# log2(e). The K_n power series loses about 2 x log2(e) bits to cancellation at x,
# and K_n(x) / I_n(x) falls like 2^(-2 x log2(e)).
LOG2_E = math.log2(math.e)
TWO_LOG2_E = 2.0 * LOG2_E

# Bits of relative accuracy required above mp.prec.
TARGET_GUARD_BITS = 10
# Bits added to the first working precision above the target.
INITIAL_GUARD_BITS = 16
# Bits added on each retry beyond the measured deficit.
RETRY_GUARD_BITS = 32
# Retries of an adaptive precision loop before it gives up (the caller falls back).
MAX_ATTEMPTS = 8


def positive_real(x: Any) -> Optional[Any]:
    """
    ``x`` as an mpf (rounded to ``mp.prec``) when it is a finite real number
    greater than zero; None otherwise.

    The ARB and gmpy2 backends evaluate only this domain. None (complex ``x``,
    as on the Talbot contour; zero; negative; inf; nan) tells the caller to
    fall back to the next backend and finally to mpmath.
    """
    if isinstance(x, mp.mpf) and x._mpf_[3] <= mp.prec:
        value = x  # already at working precision: mp.mpf(x) would be an equal copy
    else:
        try:
            value = mp.mpf(x)
        except TypeError:
            return None
    sign, mantissa, _, _ = value._mpf_
    if sign or not mantissa:  # negative, zero, inf or nan
        return None
    return value
