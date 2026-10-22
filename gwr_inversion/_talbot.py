"""Fixed Talbot method via mpmath.invertlaplace."""

from collections.abc import Callable
from typing import Any

from mpmath import invertlaplace


def talbot_single(fn: Callable[..., Any], time: float, degree: int = 32) -> Any:
    """Invert fn at a single time point using Fixed Talbot."""
    return invertlaplace(fn, time, method='talbot', degree=degree)
