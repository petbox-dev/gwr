"""
Generate diagnostic figures for gwr_inversion.transforms.

Produces 10 figures:
  - 6 boundary-comparison figures (pressure + rate per geometry)
  - 3 modifier figures (one per geometry, 4 subplots each)
  - 1 hybrid rate analysis figure (3 subplots, no-flow boundary)

Also validates all solutions by comparing M=8 vs M=16.

Usage:
    PYTHONPATH=. python test/generate_figures.py
"""

import time as timer
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from mpmath import mp

from gwr_inversion import (
    frac_bounded_cp,
    frac_bounded_nf,
    frac_inf,
    gwr,
    hybrid_rate,
    lap_transform,
    linear_bounded_cp,
    linear_bounded_nf,
    linear_inf,
    radial_bounded_cp,
    radial_bounded_nf,
    radial_inf,
)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OUTDIR = Path("test/figures")
DPS = 50

# Time arrays per geometry
RADIAL_TIME = 10 ** np.linspace(-2, 6, 81)
LINEAR_TIME = 10 ** np.linspace(-3, 3, 61)
FRAC_TIME = 10 ** np.linspace(-2, 3, 41)  # fewer points (quadrature is expensive)

# Geometry parameters
RADIAL_RED = 1000.0
FRAC_RED = 100.0

# Modifier parameters
SKIN = 3.0
CD_WBS = 1e3
CD_VWBS, CDD_VWBS, AD_VWBS = 1e3, -5.75, 1e7
LAM_NF, OMEGA_NF = 1e-6, 0.1

# Hybrid rate Arps b
ARPS_B = 0.5

# ---------------------------------------------------------------------------
# Plot limits: (xlim_lo, xlim_hi, ylim_lo, ylim_hi)
# Set ylim_lo = None to auto-compute from xlim span (square aspect).
# ---------------------------------------------------------------------------

# Boundary comparison — pressure (3 subplots: inf, nf, cp)
RADIAL_PRESSURE_LIMITS = [
    (1e-2, 1e6, None, None),   # inf
    (1e-2, 1e6, None, None),   # nf
    (1e-2, 1e6, None, None),   # cp
]
RADIAL_RATE_LIMITS = [
    (1e-2, 1e6, None, None),
    (1e-2, 1e6, None, None),
    (1e-2, 1e6, None, None),
]
LINEAR_PRESSURE_LIMITS = [
    (1e-3, 1e3, None, None),
    (1e-3, 1e3, None, None),
    (1e-3, 1e3, None, None),
]
LINEAR_RATE_LIMITS = [
    (1e-3, 1e3, None, None),
    (1e-3, 1e3, None, None),
    (1e-3, 1e3, None, None),
]
FRAC_PRESSURE_LIMITS = [
    (1e-2, 1e3, None, None),
    (1e-2, 1e3, None, None),
    (1e-2, 1e3, None, None),
]
FRAC_RATE_LIMITS = [
    (1e-2, 1e3, None, None),
    (1e-2, 1e3, None, None),
    (1e-2, 1e3, None, None),
]

# Modifier figures (4 subplots: skin, WBS, vWBS, dual-porosity)
RADIAL_MODIFIER_LIMITS = [
    (1e-2, 1e6, None, None),
    (1e-2, 1e6, None, None),
    (1e-2, 1e6, None, None),
    (1e-2, 1e6, None, None),
]
LINEAR_MODIFIER_LIMITS = [
    (1e-3, 1e3, None, None),
    (1e-3, 1e3, None, None),
    (1e-3, 1e3, None, None),
    (1e-3, 1e3, None, None),
]
FRAC_MODIFIER_LIMITS = [
    (1e-2, 1e3, None, None),
    (1e-2, 1e3, None, None),
    (1e-2, 1e3, None, None),
    (1e-2, 1e3, None, None),
]

# Hybrid rate analysis (3 subplots: radial, linear, frac)
HYBRID_LIMITS = [
    (1e4, 1e6, None, None),    # radial
    (1e-3, 1e3, None, None),   # linear
    (1e-2, 1e3, None, None),   # frac
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


REFINE = True


def _pressure(geometry_fn: Callable, time: np.ndarray, M: int, **kwargs: Any) -> np.ndarray:
    """Compute dimensionless pressure via GWR."""
    return gwr(lambda s: lap_transform(s, geometry_fn, **kwargs), time, M, refine=REFINE)


def _derivative(geometry_fn: Callable, time: np.ndarray, M: int, **kwargs: Any) -> np.ndarray:
    """Compute Bourdet derivative (t * dpD/dt) via GWR."""
    return time * gwr(lambda s: s * lap_transform(s, geometry_fn, **kwargs), time, M, refine=REFINE)


def _rate(geometry_fn: Callable, time: np.ndarray, M: int, **kwargs: Any) -> np.ndarray:
    """Compute dimensionless rate via GWR."""
    return gwr(lambda s: 1 / (s**2 * lap_transform(s, geometry_fn, **kwargs)), time, M, refine=REFINE)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def validate_case(
    label: str,
    fn: Callable[[Any], Any],
    time: np.ndarray,
    M_low: int = 8,
    M_high: int = 16,
) -> float:
    """Compare M_low vs M_high and return max relative error.

    Uses a combined relative/absolute tolerance: the error metric at each
    point is |low - high| / max(|high|, atol), which avoids false positives
    when the solution passes through zero (e.g., rate in bounded systems
    at late time).
    """
    r_low = gwr(fn, time, M_low)
    r_high = gwr(fn, time, M_high)

    # Scale-aware error: use fraction of peak magnitude as absolute tolerance
    # so near-zero values (e.g., rate crossing zero in bounded systems) don't
    # dominate the error metric
    peak = float(np.max(np.abs(r_high)))
    atol = peak * 1e-6 if peak > 0 else 1e-30
    denom = np.maximum(np.abs(r_high), atol)
    rel_err = float(np.max(np.abs(r_low - r_high) / denom))

    status = "OK" if rel_err < 1e-6 else "CHECK" if rel_err < 1e-3 else "WARN"
    print(f"  [{status}] {label:<55} max_rel_err = {rel_err:.2e}")
    return rel_err


def run_all_validation() -> Dict[str, float]:
    """Validate all cases at M=8 vs M=16."""
    print("=" * 78)
    print("  Validation: M=8 vs M=16")
    print("=" * 78)
    results: Dict[str, float] = {}

    # --- Radial ---
    print("\nRadial geometry:")
    cases: List[Tuple[str, Callable, np.ndarray]] = [
        ("radial_inf pressure", lambda s: lap_transform(s, radial_inf), RADIAL_TIME),
        ("radial_inf rate", lambda s: 1 / (s**2 * lap_transform(s, radial_inf)), RADIAL_TIME),
        (
            "radial_bounded_nf pressure",
            lambda s: lap_transform(s, radial_bounded_nf, reD=RADIAL_RED),
            RADIAL_TIME,
        ),
        (
            "radial_bounded_nf rate",
            lambda s: 1 / (s**2 * lap_transform(s, radial_bounded_nf, reD=RADIAL_RED)),
            RADIAL_TIME,
        ),
        (
            "radial_bounded_cp pressure",
            lambda s: lap_transform(s, radial_bounded_cp, reD=RADIAL_RED),
            RADIAL_TIME,
        ),
        (
            "radial_bounded_cp rate",
            lambda s: 1 / (s**2 * lap_transform(s, radial_bounded_cp, reD=RADIAL_RED)),
            RADIAL_TIME,
        ),
    ]
    for label, fn, t in cases:
        results[label] = validate_case(label, fn, t)

    # --- Linear ---
    print("\nLinear geometry:")
    cases = [
        ("linear_inf pressure", lambda s: lap_transform(s, linear_inf, xD=0.0), LINEAR_TIME),
        (
            "linear_inf rate",
            lambda s: 1 / (s**2 * lap_transform(s, linear_inf, xD=0.0)),
            LINEAR_TIME,
        ),
        (
            "linear_bounded_nf pressure",
            lambda s: lap_transform(s, linear_bounded_nf, xD=0.0),
            LINEAR_TIME,
        ),
        (
            "linear_bounded_nf rate",
            lambda s: 1 / (s**2 * lap_transform(s, linear_bounded_nf, xD=0.0)),
            LINEAR_TIME,
        ),
        (
            "linear_bounded_cp pressure",
            lambda s: lap_transform(s, linear_bounded_cp, xD=0.0),
            LINEAR_TIME,
        ),
        (
            "linear_bounded_cp rate",
            lambda s: 1 / (s**2 * lap_transform(s, linear_bounded_cp, xD=0.0)),
            LINEAR_TIME,
        ),
    ]
    for label, fn, t in cases:
        results[label] = validate_case(label, fn, t)

    # --- Fractured ---
    print("\nFractured geometry:")
    cases = [
        ("frac_inf pressure", lambda s: lap_transform(s, frac_inf), FRAC_TIME),
        ("frac_inf rate", lambda s: 1 / (s**2 * lap_transform(s, frac_inf)), FRAC_TIME),
        (
            "frac_bounded_nf pressure",
            lambda s: lap_transform(s, frac_bounded_nf, reD=FRAC_RED),
            FRAC_TIME,
        ),
        (
            "frac_bounded_nf rate",
            lambda s: 1 / (s**2 * lap_transform(s, frac_bounded_nf, reD=FRAC_RED)),
            FRAC_TIME,
        ),
        (
            "frac_bounded_cp pressure",
            lambda s: lap_transform(s, frac_bounded_cp, reD=FRAC_RED),
            FRAC_TIME,
        ),
        (
            "frac_bounded_cp rate",
            lambda s: 1 / (s**2 * lap_transform(s, frac_bounded_cp, reD=FRAC_RED)),
            FRAC_TIME,
        ),
    ]
    for label, fn, t in cases:
        results[label] = validate_case(label, fn, t)

    # --- Modifiers (radial_inf only, representative) ---
    print("\nModifiers (radial_inf):")
    cases = [
        ("radial_inf + skin", lambda s: lap_transform(s, radial_inf, skin=SKIN), RADIAL_TIME),
        ("radial_inf + WBS", lambda s: lap_transform(s, radial_inf, cD=CD_WBS), RADIAL_TIME),
        (
            "radial_inf + variable WBS",
            lambda s: lap_transform(s, radial_inf, cD=CD_VWBS, cDd=CDD_VWBS, aD=AD_VWBS),
            RADIAL_TIME,
        ),
        (
            "radial_inf + dual porosity",
            lambda s: lap_transform(s, radial_inf, lam=LAM_NF, omega=OMEGA_NF),
            RADIAL_TIME,
        ),
    ]
    for label, fn, t in cases:
        results[label] = validate_case(label, fn, t)

    # Summary
    print(f"\n{'-' * 78}")
    worst = max(results.values())
    n_warn = sum(1 for v in results.values() if v >= 1e-3)
    n_check = sum(1 for v in results.values() if 1e-6 <= v < 1e-3)
    print(f"  {len(results)} cases validated. Worst: {worst:.2e}. WARN: {n_warn}, CHECK: {n_check}")
    print(f"{'-' * 78}\n")

    return results


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------


def _setup_loglog(
    ax: Axes,
    title: str,
    xlabel: str,
    ylabel: str,
    limits: Optional[Tuple[float, float, Optional[float], Optional[float]]] = None,
) -> None:
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(True, which="major", ls=":", alpha=0.6)

    if limits is not None:
        xl, xh, yl, yh = limits
        ax.set_xlim(xl, xh)
        n_decades = np.log10(xh) - np.log10(xl)
        if yh is None:
            yh = ax.get_ylim()[1]
        if yl is None:
            yl = yh / 10 ** n_decades
        ax.set_ylim(yl, yh)
    else:
        xlim = ax.get_xlim()
        n_decades = np.log10(xlim[1]) - np.log10(xlim[0])
        ylim_upper = ax.get_ylim()[1]
        ylim_lower = ylim_upper / 10 ** n_decades
        ax.set_ylim(ylim_lower, ylim_upper)

    ax.set_aspect("equal")


# ---------------------------------------------------------------------------
# Time padding
# ---------------------------------------------------------------------------


def _padded_time(
    time: np.ndarray,
    pad_before: float = 2.0,
    pad_after: float = 1.0,
) -> Tuple[np.ndarray, int, int]:
    """Create a time array padded before and after the display range.

    Early padding stabilizes derivatives/integrals at the first display point.
    Late padding extends the computation beyond the plot limits.

    Returns (padded_array, n_pad_before, n_display).
    """
    log_min = np.log10(time[0])
    log_max = np.log10(time[-1])
    n_display = len(time)
    ppd = (n_display - 1) / (log_max - log_min)
    n_before = int(round(pad_before * ppd))
    n_after = int(round(pad_after * ppd))
    total = n_before + n_display + n_after
    padded = 10 ** np.linspace(log_min - pad_before, log_max + pad_after, total)
    return padded, n_before, n_display


# ---------------------------------------------------------------------------
# Figure 1-6: Boundary comparison (pressure + rate per geometry)
# ---------------------------------------------------------------------------


def _boundary_figures(
    geom_name: str,
    geom_inf: Callable,
    geom_nf: Callable,
    geom_cp: Callable,
    display_time: np.ndarray,
    M: int,
    geom_kwargs_nf: dict,
    geom_kwargs_cp: dict,
    geom_kwargs_inf: Optional[dict] = None,
    pressure_limits: Optional[List[Tuple[float, float, Optional[float], Optional[float]]]] = None,
    rate_limits: Optional[List[Tuple[float, float, Optional[float], Optional[float]]]] = None,
) -> Tuple[plt.Figure, plt.Figure]:
    """Generate pressure and rate boundary comparison figures for one geometry."""
    kw_inf = geom_kwargs_inf or {}

    # Pad 1 decade before and after for edge stability
    full_time, n_before, n_display = _padded_time(display_time, pad_before=1.0, pad_after=1.0)
    sl = slice(n_before, n_before + n_display)

    # --- Pressure figure ---
    fig_p, axes_p = plt.subplots(1, 3, figsize=(16, 5))
    fig_p.suptitle(f"{geom_name} — Pressure & Derivative", fontsize=12)

    for i, (ax, (label, geom, kw)) in enumerate(zip(
        axes_p,
        [
            ("Infinite-Acting", geom_inf, kw_inf),
            ("No-Flow Boundary", geom_nf, geom_kwargs_nf),
            ("Constant-Pressure Boundary", geom_cp, geom_kwargs_cp),
        ],
    )):
        pD = _pressure(geom, full_time, M, **kw)
        dpD = _derivative(geom, full_time, M, **kw)
        t = full_time[sl]
        ax.loglog(t, pD[sl], "o-", ms=2, lw=0.8, label=r"$p_D$")
        ax.loglog(t, np.abs(dpD[sl]), "s-", ms=2, lw=0.8, label=r"$t_D \, p_D'$")
        lim = pressure_limits[i] if pressure_limits else None
        _setup_loglog(ax, label, r"$t_D$", r"$p_D$", limits=lim)

    fig_p.tight_layout()

    # --- Rate figure ---
    fig_q, axes_q = plt.subplots(1, 3, figsize=(16, 5))
    fig_q.suptitle(f"{geom_name} — Rate & Derivative", fontsize=12)

    for i, (ax, (label, geom, kw)) in enumerate(zip(
        axes_q,
        [
            ("Infinite-Acting", geom_inf, kw_inf),
            ("No-Flow Boundary", geom_nf, geom_kwargs_nf),
            ("Constant-Pressure Boundary", geom_cp, geom_kwargs_cp),
        ],
    )):
        qD = _rate(geom, full_time, M, **kw)
        dqD = full_time * np.abs(np.gradient(qD, full_time))
        t = full_time[sl]
        ax.loglog(t, np.abs(qD[sl]), "o-", ms=2, lw=0.8, label=r"$q_D$")
        ax.loglog(t, dqD[sl], "s-", ms=2, lw=0.8, label=r"$t_D \, |q_D'|$")
        lim = rate_limits[i] if rate_limits else None
        _setup_loglog(ax, label, r"$t_D$", r"$q_D$", limits=lim)

    fig_q.tight_layout()
    return fig_p, fig_q


# ---------------------------------------------------------------------------
# Figure 7-9: Modifier figures (one per geometry, 4 subplots)
# ---------------------------------------------------------------------------


def _modifier_figure(
    geom_name: str,
    geom_inf: Callable,
    display_time: np.ndarray,
    M: int,
    geom_kwargs_inf: Optional[dict] = None,
    modifier_limits: Optional[List[Tuple[float, float, Optional[float], Optional[float]]]] = None,
) -> plt.Figure:
    """Generate modifier comparison figure for one geometry (infinite-acting)."""
    kw = geom_kwargs_inf or {}
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"{geom_name} — Physics Modifiers (Infinite-Acting)", fontsize=12)

    # Pad 1 decade before and after
    full_time, n_before, n_display = _padded_time(display_time, pad_before=1.0, pad_after=1.0)
    sl = slice(n_before, n_before + n_display)
    t = full_time[sl]

    # Base case
    pD_base = _pressure(geom_inf, full_time, M, **kw)
    dpD_base = _derivative(geom_inf, full_time, M, **kw)

    modifiers = [
        (f"Skin (S={SKIN})", dict(skin=SKIN)),
        (f"WBS ($C_D$={CD_WBS:.0e})", dict(cD=CD_WBS)),
        (f"Variable WBS ($C_D$={CD_VWBS:.0e})", dict(cD=CD_VWBS, cDd=CDD_VWBS, aD=AD_VWBS)),
        (r"Dual Porosity ($\lambda$=1e-6, $\omega$=0.1)", dict(lam=LAM_NF, omega=OMEGA_NF)),
    ]

    for i, (ax, (mod_label, mod_kw)) in enumerate(zip(axes.flat, modifiers)):
        # Base
        ax.loglog(t, pD_base[sl], "-", lw=0.8, color="gray", alpha=0.6, label=r"$p_D$ (base)")
        ax.loglog(
            t,
            np.abs(dpD_base[sl]),
            "--",
            lw=0.8,
            color="gray",
            alpha=0.6,
            label=r"$t_D p_D'$ (base)",
        )

        # Modified
        combined_kw = {**kw, **mod_kw}
        pD_mod = _pressure(geom_inf, full_time, M, **combined_kw)
        dpD_mod = _derivative(geom_inf, full_time, M, **combined_kw)
        ax.loglog(t, pD_mod[sl], "o-", ms=2, lw=0.8, label=r"$p_D$")
        ax.loglog(t, np.abs(dpD_mod[sl]), "s-", ms=2, lw=0.8, label=r"$t_D p_D'$")

        lim = modifier_limits[i] if modifier_limits else None
        _setup_loglog(ax, mod_label, r"$t_D$", r"$p_D$", limits=lim)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Figure 10: Hybrid rate analysis (no-flow boundary, all geometries)
# ---------------------------------------------------------------------------


def _hybrid_figure(
    M: int,
    hybrid_limits: Optional[List[Tuple[float, float, Optional[float], Optional[float]]]] = None,
) -> plt.Figure:
    """Generate hybrid rate analysis figure for all geometries (no-flow)."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("Hybrid Rate Analysis — No-Flow Boundary", fontsize=12)

    # Display time arrays per geometry
    radial_display = 10 ** np.linspace(4, 6, 81)
    linear_display = LINEAR_TIME
    frac_display = FRAC_TIME

    configs = [
        ("Radial", radial_bounded_nf, radial_display, dict(reD=RADIAL_RED)),
        ("Linear", linear_bounded_nf, linear_display, dict(xD=0.0)),
        ("Fractured", frac_bounded_nf, frac_display, dict(reD=FRAC_RED)),
    ]

    for i, (ax, (geom_name, geom_fn, display_time, kw)) in enumerate(zip(axes, configs)):
        # Compute on padded time: 2 decades before, 1 decade after display range
        full_time, n_before, n_display = _padded_time(display_time)
        qD_raw = _rate(geom_fn, full_time, M, **kw)
        qD_hybrid, qDd, qDdr, qDdi, qDdid = hybrid_rate(qD_raw, full_time, b=ARPS_B)

        # Plot only the display range
        sl = slice(n_before, n_before + n_display)
        t = full_time[sl]
        ax.loglog(t, np.abs(qD_hybrid[sl]), "o-", ms=2, lw=0.8, label=r"$q_D$")
        ax.loglog(t, np.abs(qDd[sl]), "s-", ms=2, lw=0.8, label=r"$q_{Dd}$")
        ax.loglog(t, np.abs(qDdr[sl]), "d-", ms=2, lw=0.8, label=r"$q_{Ddr}$")
        ax.loglog(t, np.abs(qDdi[sl]), "^-", ms=2, lw=0.8, label=r"$q_{Ddi}$")
        ax.loglog(t, np.abs(qDdid[sl]), "v-", ms=2, lw=0.8, label=r"$q_{Ddid}$")

        lim = hybrid_limits[i] if hybrid_limits else None
        _setup_loglog(ax, geom_name, r"$t_D$", r"$q_D$", limits=lim)

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(skip_validation: bool = False) -> None:
    mp.dps = DPS
    M = 8

    OUTDIR.mkdir(parents=True, exist_ok=True)

    # --- Validation ---
    if not skip_validation:
        t0 = timer.perf_counter()
        run_all_validation()
        print(f"Validation completed in {timer.perf_counter() - t0:.1f}s\n")

    # --- Figures ---
    print("=" * 78)
    print("  Generating figures (M=8)")
    print("=" * 78)

    def _save(name: str, fig: plt.Figure) -> None:
        path = OUTDIR / f"{name}.png"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"    saved {path}")

    # 1-2: Radial
    print("\n  Radial boundary comparison...")
    t0 = timer.perf_counter()
    fig_p, fig_q = _boundary_figures(
        "Radial",
        radial_inf,
        radial_bounded_nf,
        radial_bounded_cp,
        RADIAL_TIME,
        M,
        geom_kwargs_nf=dict(reD=RADIAL_RED),
        geom_kwargs_cp=dict(reD=RADIAL_RED),
        pressure_limits=RADIAL_PRESSURE_LIMITS,
        rate_limits=RADIAL_RATE_LIMITS,
    )
    _save("radial_pressure", fig_p)
    _save("radial_rate", fig_q)
    print(f"    {timer.perf_counter() - t0:.1f}s")

    # 3-4: Linear
    print("  Linear boundary comparison...")
    t0 = timer.perf_counter()
    fig_p, fig_q = _boundary_figures(
        "Linear",
        linear_inf,
        linear_bounded_nf,
        linear_bounded_cp,
        LINEAR_TIME,
        M,
        geom_kwargs_inf=dict(xD=0.0),
        geom_kwargs_nf=dict(xD=0.0),
        geom_kwargs_cp=dict(xD=0.0),
        pressure_limits=LINEAR_PRESSURE_LIMITS,
        rate_limits=LINEAR_RATE_LIMITS,
    )
    _save("linear_pressure", fig_p)
    _save("linear_rate", fig_q)
    print(f"    {timer.perf_counter() - t0:.1f}s")

    # 5-6: Fractured
    print("  Fractured boundary comparison...")
    t0 = timer.perf_counter()
    fig_p, fig_q = _boundary_figures(
        "Fractured Well",
        frac_inf,
        frac_bounded_nf,
        frac_bounded_cp,
        FRAC_TIME,
        M,
        geom_kwargs_nf=dict(reD=FRAC_RED),
        geom_kwargs_cp=dict(reD=FRAC_RED),
        pressure_limits=FRAC_PRESSURE_LIMITS,
        rate_limits=FRAC_RATE_LIMITS,
    )
    _save("frac_pressure", fig_p)
    _save("frac_rate", fig_q)
    print(f"    {timer.perf_counter() - t0:.1f}s")

    # 7-9: Modifiers
    print("  Modifier figures...")
    t0 = timer.perf_counter()
    _save("radial_modifiers", _modifier_figure(
        "Radial", radial_inf, RADIAL_TIME, M,
        modifier_limits=RADIAL_MODIFIER_LIMITS,
    ))
    _save("linear_modifiers", _modifier_figure(
        "Linear", linear_inf, LINEAR_TIME, M,
        geom_kwargs_inf=dict(xD=0.0),
        modifier_limits=LINEAR_MODIFIER_LIMITS,
    ))
    _save("frac_modifiers", _modifier_figure(
        "Fractured Well", frac_inf, FRAC_TIME, M,
        modifier_limits=FRAC_MODIFIER_LIMITS,
    ))
    print(f"    {timer.perf_counter() - t0:.1f}s")

    # 10: Hybrid rate analysis
    print("  Hybrid rate analysis...")
    t0 = timer.perf_counter()
    _save("hybrid_rate_analysis", _hybrid_figure(M, hybrid_limits=HYBRID_LIMITS))
    print(f"    {timer.perf_counter() - t0:.1f}s")

    print("\nDone.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate diagnostic figures.")
    parser.add_argument(
        "--skip-validation", action="store_true", help="Skip M=8 vs M=16 validation (faster)"
    )
    args = parser.parse_args()
    main(skip_validation=args.skip_validation)
