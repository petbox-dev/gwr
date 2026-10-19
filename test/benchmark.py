"""
Petroleum engineering Laplace transform benchmarks and figures.

Times the library transforms in ``gwr_inversion.transforms`` under three
configurations and draws one diagnostic panel per case:

  A  backend="mpmath", Bessel functions from mpmath   (pure-mpmath baseline)
  B  backend="gmpy2",  Bessel functions from mpmath   (gmpy2 GWR core only)
  C  backend="gmpy2",  library Bessel dispatch         (default install, [all])

"Bessel functions from mpmath" runs the library as if python-flint and the
gmpy2 Bessel code were not installed: K_n / I_n, the outer-boundary K/I
ratio and the fractured-well K_0 / I_0 integrals all take their mpmath
fallbacks, and lap_transform skips its ARB ball-arithmetic path.

Usage:
    PYTHONPATH=. python test/benchmark.py              # M = 8 and 32
    PYTHONPATH=. python test/benchmark.py --orders 8   # M = 8 only
"""

import argparse
import time as timer
from pathlib import Path
from typing import Any, Callable, Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from mpmath import mp

import gwr_inversion._bessel as _bessel
import gwr_inversion._flint as _flint
import gwr_inversion.transforms as transforms
from gwr_inversion import (
    HAS_FLINT, cache_clear, frac_bounded_nf, gwr, hybrid_rate, lap_transform,
    linear_bounded_nf, output_derivative, output_rate, radial_bounded_cp, radial_bounded_nf,
    radial_inf,
)
from gwr_inversion._gmpy2 import HAS_GMPY2


LaplaceFunction = Callable[[Any], Any]


# ---------------------------------------------------------------------------
# Bessel configurations
# ---------------------------------------------------------------------------


LIBRARY_ARB_GEOMETRIES = dict(transforms._ARB_GEOMETRIES)
LIBRARY_FLAGS = (_flint.HAS_FLINT, _bessel.HAS_FLINT, _bessel.HAS_GMPY2_BESSEL)


def _use_bessel(library: bool) -> None:
    """Library dispatch (True) or the mpmath fallbacks everywhere (False)."""
    flint_flag, dispatch_flint_flag, gmpy2_bessel_flag = (
        LIBRARY_FLAGS if library else (False, False, False))
    _flint.HAS_FLINT = flint_flag
    _bessel.HAS_FLINT = dispatch_flint_flag
    _bessel.HAS_GMPY2_BESSEL = gmpy2_bessel_flag
    transforms._ARB_GEOMETRIES = LIBRARY_ARB_GEOMETRIES if library else {}


CONFIGURATIONS: List[Tuple[str, str, bool]] = [
    ("A mpmath core + mpmath Bessel", "mpmath", False),
    ("B gmpy2 core + mpmath Bessel", "gmpy2", False),
    ("C gmpy2 core + library Bessel", "gmpy2", True),
]


# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------

LAMBDA, OMEGA = 1e-6, 0.1
RADIAL_TIME = 10 ** np.linspace(1, 7, 31)
FRAC_TIME = 10 ** np.linspace(-2, 2, 21)

TIMING_CASES: List[Tuple[str, LaplaceFunction, "np.ndarray[Any, Any]"]] = [
    ("Radial inf + WBS (cD=1e3), pressure",
     lambda s: lap_transform(s, radial_inf, cD=1e3), RADIAL_TIME),
    ("Radial inf + variable WBS (1e3, 10, 0.01), pressure",
     lambda s: lap_transform(s, radial_inf, cD=1e3, cDd=10.0, aD=0.01), RADIAL_TIME),
    ("Radial inf + dual porosity + WBS (cD=100), pressure",
     lambda s: lap_transform(s, radial_inf, lam=LAMBDA, omega=OMEGA, cD=100.0), RADIAL_TIME),
    ("Radial no-flow (reD=1000), rate",
     lambda s: output_rate(lap_transform(s, radial_bounded_nf, reD=1000.0), s), RADIAL_TIME),
    ("Radial const-P (reD=1000), pressure",
     lambda s: lap_transform(s, radial_bounded_cp, reD=1000.0), RADIAL_TIME),
    ("Radial no-flow (reD=1000) + dual porosity, pressure",
     lambda s: lap_transform(s, radial_bounded_nf, reD=1000.0, lam=LAMBDA, omega=OMEGA),
     RADIAL_TIME),
    ("Linear no-flow, pressure (no Bessel)",
     lambda s: lap_transform(s, linear_bounded_nf, xD=0.0), RADIAL_TIME),
    ("Fractured well no-flow (reD=100), pressure",
     lambda s: lap_transform(s, frac_bounded_nf, reD=100.0), FRAC_TIME),
]


def run_timings(orders: Sequence[int]) -> None:
    """Time every case under every configuration; report max rel diff vs A."""
    for M in orders:
        print(f"\nM = {M}")
        print(f"  {'Case':<54} {'A (s)':>8} {'B (s)':>8} {'C (s)':>8} {'A/C':>6} "
              f"{'C ms/pt':>8} {'max rel C-A':>12}")
        print(f"  {'-' * 110}")
        for label, laplace_function, time_array in TIMING_CASES:
            results: Dict[str, "np.ndarray[Any, Any]"] = {}
            seconds: Dict[str, float] = {}
            for name, backend, library in CONFIGURATIONS:
                if backend == "gmpy2" and not HAS_GMPY2:
                    continue
                _use_bessel(library)
                cache_clear()
                start = timer.perf_counter()
                results[name[0]] = np.asarray(gwr(laplace_function, time_array, M,
                                                  backend=backend))
                seconds[name[0]] = timer.perf_counter() - start
            _use_bessel(True)
            if "C" not in seconds:
                print(f"  {label:<54} {seconds['A']:8.2f}  (gmpy2 not installed)")
                continue
            relative = np.max(np.abs((results["C"] - results["A"]) / results["A"]))
            print(f"  {label:<54} {seconds['A']:8.2f} {seconds['B']:8.2f} {seconds['C']:8.2f} "
                  f"{seconds['A'] / seconds['C']:6.1f} "
                  f"{seconds['C'] / len(time_array) * 1e3:8.1f} {relative:12.2e}", flush=True)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def _square_decade_limits(
    axis: Axes,
    x_limits: Tuple[float, float],
    series: Sequence["np.ndarray[Any, Any]"],
    x_values: "np.ndarray[Any, Any]",
) -> None:
    """Fix x to ``x_limits`` and y to the same number of decades, ending one
    half-decade above the largest plotted value inside ``x_limits``."""
    inside = (x_values >= x_limits[0]) & (x_values <= x_limits[1])
    largest = max(float(np.nanmax(np.abs(values[inside]))) for values in series)
    y_high = 10 ** (np.ceil(np.log10(largest) * 2) / 2 + 0.5)
    decades = np.log10(x_limits[1]) - np.log10(x_limits[0])
    axis.set_xlim(*x_limits)
    axis.set_ylim(y_high / 10 ** decades, y_high)
    axis.set_aspect("equal", adjustable="box")


def _plot_panel(
    axis: Axes,
    title: str,
    x_values: "np.ndarray[Any, Any]",
    x_label: str,
    series: Sequence["np.ndarray[Any, Any]"],
    labels: Sequence[str],
    x_limits: Tuple[float, float],
) -> None:
    markers = ["o-", "s-", "d-", "^-", "v-"]
    for values, label, marker in zip(series, labels, markers):
        axis.loglog(x_values, np.abs(values), marker, ms=2, lw=0.8, label=label)
    axis.set_title(title, fontsize=10)
    axis.set_xlabel(x_label, fontsize=9)
    axis.set_ylabel("Dimensionless pressure / rate", fontsize=9)
    axis.legend(fontsize=8)
    axis.grid(True, which="major", ls=":", alpha=0.6)
    _square_decade_limits(axis, x_limits, series, x_values)


def _pressure_and_derivative(
    geometry: Callable[..., Any],
    time_array: "np.ndarray[Any, Any]",
    M: int,
    **kwargs: Any,
) -> Tuple["np.ndarray[Any, Any]", "np.ndarray[Any, Any]"]:
    pressure = gwr(lambda s: lap_transform(s, geometry, **kwargs), time_array, M)
    derivative = time_array * gwr(
        lambda s: output_derivative(lap_transform(s, geometry, **kwargs), s), time_array, M)
    return pressure, derivative


def generate_figures(outdir: str) -> None:
    """Draw one log-log diagnostic panel per benchmark case."""
    outpath = Path(outdir)
    outpath.mkdir(parents=True, exist_ok=True)
    M = 8
    pressure_labels = [r"$p_D$", r"$t_D \, p_D'$"]
    figure, axes = plt.subplots(3, 3, figsize=(18, 18))

    # 1. Infinite-acting radial flow
    time_array = 10 ** np.linspace(-3, 5, 81)
    _plot_panel(axes[0, 0], "Infinite-acting Radial Flow", time_array, r"$t_D$",
                _pressure_and_derivative(radial_inf, time_array, M), pressure_labels,
                (1e-3, 1e5))

    # 2. Infinite-acting radial flow + WBS
    wellbore_storage = 1e13
    time_array = 10 ** np.linspace(8, 18, 101)
    _plot_panel(axes[0, 1], r"Radial Flow + WBS ($C_D = 10^{13}$)",
                time_array / wellbore_storage, r"$t_D / C_D$",
                _pressure_and_derivative(radial_inf, time_array, M, cD=wellbore_storage),
                pressure_labels, (1e-5, 1e5))

    # 3. Fractured well, no-flow boundary
    time_array = 10 ** np.linspace(-3, 5, 81)
    _plot_panel(axes[0, 2], r"Fractured Well, No-Flow ($r_{eD} = 100$)", time_array,
                r"$t_{Dxf}$",
                _pressure_and_derivative(frac_bounded_nf, time_array, M, reD=100.0),
                pressure_labels, (1e-3, 1e5))

    # 4. Variable wellbore storage
    wellbore_storage = 1e7
    time_array = 10 ** np.linspace(2, 12, 101)
    _plot_panel(axes[1, 0], "Radial Flow, Variable WBS", time_array / wellbore_storage,
                r"$t_D / C_D$",
                _pressure_and_derivative(radial_inf, time_array, M, cD=wellbore_storage,
                                         cDd=1.0, aD=-0.575),
                pressure_labels, (1e-5, 1e5))

    # 5. Naturally fractured, infinite-acting, WBS
    wellbore_storage = 100.0
    time_array = 10 ** np.linspace(-1, 7, 81)
    _plot_panel(axes[1, 1], r"Radial Flow, Dual Porosity ($\lambda=10^{-6}, \omega=0.1$)",
                time_array / wellbore_storage, r"$t_D / C_D$",
                _pressure_and_derivative(radial_inf, time_array, M, lam=LAMBDA, omega=OMEGA,
                                         cD=wellbore_storage),
                pressure_labels, (1e-3, 1e5))

    # 6. Hybrid rate analysis, bounded radial flow
    outer_radius = 24.0
    time_array = 10 ** np.linspace(-1, 6, 141)
    rate = gwr(lambda s: output_rate(lap_transform(s, radial_bounded_nf, reD=outer_radius), s),
               time_array, M)
    hybrid_series = hybrid_rate(rate, time_array, b=0.5)
    display = time_array >= 1e1
    _plot_panel(axes[1, 2], r"Hybrid Rate, No-Flow ($r_{eD} = 24$)", time_array[display],
                r"$t_D$", [values[display] for values in hybrid_series],
                [r"$q_D$", r"$q_{Dd}$", r"$q_{Ddr}$", r"$q_{Ddi}$", r"$q_{Ddid}$"],
                (1e1, 1e5))

    # 7. Constant-pressure boundary
    time_array = 10 ** np.linspace(-3, 3, 61)
    _plot_panel(axes[2, 0], r"Radial Flow, Const-P ($r_{eD} = 24$)", time_array, r"$t_D$",
                _pressure_and_derivative(radial_bounded_cp, time_array, 16,
                                         reD=outer_radius),
                pressure_labels, (1e-3, 1e3))

    # 8. Naturally fractured, no-flow boundary
    _plot_panel(axes[2, 1], r"Dual Porosity, No-Flow ($r_{eD} = 24$)", time_array, r"$t_D$",
                _pressure_and_derivative(radial_bounded_nf, time_array, M, reD=outer_radius,
                                         lam=LAMBDA, omega=OMEGA),
                pressure_labels, (1e-3, 1e3))

    # 9. Linear flow, no-flow boundary (no Bessel functions)
    time_array = 10 ** np.linspace(-3, 3, 61)
    _plot_panel(axes[2, 2], "Linear Flow, No-Flow", time_array, r"$t_D$",
                _pressure_and_derivative(linear_bounded_nf, time_array, M, xD=0.0),
                pressure_labels, (1e-3, 1e3))

    figure.tight_layout()
    figure.savefig(outpath / "benchmarks.png", dpi=150)
    print(f"\nSaved figure to {outpath / 'benchmarks.png'}")
    plt.close(figure)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark the library Laplace transforms.")
    parser.add_argument("--orders", type=int, nargs="+", default=[8, 32],
                        help="GWR M values to time (default: 8 32)")
    parser.add_argument("--skip-figures", action="store_true", help="Timings only")
    arguments = parser.parse_args()

    print("=" * 72)
    print("  Laplace Transform Benchmarks (gwr_inversion.transforms)")
    print("=" * 72)
    print(f"\n  gmpy2: {'installed' if HAS_GMPY2 else 'not installed'}")
    print(f"  flint: {'installed' if HAS_FLINT else 'not installed'}")
    run_timings(arguments.orders)

    if not arguments.skip_figures:
        generate_figures("test")


if __name__ == "__main__":
    main()
