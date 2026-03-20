"""
Petroleum engineering Laplace transform benchmarks and figures.

Laplace-domain functions extracted from SPE notebooks:
  - SPE-011463: Infinite-acting radial flow with variable wellbore storage
  - SPE-011964/018160: Naturally fractured reservoir (dual-porosity)
  - SPE-013169: Naturally fractured with bounded reservoir
  - ETX gas well: Fractured well with no-flow outer boundary

Usage:
    PYTHONPATH=. python test/benchmark_petroleum.py
"""

import time as timer
from math import exp, log10
from pathlib import Path
from typing import Any, Callable

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from mpmath import mp

from gwr_inversion import gwr, besselk, besseli, HAS_FLINT
from gwr_inversion._gmpy2 import HAS_GMPY2


# ---------------------------------------------------------------------------
# Laplace-domain functions (module-level for pickling in parallel mode)
# ---------------------------------------------------------------------------


def lap_rad_inf(s: Any) -> Any:
    """Infinite-acting radial flow.

    SPE-011463. The simplest reservoir model.
    """
    sqrt_s = mp.sqrt(s)
    k0 = besselk(0, sqrt_s)
    k1 = besselk(1, sqrt_s)
    return k0 / (sqrt_s * k1) / s


def lap_rad_inf_wbs(s: Any, cD: float = 0.0) -> Any:
    """Infinite-acting radial flow with wellbore storage.

    SPE-011463. The simplest reservoir model.
    """
    sqrt_s = mp.sqrt(s)
    k0 = mp.besselk(0, sqrt_s)
    k1 = mp.besselk(1, sqrt_s)
    psd = k0 / (sqrt_s * k1) / s
    return 1.0 / (1.0 / psd + s ** 2 * cD)



def lap_rad_inf_wbs_arb(s: Any, cD: float = 0.0) -> Any:
    """Same as lap_rad_inf_wbs but using ARB Bessel functions."""
    sqrt_s = mp.sqrt(s)
    k0 = besselk(0, sqrt_s)
    k1 = besselk(1, sqrt_s)
    psd = k0 / (sqrt_s * k1) / s
    return 1.0 / (1.0 / psd + s ** 2 * cD)



def lap_rad_inf_vwbs(s: Any, cD: float = 0.0, cDd: float = 1.0, aD: float = 1.0) -> Any:
    """Infinite-acting radial flow with variable wellbore storage.

    SPE-011463. Models changing wellbore storage coefficient over time
    (e.g., phase segregation, changing compressibility).
    """
    sqrt_s = mp.sqrt(s)
    k0 = besselk(0, sqrt_s)
    k1 = besselk(1, sqrt_s)
    psd = k0 / (sqrt_s * k1) / s
    wbs = 1 + s ** 2 * cD * cDd * (1 / s - 1 / (s + 1 / aD))
    return wbs / (1.0 / psd + s ** 2 * cD)



def lap_rad_inf_natf(s: Any, lam: float = 1e-6, omega: float = 0.1, cD: float = 0.0) -> Any:
    """Naturally fractured reservoir (dual-porosity), infinite-acting.

    SPE-011964, SPE-018160. Warren & Root (1963) model.
    lam = interporosity flow coefficient (matrix-fracture coupling)
    omega = storativity ratio (fracture storage / total storage)
    """
    fs = (lam + (omega * (1.0 - omega) * s)) / (lam + ((1.0 - omega) * s))
    sqrt_fs = mp.sqrt(fs * s)
    k0 = besselk(0, sqrt_fs)
    k1 = besselk(1, sqrt_fs)
    psd = k0 / (sqrt_fs * k1) / s
    return 1.0 / (1.0 / psd + s ** 2 * cD)



def lap_rad_bounded_qD(s: Any, reD: float = 1000.0, cD: float = 0.0) -> Any:
    """
    Bounded radial flow (no-flow outer boundary).
    """
    sqrt_s = mp.sqrt(s)

    k0 = besselk(0, sqrt_s)
    i0 = besseli(0, sqrt_s)
    i1_reD = besseli(1, sqrt_s * reD)
    k1_reD = besselk(1, sqrt_s * reD)
    k1 = besselk(1, sqrt_s)
    i1 = besseli(1, sqrt_s)

    psd = (
        (k0 * i1_reD + i0 * k1_reD)
        / (sqrt_s * k1 * i1_reD - sqrt_s * i1 * k1_reD)
    ) / s
    # pwd = 1.0 / (1.0 / psd + s ** 2 * cD)
    qd = 1.0 / (s ** 2 * psd)
    return qd



def lap_rad_constant_press(s: Any, reD: float = 1000.0, cD: float = 0.0) -> Any:
    """
    Radial flow, constant-pressure outer boundary.
    """
    sqrt_s = mp.sqrt(s)

    k0 = besselk(0, sqrt_s)
    i0 = besseli(0, sqrt_s)
    k0_reD = besselk(0, sqrt_s * reD)
    i0_reD = besseli(0, sqrt_s * reD)
    k1 = besselk(1, sqrt_s)
    i1 = besseli(1, sqrt_s)

    psd = (
        (k0 * i0_reD - i0 * k0_reD)
        / (sqrt_s * k1 * i0_reD + sqrt_s * i1 * k0_reD)
    ) / s
    # pwd = 1.0 / (1.0 / psd + s ** 2 * cD)
    # qd = 1.0 / (s ** 2 * psd)
    return psd



def lap_rad_nat_frac_bounded(s: Any, lam: float = 1e-6, omega: float = 0.1, reD: float = 1000.0, cD: float = 0.0) -> Any:
    """
    Radial flow, naturally fractured, constant-pressure outer boundary.
    """
    fs = (lam + (omega * (1.0 - omega) * s)) / (lam + ((1.0 - omega) * s))
    sqrt_fs = mp.sqrt(fs * s)

    k0 = besselk(0, sqrt_fs)
    i0 = besseli(0, sqrt_fs)
    i1_reD = besseli(1, sqrt_fs * reD)
    k1_reD = besselk(1, sqrt_fs * reD)
    k1 = besselk(1, sqrt_fs)
    i1 = besseli(1, sqrt_fs)

    psd = (
        (k0 * i1_reD + i0 * k1_reD)
        / (sqrt_fs * k1 * i1_reD - sqrt_fs * i1 * k1_reD)
    ) / s
    pwd = 1.0 / (1.0 / psd + s ** 2 * cD)
    return psd



def lap_linear_nf(s: Any, xD: float = 0.0) -> Any:
    """Linear flow with no-flow outer boundary.

    pwd = 1/(s * sqrt(s)) * cosh(sqrt(s)*(1-xD)) / sinh(sqrt(s))
    xD = dimensionless distance from wellbore (0 = at well, 1 = at boundary)
    """
    sqrt_s = mp.sqrt(s)
    return mp.cosh(sqrt_s * (1 - xD)) / (s * sqrt_s * mp.sinh(sqrt_s))


def _lap_frac_rad_inf(sqrt_s: Any, xD: float) -> Any:
    """Infinite-acting component of fractured well solution."""
    k0_int_1 = mp.quad(lambda x: besselk(0, x), [0.0, sqrt_s * (1 - xD)])
    k0_int_2 = mp.quad(lambda x: besselk(0, x), [0.0, sqrt_s * (1 + xD)])
    return k0_int_1 + k0_int_2



def _lap_frac_rad_nf(sqrt_s: Any, xD: float) -> Any:
    """No-flow boundary component of fractured well solution."""
    i0_int_1 = mp.besseli(0, sqrt_s * (1 - xD), -1)
    i0_int_2 = mp.besseli(0, sqrt_s * (1 + xD), -1)
    return i0_int_1 + i0_int_2



def lap_frac_rad(s: Any, reD: float = 100.0, xD: float = 0.732) -> Any:
    """Fractured well with no-flow outer boundary.

    ETX gas well RTA. Vertically fractured well model.
    reD = dimensionless external radius
    xD = dimensionless fracture half-length position (default 0.732 for uniform flux)

    Uses mp.quad (numerical quadrature) inside the Laplace function,
    making this the most computationally expensive transform.
    """
    sqrt_s = mp.sqrt(s)
    c_ = 1 / (2 * s * sqrt_s)

    inf = c_ * _lap_frac_rad_inf(sqrt_s, xD)

    k1 = besselk(1, sqrt_s * reD)
    i1 = besseli(1, sqrt_s * reD)
    nf = c_ * k1 / i1 * _lap_frac_rad_nf(sqrt_s, xD)

    return inf + nf


# ---------------------------------------------------------------------------
# Benchmark helpers
# ---------------------------------------------------------------------------

def bench_row(method: str, elapsed: float, n_points: int) -> None:
    print(f"  {method:<45} {elapsed:>8.3f}s  ({elapsed/n_points*1000:>6.1f} ms/pt)")


def run_benchmark(name: str, fn: Callable[..., Any], time_array: "np.ndarray[Any, Any]", M: int = 8) -> "np.ndarray[Any, Any]":
    """Benchmark a single Laplace function across backends."""
    n = len(time_array)
    print(f"\n{name} ({n} points, M={M}):")
    print(f"  {'Method':<45} {'Total':>8}   {'Per point':>10}")
    print(f"  {'-'*68}")

    # mpmath
    t0 = timer.perf_counter()
    ref = np.array([float(v) for v in gwr(fn, time_array, M, backend="mpmath")])
    bench_row("mpmath backend", timer.perf_counter() - t0, n)

    if HAS_GMPY2:
        t0 = timer.perf_counter()
        r = np.array([float(v) for v in gwr(fn, time_array, M, backend="gmpy2")])
        bench_row("gmpy2 backend", timer.perf_counter() - t0, n)

    return ref


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_pressure_derivative(
    ax: Axes, time: "np.ndarray[Any, Any]",
    title: str,
    series: list["np.ndarray[Any, Any]"],
    labels: list[str] = ["p_D", "p_Dd"],
    xlims: tuple[float | None, float | None] = (None, None),
    ylims: tuple[float | None, float | None] = (None, None),
) -> None:
    """Standard log-log pressure/derivative diagnostic plot."""
    ax.loglog(time, series[0], 'o-', ms=3, lw=0.8, label=labels[0])
    ax.loglog(time, np.abs(series[1]), 's-', ms=3, lw=0.8, label=labels[1])
    if len(series) > 2:
        ax.loglog(time, np.abs(series[2]), 'd-', ms=3, lw=0.8, label=labels[2])
    if len(series) > 3:
        ax.loglog(time, np.abs(series[3]), 'x-', ms=3, lw=0.8, label=labels[3])
    if len(series) > 4:
        ax.loglog(time, np.abs(series[4]), 'v-', ms=3, lw=0.8, label=labels[4])

    ax.set_xlabel("Dimensionless time, $t_D / C_D$")
    ax.set_ylabel("Dimensionless pressure / derivative")
    ax.set_title(title)
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.legend(fontsize=8)
    ax.set_aspect(1, 'datalim')
    ax.grid(True, which='major', ls=':', alpha=0.75)


def generate_figures(outdir: str) -> None:
    """Generate diagnostic plots for all petroleum Laplace transforms."""
    outpath = Path(outdir)
    outpath.mkdir(parents=True, exist_ok=True)

    mp.dps = 50
    M = 8

    fig, axes = plt.subplots(3, 3, figsize=(18, 18))

    # 1. Infinite-acting radial flow (cD = 1e3)
    time = 10 ** np.linspace(-5, 5, 61)

    pD = np.array([float(v) for v in gwr(
        lambda s: lap_rad_inf(s), time, M)])
    dpD = time * np.array([float(v) for v in gwr(
        lambda s: s * lap_rad_inf(s), time, M)])

    plot_pressure_derivative(
        axes[0, 0], time,
        "Infinite-acting Radial Flow",
        [pD, dpD]
    )

    # 2. Infinite-acting radial flow w/ wellbore storage (cD = 1e3)
    cD = 1e13
    _cD = round(log10(cD), 0)
    time = 10 ** np.linspace(_cD - 5, _cD + 5, 61)

    pD = np.array([float(v) for v in gwr(
        lambda s: lap_rad_inf_wbs(s, cD=cD), time, M)])
    dpD = time * np.array([float(v) for v in gwr(
        lambda s: s * lap_rad_inf_wbs(s, cD=cD), time, M)])

    plot_pressure_derivative(
        axes[0, 1], time / cD,
        "Infinite-acting Radial Flow + WBS (mpmath)",
        [pD, dpD]
    )

    # 3. Infinite-acting radial flow w/ wellbore storage (arb) (cD = 1e3)
    cD = 1e13
    _cD = round(log10(cD), 0)
    time = 10 ** np.linspace(_cD - 5, _cD + 5, 61)

    pD = np.array([float(v) for v in gwr(
        lambda s: lap_rad_inf_wbs_arb(s, cD=cD), time, M)])
    dpD = time * np.array([float(v) for v in gwr(
        lambda s: s * lap_rad_inf_wbs_arb(s, cD=cD), time, M)])

    plot_pressure_derivative(
        axes[0, 2], time / cD,
        "Infinite-acting Radial Flow + WBS (ARB)",
        [pD, dpD]
    )

    # 4. Infinite-acting radial flow w/ Variable wellbore storage
    cD, cDd, aD = 1e7, 1.0, -.575
    _cD = round(log10(cD), 0)
    time = 10 ** np.linspace(_cD - 5, _cD + 5, 61)

    pD = np.array([float(v) for v in gwr(
        lambda s: lap_rad_inf_vwbs(s, cD=cD, cDd=cDd, aD=aD), time, M)])
    dpD = time * np.array([float(v) for v in gwr(
        lambda s: s * lap_rad_inf_vwbs(s, cD=cD, cDd=cDd, aD=aD), time, M)])

    plot_pressure_derivative(
        axes[1, 0], time / cD,
        "Infinite-acting Radial Flow w/ Variable WBS\n",
        [pD, dpD]
    )

    # 5. Naturally fractured, infinite-acting
    lam, omega, cD = 1e-6, 0.1, 100.0
    time = 10 ** np.linspace(-1, 7, 61)

    pD = np.array([float(v) for v in gwr(
        lambda s: lap_rad_inf_natf(s, lam=lam, omega=omega, cD=cD), time, M)])
    dpD = time * np.array([float(v) for v in gwr(
        lambda s: s * lap_rad_inf_natf(s, lam=lam, omega=omega, cD=cD), time, M)])

    plot_pressure_derivative(
        axes[1, 1], time / cD,
        "Infinite-acting Radial Flow, Naturally Fractured\n",
        [pD, dpD]
    )

    # 6. Bounded radial flow (homogeneous, no-flow)
    from scipy.integrate import simpson
    M = 8
    skin = -4.75
    rw = 0.354
    reD = 1000.0 / (rw * exp(-skin))
    cD = .0001 * exp(2 * -skin)
    reD = round(reD, 0)
    cD = exp(round(log10(cD), 0))
    arps_b = 0.5

    time = 10 ** np.linspace(-2, 6, 101)

    qD = np.array([float(v) for v in gwr(
        lambda s: lap_rad_bounded_qD(s, reD=reD, cD=cD), time, M)])
    qD_hyp = 1 / (1 + time * arps_b) ** (1 / arps_b) - np.exp(-time)
    qD_hybrid = qD + qD_hyp

    qDd = time * np.abs(np.gradient(qD_hybrid, time))

    qDdr = 1 / (time * np.gradient(1 / qD_hybrid, time))

    qDdi = np.array([simpson(qD_hybrid[:i + 1], time[:i + 1]) for i in range(len(time))]) / time
    B = mp.log(qD_hybrid[1] / qD_hybrid[0]) / mp.log(time[1] / time[0])
    A = qD_hybrid[0] / time[0] ** B
    qDdi[0] = A / (B + 1) * time[0] ** (B + 1)

    qDdid = qDdi - qD_hybrid

    plot_pressure_derivative(
        axes[1, 2], time,
        "Radial flow Hybrid Model (No-Flow Boundary)\n",
        [qD, qDd, qDdr, qDdi, qDdid],
        labels=["q_D", "q_Dd", "q_Ddr", "q_Ddi", "q_Ddid"],
        xlims=(1e1, 1e5),
        ylims=(1e-3, 1e1)
    )

    # 7. Radial Flow, constant-pressure boundary
    M = 16
    cD - 1e13

    time = 10 ** np.linspace(-3, 3, 101)

    pD = np.array([float(v) for v in gwr(
        lambda s: lap_rad_constant_press(s, reD=reD, cD=cD), time, M)])
    pDd = time * np.array([float(v) for v in gwr(
        lambda s: s * lap_rad_constant_press(s, reD=reD, cD=cD), time, M)])

    plot_pressure_derivative(
        axes[2, 0], time / cD,
        "Radial flow (Const-P Boundary)\n",
        [pD, pDd]
    )

    # 8. Radial Flow, naturally fractured, constant-pressure boundary
    M = 8
    lam, omega = 1e-6, 0.1
    pD = np.array([float(v) for v in gwr(
        lambda s: lap_rad_nat_frac_bounded(s, lam=lam, omega=omega, reD=reD, cD=cD),
        time, M)])
    dpD = time * np.array([float(v) for v in gwr(
        lambda s: s * lap_rad_nat_frac_bounded(s, lam=lam, omega=omega, reD=reD, cD=cD),
        time, M)])

    plot_pressure_derivative(
        axes[2, 1], time / cD,
        "Radial flow, naturally fractured (No-Flow Boundary)\n",
        [np.abs(pD), dpD], labels=["q_D", "dq_D"]
    )

    axes[2, 2].set_visible(False)

    plt.tight_layout()
    fig.savefig(outpath / "benchmarks.png", dpi=150)
    print(f"\nSaved figure to {outpath / 'benchmarks.png'}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    mp.dps = 50

    print("=" * 72)
    print("  Laplace Transform Benchmarks")
    print("=" * 72)
    print(f"\n  gmpy2: {'installed' if HAS_GMPY2 else 'not installed'}")
    print(f"  flint: {'installed' if HAS_FLINT else 'not installed'}")

    # 1. Infinite-acting radial flow (simple — besselk only)
    cD = 1e3
    time = 10 ** np.linspace(1, 7, 31)
    run_benchmark(
        "Infinite-acting Radial Flow + WBS (ARB)",
        lambda s: lap_rad_inf_wbs_arb(s, cD=cD), time, M=8)

    # 2. Variable wellbore storage
    cD, cDd, aD = 1e3, 10.0, 0.01
    run_benchmark(
        "Infinite-acting Radial Flow w/ Variable WBS",
        lambda s: lap_rad_inf_vwbs(s, cD=cD, cDd=cDd, aD=aD), time, M=8)

    # 3. Naturally fractured, infinite-acting
    lam, omega, cD = 1e-6, 0.1, 100.0
    run_benchmark(
        "Infinite-acting Radial Flow, Naturally Fractured",
        lambda s: lap_rad_inf_natf(s, lam=lam, omega=omega, cD=cD), time, M=8)

    # 4. Radial flow, no-flow boundary (homogeneous)
    reD, cD = 1000.0, 100.0
    run_benchmark(
        "Radial flow, no-flow boundary (reD=1000)",
        lambda s: lap_rad_bounded_qD(s, reD=reD, cD=cD), time, M=8)

    # 5. Radial flow, constant-pressure boundary
    run_benchmark(
        "Radial flow, const-P boundary (reD=1000)",
        lambda s: lap_rad_constant_press(s, reD=reD, cD=cD),
        time, M=8)

    # 6. Radial flow, naturally fractured, no-flow boundary
    run_benchmark(
        "Radial flow, naturally fractured, no-flow boundary\n"
        r"$\lambda = 10^{-6}, \omega = 0.1, r_{eD} = 1000$",
        lambda s: lap_rad_nat_frac_bounded(s, lam=lam, omega=omega, reD=reD, cD=cD),
        time, M=8)

    # 7. Fractured well (expensive — uses quadrature)
    reD = 100.0
    time = 10 ** np.linspace(-2, 2, 21)
    run_benchmark(
        "Fractured well, no-flow boundary (reD=100)",
        lambda s: lap_frac_rad(s, reD=reD), time, M=8)

    # Generate figures
    print("\n" + "=" * 72)
    print("  Generating figures...")
    print("=" * 72)
    generate_figures("test")


if __name__ == '__main__':
    main()
