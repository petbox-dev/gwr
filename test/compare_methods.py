"""
Compare numerical Laplace inversion methods and backends.

Tests four problems from the test suite:
  1. Unit step: H(t-1) = L^-1{ e^(-s)/s }
  2. Ln(t) = L^-1{ -(ln(s) + 0.577216)/s }
  3. Radial flow with wellbore storage (besselk)
  4. Hard oscillatory: sin(t) = L^-1{ 1/(1+s^2) } at t=1000

Usage:
    PYTHONPATH=. python test/compare_methods.py
"""

import time as timer
import numpy as np
from mpmath import mp, invertlaplace

from gwr_inversion import gwr, talbot, besselk
from gwr_inversion._gmpy2 import HAS_GMPY2
from gwr_inversion._flint import HAS_FLINT


# ---------------------------------------------------------------------------
# Laplace-domain functions (module-level for pickling in parallel mode)
# ---------------------------------------------------------------------------

def lap_unit(s):
    """Unit step H(t-1): e^(-s)/s."""
    return mp.exp(-s) / s


def lap_ln(s):
    """Ln(t): -(ln(s) + 0.577216)/s."""
    return -mp.log(s) / s - 0.577216 / s


def lap_rad(s):
    """Radial flow with wellbore storage (cD=1e3, no skin), mpmath besselk."""
    b = mp.besselk(0, mp.sqrt(s))
    return 1.0 / (s / b + (s ** 2) * 1e3)


def lap_rad_arb(s):
    """Radial flow with wellbore storage (cD=1e3, no skin), ARB besselk."""
    b = besselk(0, mp.sqrt(s))
    return 1.0 / (s / b + (s ** 2) * 1e3)


def lap_sin(s, *args):
    """sin(t): 1/(1+s^2)."""
    return 1.0 / (1.0 + s * s)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def section(title):
    print(f"\n{'=' * 90}")
    print(f"  {title}")
    print(f"{'=' * 90}\n")


def row(method, param, result_f, error, elapsed):
    print(f"{method:<30} {param:>10} {result_f:>25.15e} {error:>15.3e} {elapsed:>10.3f}")


def header():
    print(f"{'Method':<30} {'M/degree':>10} {'Result':>25} "
          f"{'Error':>15} {'Time (s)':>10}")
    print("-" * 95)


def array_header():
    print(f"{'Method':<45} {'Max Error':>15} {'Time (s)':>10}")
    print("-" * 75)


def array_row(method, max_err, elapsed):
    print(f"{method:<45} {max_err:>15.3e} {elapsed:>10.3f}")


def bench_gwr_scalar(fn, t_val, expected, m_values, label):
    """Benchmark GWR on a single time point at various M values (both backends)."""
    for M in m_values:
        t0 = timer.perf_counter()
        result = float(gwr(fn, t_val, M, backend="mpmath"))
        elapsed = timer.perf_counter() - t0
        row(f"GWR (mpmath)", M, result, abs(result - expected), elapsed)

    if HAS_GMPY2:
        print()
        for M in m_values:
            t0 = timer.perf_counter()
            result = float(gwr(fn, t_val, M, backend="gmpy2"))
            elapsed = timer.perf_counter() - t0
            row(f"GWR (gmpy2)", M, result, abs(result - expected), elapsed)


def bench_mpmath_scalar(fn, t_val, expected, degrees):
    """Benchmark mpmath.invertlaplace methods on a single time point."""
    for method_name in ['talbot', 'stehfest', 'dehoog', 'cohen']:
        for degree in degrees:
            try:
                t0 = timer.perf_counter()
                result = float(invertlaplace(fn, t_val, method=method_name, degree=degree))
                elapsed = timer.perf_counter() - t0
                row(f"mpmath {method_name}", degree, result,
                    abs(result - expected), elapsed)
            except Exception:
                print(f"mpmath {method_name:<24} {degree:>10} {'FAILED':>25} "
                      f"{'':>15} {'':>10}")


def bench_array(fn, time_array, expected, label):
    """Benchmark GWR and Talbot on a time array."""
    # GWR mpmath sequential
    t0 = timer.perf_counter()
    result = np.array([float(v) for v in gwr(fn, time_array, 32, backend="mpmath", workers=1)])
    elapsed = timer.perf_counter() - t0
    array_row(f"GWR M=32 mpmath sequential", np.max(np.abs(result - expected)), elapsed)

    # GWR mpmath parallel
    t0 = timer.perf_counter()
    result = np.array([float(v) for v in gwr(fn, time_array, 32, backend="mpmath", workers=4)])
    elapsed = timer.perf_counter() - t0
    array_row(f"GWR M=32 mpmath parallel (4 workers)", np.max(np.abs(result - expected)), elapsed)

    if HAS_GMPY2:
        # GWR gmpy2 sequential
        t0 = timer.perf_counter()
        result = np.array([float(v) for v in gwr(fn, time_array, 32, backend="gmpy2", workers=1)])
        elapsed = timer.perf_counter() - t0
        array_row(f"GWR M=32 gmpy2 sequential", np.max(np.abs(result - expected)), elapsed)

        # GWR gmpy2 parallel
        t0 = timer.perf_counter()
        result = np.array([float(v) for v in gwr(fn, time_array, 32, backend="gmpy2", workers=4)])
        elapsed = timer.perf_counter() - t0
        array_row(f"GWR M=32 gmpy2 parallel (4 workers)", np.max(np.abs(result - expected)), elapsed)

    # Talbot
    t0 = timer.perf_counter()
    result = np.array([float(v) for v in talbot(fn, time_array, degree=32)])
    elapsed = timer.perf_counter() - t0
    array_row(f"Talbot degree=32", np.max(np.abs(result - expected)), elapsed)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    mp.dps = 50

    # ========================================================================
    #  1. Unit step H(t-1)
    # ========================================================================
    section("Unit Step: H(t-1) evaluated at t=2.0")
    expected = 1.0
    print(f"Expected: {expected}\n")
    header()

    bench_gwr_scalar(lap_unit, 2.0, expected, [8, 16, 32, 64, 128], "(mpmath)")
    print()
    bench_mpmath_scalar(lap_unit, 2.0, expected, [8, 16, 32, 64])

    section("Unit Step: H(t-1) across 101-point time array")
    time_unit = 10 ** np.linspace(-1, 1, 101)
    expected_unit = np.where(time_unit >= 1.0, 1.0, 0.0)
    array_header()
    bench_array(lap_unit, time_unit, expected_unit, "unit step")

    # ========================================================================
    #  2. Ln(t)
    # ========================================================================
    section("Ln(t) across 31-point time array")
    time_ln = 10 ** np.linspace(-2, 5, 31)
    expected_ln = np.log(time_ln)
    print(f"Time range: {time_ln[0]:.2e} to {time_ln[-1]:.2e}, "
          f"{len(time_ln)} points\n")
    array_header()
    bench_array(lap_ln, time_ln, expected_ln, "ln(t)")

    # ========================================================================
    #  3. Radial flow with wellbore storage
    # ========================================================================
    section("Radial Flow (cD=1e3) across 31-point time array")
    time_rad = 10 ** np.linspace(1, 7, 31)
    print(f"Time range: {time_rad[0]:.2e} to {time_rad[-1]:.2e}, "
          f"{len(time_rad)} points\n")
    print("No closed-form solution — comparing GWR vs Talbot\n")

    # GWR M=32 mpmath (reference)
    t0 = timer.perf_counter()
    ref_gwr = np.array([float(v) for v in gwr(lap_rad, time_rad, 32, backend="mpmath")])
    elapsed_gwr = timer.perf_counter() - t0

    print(f"{'Method':<40} {'Time (s)':>10}   {'Max Rel Diff vs GWR M=32':>25}")
    print("-" * 80)
    print(f"{'GWR M=32 mpmath':<40} {elapsed_gwr:>10.3f}   {'(reference)':>25}")

    if HAS_FLINT:
        t0 = timer.perf_counter()
        res = np.array([float(v) for v in gwr(lap_rad_arb, time_rad, 32, backend="mpmath")])
        elapsed = timer.perf_counter() - t0
        rel = np.max(np.abs((res - ref_gwr) / ref_gwr))
        print(f"{'GWR M=32 ARB besselk + mpmath':<40} {elapsed:>10.3f}   {rel:>25.3e}")

    if HAS_GMPY2:
        for M in [8, 16, 32]:
            t0 = timer.perf_counter()
            res = np.array([float(v) for v in gwr(lap_rad, time_rad, M, backend="gmpy2")])
            elapsed = timer.perf_counter() - t0
            rel = np.max(np.abs((res - ref_gwr) / ref_gwr))
            print(f"{'GWR M=' + str(M) + ' gmpy2':<40} {elapsed:>10.3f}   {rel:>25.3e}")

    if HAS_FLINT and HAS_GMPY2:
        t0 = timer.perf_counter()
        res = np.array([float(v) for v in gwr(lap_rad_arb, time_rad, 32, backend="gmpy2")])
        elapsed = timer.perf_counter() - t0
        rel = np.max(np.abs((res - ref_gwr) / ref_gwr))
        print(f"{'GWR M=32 ARB besselk + gmpy2':<40} {elapsed:>10.3f}   {rel:>25.3e}")

    for deg in [8, 16, 32, 64]:
        t0 = timer.perf_counter()
        res = np.array([float(v) for v in talbot(lap_rad, time_rad, degree=deg)])
        elapsed = timer.perf_counter() - t0
        rel = np.max(np.abs((res - ref_gwr) / ref_gwr))
        print(f"{'Talbot degree=' + str(deg):<40} {elapsed:>10.3f}   {rel:>25.3e}")

    # ========================================================================
    #  4. Hard oscillatory: sin(t) at t=1000
    # ========================================================================
    section("Hard oscillatory: sin(t) at t=1000")
    expected_sin = float(mp.sin(1000))
    print(f"Expected: {expected_sin}\n")
    header()

    bench_gwr_scalar(lap_sin, 1000.0, expected_sin,
                     [32, 128, 512], "")

    if not HAS_GMPY2:
        print("\n  [gmpy2 not installed -- skipping gmpy2 backend]")

    print()
    bench_mpmath_scalar(lap_sin, 1000.0, expected_sin, [32, 128, 256])


if __name__ == '__main__':
    main()
