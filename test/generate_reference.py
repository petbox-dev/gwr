"""
Run once to generate reference data for regression tests.

    python test/generate_reference.py

Writes test/reference_data.npz.
"""

from functools import lru_cache

import numpy as np
from mpmath import mp

from gwr_inversion import gwr


# ---------------------------------------------------------------------------
# Laplace-domain functions
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1024)
def lap_unit(s, a):
    return 1 / s * mp.exp(-a * s)


@lru_cache(maxsize=1024)
def lap_ln(s):
    return -mp.log(s) / s - 0.577216 / s


@lru_cache(maxsize=1024)
def lap_rad(s, cD=0.0, skin=0.0):
    b = mp.besselk(0, mp.sqrt(s))
    return 1.0 / (s / (b + skin) + (s ** 2) * cD)


@lru_cache(maxsize=None)
def lap_sin(s, *args):
    return 1.0 / (1.0 + s * s)


def to_float_array(result):
    if isinstance(result, np.ndarray):
        return np.array([float(v) for v in result])
    if isinstance(result, list):
        return np.array([float(v) for v in result])
    return np.array([float(result)])


# ---------------------------------------------------------------------------
# Generate all reference data
# ---------------------------------------------------------------------------

def main():
    data = {}

    # Unit step — M=8
    print('Unit step M=8 ...')
    time_unit = 10 ** np.linspace(-1, 1, 101)
    data['time_unit'] = time_unit
    data['unit_m8'] = to_float_array(
        gwr(lambda s: lap_unit(s, 1.0), time_unit, 8)
    )

    # Unit step — M=128
    print('Unit step M=128 ...')
    data['unit_m128'] = to_float_array(
        gwr(lambda s: lap_unit(s, 1.0), time_unit, 128)
    )

    # Ln(t) — M=4
    print('Ln(t) M=4 ...')
    time_ln = 10 ** np.linspace(-2, 5, 31)
    data['time_ln'] = time_ln
    data['ln_m4'] = to_float_array(gwr(lap_ln, time_ln, 4))

    # Radial flow — M=8, cD=1e3
    print('Radial flow M=8 ...')
    time_rad = 10 ** np.linspace(1, 7, 31)
    cD = 1e3
    data['time_rad'] = time_rad
    data['rad_pressure'] = to_float_array(
        gwr(lambda s: lap_rad(s, cD=cD), time_rad, 8)
    )
    data['rad_derivative'] = to_float_array(
        time_rad * gwr(lambda s: s * lap_rad(s, cD=cD), time_rad, 8)
    )
    data['rad_integral'] = to_float_array(
        1 / time_rad * gwr(lambda s: 1 / s * lap_rad(s, cD=cD), time_rad, 8)
    )

    # Hard oscillatory — sin(t) at t=1000
    print('sin(t) at t=1000 ...')
    t0 = np.array([1000.0])
    data['t0_hard'] = t0
    data['sin_expected'] = np.array([float(mp.sin(1000.0))])
    for m in [256, 512, 768, 1024]:
        print(f'  M={m} ...')
        data[f'sin_m{m}'] = to_float_array(gwr(lap_sin, t0, m))

    # Scalar and list input types
    print('Input type checks ...')
    data['scalar_result'] = np.array([float(gwr(lap_ln, 5.0, 4))])
    data['list_result'] = to_float_array(gwr(lap_ln, [1.0, 5.0, 10.0], 4))

    outpath = 'test/reference_data.npz'
    np.savez(outpath, **data)
    print(f'Saved to {outpath}')


if __name__ == '__main__':
    main()
