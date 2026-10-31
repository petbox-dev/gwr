"""Tests for gwr_batch() and evaluate(): several functions or models on one time grid."""
from __future__ import annotations

from collections.abc import Iterator
from typing import Any

import numpy as np
import pytest
from mpmath import mp

from gwr_inversion import (
    Model, RadialInfinite, RadialNoFlow, Skin, WellboreStorage, _flint, evaluate, gwr,
    gwr_batch, shutdown_pool,
)

TIMES = [0.01, 0.1, 1.0, 10.0, 100.0]
M = 8


@pytest.fixture(autouse=True)
def restore_precision() -> Iterator[None]:
    with mp.workdps(mp.dps):
        yield


def _exponential(s: Any) -> Any:
    """L{e^-t}; module-level for pickling."""
    return 1 / (s + 1)


def _ln(s: Any) -> Any:
    """L{ln t + gamma}; module-level for pickling."""
    return -mp.log(s) / s


def _models() -> list[Model]:
    return [Model(RadialNoFlow(radius), Skin(1.0), WellboreStorage(10.0))
            for radius in (5.0, 20.0, 100.0)]


def _rows_equal(batch: np.ndarray[Any, Any], rows: list[Any]) -> bool:
    return all(np.array_equal(np.asarray(batch[i], dtype=object), np.asarray(row, dtype=object))
               for i, row in enumerate(rows))


# ---------------------------------------------------------------- gwr_batch

@pytest.mark.parametrize("time", [TIMES, np.array(TIMES), tuple(TIMES)])
def test_rows_are_bit_identical_to_gwr(time: Any) -> None:
    fns = [_exponential, _ln]
    batch = gwr_batch(fns, time, M)
    assert batch.shape == (2, len(TIMES)) and batch.dtype == np.float64
    assert _rows_equal(batch, [gwr(fn, list(TIMES), M) for fn in fns])


def test_mpf_rows_are_bit_identical_to_gwr() -> None:
    batch = gwr_batch([_exponential, _ln], TIMES, M, dtype="mpf")
    assert batch.dtype == object and batch.shape == (2, len(TIMES))
    for row, fn in zip(batch, [_exponential, _ln]):
        expected = gwr(fn, TIMES, M, dtype="mpf")
        assert all(isinstance(v, mp.mpf) and v == e for v, e in zip(row, expected))


def test_scalar_and_zero_dimensional_time() -> None:
    scalar = gwr_batch([_exponential, _ln], 1.0, M)
    assert scalar.shape == (2,)
    assert list(scalar) == [gwr(_exponential, 1.0, M), gwr(_ln, 1.0, M)]
    zero_d = gwr_batch([_exponential, _ln], np.array(1.0), M)
    assert zero_d.shape == (2, 1)
    assert np.array_equal(zero_d[:, 0], scalar)


def test_parallel_rows_are_bit_identical_to_sequential() -> None:
    try:
        sequential = gwr_batch([_exponential, _ln], TIMES, M)
        parallel = gwr_batch([_exponential, _ln], TIMES, M, workers=2)
    finally:
        shutdown_pool()
    assert np.array_equal(sequential, parallel)


def test_refine_matches_gwr_row_by_row() -> None:
    times = list(np.linspace(0.5, 12.0, 12))   # e^-t falls below 1e-4 of its peak
    batch = gwr_batch([_exponential, _ln], times, M, refine=True)
    for row, fn in zip(batch, [_exponential, _ln]):
        assert np.array_equal(row, np.asarray(gwr(fn, times, M, refine=True)))
    assert not np.array_equal(batch[0], gwr_batch([_exponential], times, M)[0])


def test_empty_time_grid() -> None:
    assert gwr_batch([_exponential, _ln], [], M).shape == (2, 0)


@pytest.mark.parametrize("fns, error, message", [
    (_exponential, TypeError, "for one function use gwr"),
    ([], ValueError, "at least one function"),
    ([_exponential, 3.0], TypeError, "must be callable"),
    (42, TypeError, "sequence of functions"),
])
def test_function_validation(fns: Any, error: type, message: str) -> None:
    with pytest.raises(error, match=message):
        gwr_batch(fns, TIMES, M)


def test_dtype_validation() -> None:
    with pytest.raises(ValueError, match="dtype"):
        gwr_batch([_exponential], TIMES, M, dtype="float32")


# ---------------------------------------------------------------- evaluate

@pytest.mark.parametrize("output", ["pressure", "rate", "cumulative", "derivative"])
def test_evaluate_rows_are_bit_identical_to_model_methods(output: str) -> None:
    models = _models()
    batch = evaluate(models, TIMES, output, M)
    assert batch.shape == (len(models), len(TIMES))
    assert _rows_equal(batch, [getattr(model, output)(TIMES, M) for model in models])


def test_evaluate_derivative_with_scalar_and_array_time() -> None:
    models = _models()
    scalar = evaluate(models, 2.0, "derivative", M)
    assert list(scalar) == [model.derivative(2.0, M) for model in models]
    array = evaluate(models, np.array(TIMES), "derivative", M)
    assert _rows_equal(array, [model.derivative(np.array(TIMES), M) for model in models])


def test_evaluate_reads_a_generator_once() -> None:
    models = _models()
    batch = evaluate(models, (t for t in TIMES), "derivative", M)
    assert _rows_equal(batch, [model.derivative(TIMES, M) for model in models])


def test_evaluate_in_parallel_matches_sequential() -> None:
    models = _models()
    try:
        parallel = evaluate(models, TIMES, "pressure", M, workers=2)
    finally:
        shutdown_pool()
    assert np.array_equal(parallel, evaluate(models, TIMES, "pressure", M))


@pytest.mark.parametrize("models, output, options, error, message", [
    (Model(RadialInfinite()), "pressure", {}, TypeError, "sequence of models"),
    (["not a model"], "pressure", {}, TypeError, "must be a Model"),
    (_models(), "hybrid_rate", {}, ValueError, "output must be one of"),
    (_models(), "derivative", {"dtype": "mpf"}, ValueError, "dtype='float' only"),
])
def test_evaluate_validation(models: Any, output: str, options: dict[str, Any],
                             error: type, message: str) -> None:
    with pytest.raises(error, match=message):
        evaluate(models, TIMES, output, M, **options)


@pytest.mark.skipif(not _flint.HAS_FLINT, reason="python-flint not installed")
def test_models_in_one_batch_share_the_bessel_values_at_the_nodes() -> None:
    # Bounded radial models that differ in outer radius ask for the same
    # K0, K1, I0, I1 at sqrt(s); the second and third models add no evaluations.
    models = _models()
    _flint.cache_clear()
    evaluate(models[:1], TIMES, "pressure", M)
    one = _flint._bessel_order_01_cached.cache_info().misses
    _flint.cache_clear()
    evaluate(models, TIMES, "pressure", M)
    assert _flint._bessel_order_01_cached.cache_info().misses == one


@pytest.mark.skipif(not _flint.HAS_FLINT, reason="python-flint not installed")
def test_batch_shares_nodes_with_a_cache_of_one_time_point(monkeypatch: Any) -> None:
    # Time-major order needs only one time point's 2M nodes in the cache; model by
    # model, the same cache size forgets every node before the next model needs it.
    from functools import lru_cache
    small = lru_cache(maxsize=2 * M)(_flint._bessel_order_01_cached.__wrapped__)
    monkeypatch.setattr(_flint, "_bessel_order_01_cached", small)
    models = _models()
    evaluate(models[:1], TIMES, "pressure", M)
    one = small.cache_info().misses
    small.cache_clear()
    evaluate(models, TIMES, "pressure", M)
    assert small.cache_info().misses == one
    small.cache_clear()
    for model in models:
        model.pressure(TIMES, M)
    assert small.cache_info().misses == len(models) * one
