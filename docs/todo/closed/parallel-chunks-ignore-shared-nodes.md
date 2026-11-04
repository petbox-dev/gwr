# Parallel chunks ignore GWR nodes shared between time points

Severity: low
Opened: 2026-10-31
Closed: 2026-10-31 -- settled by: measured; node-aware chunking computes fewer Bessel values but runs slower than interleaving, so interleaving stays.
If this was wrong, you will see: on an evenly spaced time grid, `gwr(..., workers=4)` loses more than about 10% of its single-process speedup to repeated Bessel evaluations, with a chunking that balances load as well as interleaving does.

## What happens

The GWR nodes are `s = k ln2 / t`, k = 1..2M. Two time points share a node when
their ratio is a fraction with numerator and denominator up to 2M. The rounded
mpf nodes coincide bit for bit, so a process that evaluates both times computes
the shared Bessel values once (order-0/1 cache). `_map_interleaved` gives each
worker `times[i::n]`, so coincident nodes on different workers are computed twice.

## Evidence

Nodes computed exactly as `gwr_single` does (gmpy2 and mpmath backends agree to
about 1%), 4 workers:

| Grid | Nodes that coincide | Duplicate work, interleaved |
|---|---|---|
| 21 points, log-spaced 0.01-200 | 0% | 0% |
| 61 points, `10 ** linspace(1, 7, 61)` | 1.3-2.3% | 1.3-2.4% |
| 60 points, linear 1-60 | 34-35% | 30-31% |
| 100 points, linear 0.5-50 | 34-36% | 31-33% |

A greedy partitioner (each time to the worker whose node set grows least, at
most ceil(n/w) times per worker) cut the duplicate work on the linear grids to
7-20% (4 workers) and to 0 on the log grids. Timed with real workers,
`RadialNoFlow(50)` plus storage, 100 linear points, best of 3:

| M | Interleaved | Greedy |
|---|---|---|
| 16 | 0.33-0.41 s | 0.36-0.41 s |
| 32 | 1.00-1.17 s | 1.32-1.36 s |

On 61 log points the two were within noise. The cost of a time point varies
along the curve (the boundary ratio and the asymptotic switch depend on s), and
interleaving balances that; the greedy chunks did not. The Bessel values at
sqrt(s) are about 40% of a call (profile, M = 32), so removing a third of them
bounds the gain near 10% even with perfect balance.

## Ruled out

- A shared cross-process cache: one order-0/1 ARB evaluation is 75-175 us; a
  manager-dict lookup costs tens of microseconds and every miss adds a write.
- A two-stage plan (distinct nodes first, then seed the workers' caches): the
  parent would have to predict each model's Bessel arguments (sqrt of s_eff
  under dual porosity), so it works only for built-in geometries.

## Files

- `gwr_inversion/_parallel.py`: `_map_interleaved`, `_interleaved_chunks`
