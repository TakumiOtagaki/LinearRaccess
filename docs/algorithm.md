# Algorithm

## Output quantity

For an RNA sequence of length `N`, a requested window length `w`, and a
zero-based start `i`, LinearRaccess reports

```text
P_w(i) = Pr(all positions in [i, i+w-1] are unpaired).
```

An unpaired window lies in exactly one unpaired segment of an exterior,
hairpin, bulge, internal, or multiloop context. Therefore

```text
P_w(i) = P_w^E(i) + P_w^H(i) + P_w^B(i) + P_w^I(i) + P_w^M(i).
```

The `-byloop` output exposes these five terms.

## Computation

1. A left-to-right inside pass computes log-space partition weights for the
   exterior, stem, stem-end, multiloop, multiloop-bifurcation, and auxiliary
   states.
2. At each position, every sparse state table retains the best `b` candidates
   under a prefix-aware score. Equal scores are resolved by state index, so
   pruning is deterministic.
3. A right-to-left outside pass traverses the retained hypergraph and computes
   posterior event weights.
4. Each unpaired segment event contributes to every requested window it
   contains. Difference arrays apply the contribution to a contiguous range of
   window starts in constant time per event, followed by one prefix sum.

`beam=0` disables pruning and is intended only for short exact tests. Unlike a
maximum-span approximation, beam pruning can retain long-range base pairs.

## Complexity

At one sequence position, the dominant bifurcation combines at most `b`
retained `M2` states with at most `b` retained `M1` states. The inside and
outside traversals therefore have an `O(N b^2)` upper bound. By default
`c_multi = c_hairpin = b`, so multiloop unpaired runs cost `O(N b^2)` and
hairpin candidates `O(N b)`; with `MAXLOOP` (30) and the number of requested
window lengths fixed, every transition family stays linear in `N`. Retained
tables use `O(N b)` memory; output storage is linear in `N` per requested
window length.

With `beam=0` the default caps become the sequence length, i.e. the exact,
uncapped model; that setting is for validation and is outside the linear
bound.

## Numerical modes

The default polynomial log-sum-exp is faster but approximate. Exact validation
uses `-no-fast-logsumexp`; raw posterior sums are preserved rather than clamped
so numerical-bound violations remain visible. Profile normalization is another
explicit option and is disabled for oracle comparisons.
