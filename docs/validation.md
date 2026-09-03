# Validation record

## Standalone CI contract

The default, dependency-free CI build runs on Linux and macOS and checks:

- frozen exact-mode outputs for both bundled energy models;
- equality of total accessibility and the five context components;
- `P(w=1) = 1 - P(stem)`;
- probability bounds and nested-window monotonicity;
- deterministic pruning when beam scores tie;
- sequence canonicalization and invalid-input rejection;
- CLI metadata and raw-probability output;
- AddressSanitizer and UndefinedBehaviorSanitizer on Linux.

Run the same checks with `make test` or CMake/CTest.

## Raccess oracle evidence

Before extraction, the optional Raccess configuration was audited in
`TakumiOtagaki/LinearRaccessCapR` through commit `14ba248`. With `beam=0`, both
loop caps at least `N`, exact log-sum-exp, and profile normalization disabled,
the implementation agreed to floating-point precision with:

- an in-process full-span Raccess oracle on crafted structures;
- twelve deterministic random sequences of lengths 12 through 45;
- an independent golden subset generated with the coauthor-maintained Raccess
  CLI at commit `ad29d098`.

The audit also checked probability bounds, context decomposition,
nested-window monotonicity, the single-position identity, whole-sequence
exterior accessibility, object reuse, and cap behavior. This is evidence for
the optional Turner-1999/Raccess configuration; it must not be silently
attributed to the default Turner-2004 standalone backend.

## Approximation and scaling boundary

The production approximations are finite beam width, `c_multi`, `c_hairpin`,
and polynomial log-sum-exp. Existing 500-nt convergence runs show that error
decreases as the beam and caps are relaxed, but no one setting is asserted to
be sufficient for every input class.

A deterministic work audit at `b=200` counted 2,083,773; 8,451,598;
22,842,738; and 54,518,887 inside bifurcation candidates at lengths 250, 500,
1,000, and 2,000. These were 20.8%, 42.3%, 57.1%, and 68.1% of `N b^2`, and no
post-prune table exceeded 200 states. The increasing fraction explains why
short timing series mix beam warm-up with steady-state work.

The planned 30,000-nt HPC runtime and peak-memory experiment remains necessary
for the manuscript's empirical scalability figure. It does not change the
recurrence-level `O(N b^2)` upper bound.
