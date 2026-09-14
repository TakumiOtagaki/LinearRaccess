# LinearRaccess

LinearRaccess computes the probability that every nucleotide in each requested
RNA window is unpaired. It combines a beam-pruned inside/outside dynamic
program with range-add posterior aggregation, avoiding a global base-pair-span
cutoff.

This repository is the standalone algorithm, command-line tool, and C++ API
extracted from `TakumiOtagaki/LinearRaccessCapR` at commit `14ba248`. It does
not contain ENBPD, Python bindings, or a vendored Raccess source tree.

## Build and test

LinearRaccess requires a C++17 compiler, CMake 3.20 or later, and Python 3 for
the dependency-free CLI smoke test.

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

The Make build is also supported:

```bash
make -j4
make test
```

## Command line

Given `input.fa`:

```text
>example
GCGAAACGC
```

run:

```bash
./build/LinRacc \
  -seqfile=input.fa \
  -outfile=accessibility.txt \
  -access_len=1,3,7 \
  -beam=200 \
  -metadata \
  -source-revision=$(git rev-parse --short HEAD)
```

The normal output is Raccess-compatible `-RT log(P)` accessibility energy.
Use `-probabilities` for raw probabilities or `-byloop` for the exterior,
hairpin, bulge, internal-loop, and multiloop contributions. Run
`./build/LinRacc --help` for all options.

For every reported result, record at least the energy backend/model, beam,
`c_multi`, `c_hairpin`, log-sum-exp mode, normalization setting, and requested
window lengths. `-metadata` writes these fields into the output.

### Build provenance

`LinRacc --build-info` prints the software version and embedded source revision
as JSON. CMake refreshes the revision header on every build, including after a
commit without source edits. A dirty checkout is marked with `-dirty`; a source
archive without Git metadata reports `unknown`. Probability output metadata uses
this embedded revision by default. `-source-revision` remains an explicit output
label override and does not change `--build-info`.

For reproducible production jobs, build from a clean commit and record both the
embedded revision and the executable's SHA-256. Do not infer an old executable's
revision from the current checkout.

## C++ API

```cpp
#include <linearraccess/api.hpp>

lcr::api::LinearRaccessConfig config;
config.beam = 200;
auto result = lcr::api::linear_raccess("GCGAAACGC", {1, 3}, config);
double p = result.accessibility[1][0];  // P(window [0, 2] is unpaired)
```

The installed CMake target is `LinearRaccess::linearraccess`.

## Energy backends

The default standalone build uses the LinearCapR energy adapter and Turner
2004 parameters. `-energy=turner1999` selects the legacy model.

The Raccess compatibility backend is optional because the upstream Raccess
license in the development repository prohibits redistribution. This
repository contains only the adapter written for LinearRaccess; it neither
downloads nor vendors Raccess. If you have an independently obtained source
tree and permission for your intended use, enable the adapter explicitly:

```bash
cmake -S . -B build-raccess \
  -DLINEARRACCESS_WITH_RACCESS=ON \
  -DRACCESS_INCLUDE_DIR=/path/to/raccess/src
cmake --build build-raccess -j4
```

The output backend is then selected with `-engine=raccess`. The default stays
`lincapr` so that identical command lines do not silently change model when an
optional dependency happens to be present.

## Algorithmic contract and limitations

For fixed beam width `b`, loop caps, and number of requested window lengths,
the production work is `O(N b^2)` and retained DP memory is `O(N b)`. The
quadratic-in-beam term is the retained `M2 x M1` bifurcation. Short prefixes
under-fill the beam, so a finite timing fit need not show the asymptotic upper
bound directly.

Four settings can change numerical results: beam pruning, `c_multi`,
`c_hairpin`, and the polynomial fast log-sum-exp. Set `beam=0`, raise both loop
caps to at least the sequence length, use `-no-fast-logsumexp`, and use
`-no-normalize` for exact short-sequence validation. Raising loop caps with
sequence length deliberately leaves the fixed-cap complexity regime and is
not intended for long inputs.

See [docs/algorithm.md](docs/algorithm.md) and
[docs/validation.md](docs/validation.md) for the method and evidence boundary.

## Citation

The LinearRaccess manuscript is in preparation. The underlying LinearCapR
method is:

> Otagaki T, Hosokawa H, Fukunaga T, Iwakiri J, Terai G, Asai K.
> LinearCapR: linear-time computation of per-nucleotide structural-context
> probabilities of RNA without base-pair span limits. *Bioinformatics* 42(6),
> btag295 (2026). https://doi.org/10.1093/bioinformatics/btag295

## License and provenance

Original LinearRaccess code is available under the MIT license in `LICENSE`.
Parameter tables and the fast log-sum-exp routine retain their upstream terms;
see `THIRD_PARTY_NOTICES.md` and `third_party/`. The ViennaRNA 2.5.1 terms
include attribution and no-fee redistribution conditions, so the repository
must not be described as uniformly MIT-licensed.
