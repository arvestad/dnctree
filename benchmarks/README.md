# Benchmarks

A repeatable way to check that a change to dnctree's internals is both
**correct** and **faster** — the two things that must both hold before
swapping any hand-rolled implementation for a call into
[fastphylo](https://pypi.org/project/fastphylo/).

Each run builds synthetic model trees (caterpillar shape, fixed sizes),
perturbs the true distances by a fixed amount, reconstructs a tree with a
given algorithm, and measures both wall-clock time and topological accuracy
(Robinson-Foulds distance and tree-matching-distance) against the known true
tree. This is the same correctness methodology `dnctree --alg-testing`
already uses (see `dnctree/algtesting.py`), with timing added and results
saved to a comparable, structured format.

## Usage

```shell
# Run and print a summary table
python benchmarks/bench.py

# Also save results for later comparison
python benchmarks/bench.py -o benchmarks/results/baseline.json --label dnctree

# Once an alternative implementation exists (e.g. calling fastphylo),
# save its results the same way and compare:
python benchmarks/bench.py -o benchmarks/results/fastphylo.json --label fastphylo
python benchmarks/compare.py benchmarks/results/baseline.json benchmarks/results/fastphylo.json
```

`benchmarks/results/baseline.json` is the committed reference point captured
against the pre-fastphylo implementation. Any future migration should be
compared against it: correctness (RF / tree-matching-distance) should stay
the same or improve, and time should decrease.

## Distance-estimator accuracy (`distance_accuracy.py`)

`bench.py` answers "given an already-noisy distance matrix, how well does
the tree-reconstruction algorithm work?" — it never touches sequence-based
distance *estimation* at all (`TestingPartialDistanceMatrix` reads true
distances straight from a model tree and adds synthetic noise). A different
question — "how accurately does a distance estimator recover the true
evolutionary distance from real sequence data?" — needs its own tool:

```shell
python benchmarks/distance_accuracy.py
python benchmarks/distance_accuracy.py -o benchmarks/results/distance_accuracy.json
```

For each protein model (JTT, WAG, LG), sequence pairs are simulated via
`modelmatcher`'s own `sample_count_matrix` at a known true distance, and
both dnctree's `ml_distance_estimate` and `fastphylo.distance_matrix` are
run on the identical simulated data. For DNA (JC, K2P, TN93) — where
dnctree has no independent estimator of its own, since it already delegates
100% to fastphylo — the standard closed-form analytical formulas
(Jukes-Cantor 1969, Kimura 1980, Tamura-Nei 1993) serve as an independent
reference instead, computed from sequences simulated via each model's own
rate matrix (matrix-exponential substitution). The DNA simulator and
formulas self-validate (recover a known distance at large `n_sites`) before
being trusted to judge fastphylo — the script aborts rather than reporting
if that check fails.

## Real-data end-to-end speed (`real_data_speed.py`)

Neither of the above times a full, real run: `bench.py` reconstructs trees
from synthetic already-computed distances (no real distance estimation
happens at all), and `distance_accuracy.py` measures estimator accuracy,
not speed. This one runs `divide_n_conquer_tree` end-to-end (both
algorithms) on the real alignments already checked into `testdata/`, at a
few different widths:

```shell
python benchmarks/real_data_speed.py
python benchmarks/real_data_speed.py -o benchmarks/results/real_data_speed.json
```

This is the one that actually caught the ~3.9x speedup from fastphylo's
`RateMatrix` caching fix (see the `fastphylo-py` repo) — a real, standing
check for future fastphylo/dnctree changes on real data, not just that one.
