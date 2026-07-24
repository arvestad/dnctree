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
