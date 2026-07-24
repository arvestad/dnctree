#!/usr/bin/env python3
"""
End-to-end speed on real biological data, using the files already checked
into testdata/. Neither of the other two benchmarks cover this: bench.py
times tree reconstruction against synthetic, already-computed distances
(TestingPartialDistanceMatrix never touches real distance estimation at
all), and distance_accuracy.py measures estimator accuracy, not speed.
This is the one that would have caught the fastphylo RateMatrix caching
win (~3.9x on divide_n_conquer_tree end-to-end) - useful as a standing
check for future fastphylo/dnctree changes, not just this one.

Usage:
    python benchmarks/real_data_speed.py
    python benchmarks/real_data_speed.py -o benchmarks/results/real_data_speed.json
"""
import argparse
import json
import statistics
import time

from alv.io import guess_format, read_alignment

from dnctree.api import divide_n_conquer_tree
from dnctree.msa import MSA
from dnctree.version import __version__ as dnctree_version

try:
    import fastphylo
    FASTPHYLO_VERSION = fastphylo.__version__
except ImportError:
    FASTPHYLO_VERSION = None

TESTDATA_DIR = "testdata"
FILES = ["s63_L50.phylip", "s6_L100.phylip", "s83_L500.phylip"]
N_REPS = 5


def load_msa(filename):
    path = f"{TESTDATA_DIR}/{filename}"
    fmt = guess_format(path)
    alv_msa, _ = read_alignment(path, "guess", fmt, None, None)
    return MSA.from_alv(alv_msa)


def time_it(fn, reps=N_REPS):
    times = [_time_once(fn) for _ in range(reps)]
    return {"mean": statistics.mean(times), "min": min(times), "max": max(times)}


def _time_once(fn):
    t0 = time.perf_counter()
    fn()
    return time.perf_counter() - t0


def run():
    results = []
    for filename in FILES:
        msa = load_msa(filename)
        n_taxa = len(msa.taxa())

        for label, simple_alg in (("core-tree", False), ("simple", True)):
            timing = time_it(lambda: divide_n_conquer_tree(msa, simple_alg=simple_alg))
            results.append({
                "file": filename,
                "n_taxa": n_taxa,
                "msa_width": msa.msa_width,
                "algorithm": label,
                "time_seconds": timing,
            })
    return results


def print_summary(results):
    header = f"{'file':18} {'n_taxa':>6} {'width':>6} {'algorithm':10} {'time (s)':>10}"
    print(header)
    print("-" * len(header))
    for r in results:
        print(f"{r['file']:18} {r['n_taxa']:>6} {r['msa_width']:>6} {r['algorithm']:10} "
              f"{r['time_seconds']['mean']:>10.4f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--output", metavar="FILE", help="Write results as JSON to this file")
    args = ap.parse_args()

    results = run()
    print_summary(results)

    if args.output:
        with open(args.output, "w") as fh:
            json.dump({
                "dnctree_version": dnctree_version,
                "fastphylo_version": FASTPHYLO_VERSION,
                "n_reps": N_REPS,
                "results": results,
            }, fh, indent=2)
        print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
