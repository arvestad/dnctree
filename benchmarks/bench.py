#!/usr/bin/env python3
"""
Benchmark harness for dnctree: measures correctness (against a known model
tree) and speed together, so that a future implementation (e.g. one calling
into fastphylo) can be compared against this baseline on both axes at once.

Usage:
    python benchmarks/bench.py                        # run and print a summary
    python benchmarks/bench.py -o results/baseline.json  # also save results to JSON
"""
import argparse
import json
import random
import statistics
import subprocess
import sys
import time
import warnings

import ete3

# tree_matching_distance's own __init__.py uses the deprecated pkg_resources
# API for its version lookup; harmless, but noisy every run.
warnings.filterwarnings("ignore", message=r"pkg_resources is deprecated.*", category=UserWarning)
from tree_matching_distance import distance as matching_distance

from dnctree.stupid_tree import caterpillar
from dnctree.algtesting import TestingPartialDistanceMatrix, testing_divide_n_conquer
from dnctree.version import __version__ as dnctree_version


# (label, n_taxa) -- caterpillar trees are deterministic given n_taxa, so the
# only randomness is in the distance perturbation added below, which we seed.
SCENARIOS = [
    ("caterpillar-20", 20),
    ("caterpillar-50", 50),
    ("caterpillar-100", 100),
]

# (label, algorithm, base_case_size) -- "nj" reuses the core-tree code path
# with base_case_size == n_taxa, which forces the pure weighted-NJ base case
# with no divide-and-conquer recursion (same trick main.py's
# --alg-testing-nj uses).
ALGORITHMS = [
    ("nj", "core-tree", None),       # base_case_size filled in with n_taxa
    ("core-tree", "core-tree", 10),
    ("simple", "simple", 10),
]

EPSILON = 0.05
ERROR_DISTRIBUTION = "uniform"
N_REPS = 5


def git_commit():
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True, cwd="."
        )
        return out.stdout.strip()
    except Exception:
        return None


def run_one(model_tree_newick, n_taxa, algorithm, base_case_size, seed):
    random.seed(seed)
    dm = TestingPartialDistanceMatrix(model_tree_newick, EPSILON, ERROR_DISTRIBUTION)

    start = time.perf_counter()
    t = testing_divide_n_conquer(
        dm,
        algorithm=algorithm,
        base_case_size=base_case_size if base_case_size is not None else n_taxa,
    )
    elapsed = time.perf_counter() - start

    estimated_tree = ete3.Tree(str(t))
    rf, rf_max, *_ = estimated_tree.robinson_foulds(dm._dt, unrooted_trees=True)
    rf_normalized = rf / rf_max if rf_max else 0.0
    tm = matching_distance(estimated_tree, dm._dt)

    return elapsed, rf_normalized, tm


def summarize(values):
    return {
        "mean": statistics.mean(values),
        "min": min(values),
        "max": max(values),
        "stdev": statistics.stdev(values) if len(values) > 1 else 0.0,
    }


def run_benchmarks():
    results = []
    for scenario_label, n_taxa in SCENARIOS:
        model_tree_newick = caterpillar(n_taxa, branchlength=0.1)
        for algo_label, algorithm, base_case_size in ALGORITHMS:
            times, rfs, tms = [], [], []
            for rep in range(N_REPS):
                elapsed, rf_normalized, tm = run_one(
                    model_tree_newick, n_taxa, algorithm, base_case_size, seed=rep
                )
                times.append(elapsed)
                rfs.append(rf_normalized)
                tms.append(tm)
            results.append({
                "scenario": scenario_label,
                "n_taxa": n_taxa,
                "algorithm": algo_label,
                "base_case_size": base_case_size if base_case_size is not None else n_taxa,
                "n_reps": N_REPS,
                "time_seconds": summarize(times),
                "rf_normalized": summarize(rfs),
                "tree_matching_distance": summarize(tms),
            })
    return results


def print_summary(results):
    header = f"{'scenario':16} {'algorithm':10} {'time (s)':>10} {'RF (norm)':>10} {'TM dist':>10}"
    print(header)
    print("-" * len(header))
    for r in results:
        print(f"{r['scenario']:16} {r['algorithm']:10} "
              f"{r['time_seconds']['mean']:>10.4f} "
              f"{r['rf_normalized']['mean']:>10.3f} "
              f"{r['tree_matching_distance']['mean']:>10.2f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--output", metavar="FILE",
                     help="Also write results as JSON to this file.")
    ap.add_argument("--label", default="dnctree",
                     help="Label for this run (e.g. 'dnctree' or 'fastphylo'). Default: %(default)s")
    args = ap.parse_args()

    results = run_benchmarks()
    print_summary(results)

    if args.output:
        payload = {
            "label": args.label,
            "dnctree_version": dnctree_version,
            "git_commit": git_commit(),
            "python_version": sys.version.split()[0],
            "epsilon": EPSILON,
            "error_distribution": ERROR_DISTRIBUTION,
            "results": results,
        }
        with open(args.output, "w") as fh:
            json.dump(payload, fh, indent=2)
        print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
