#!/usr/bin/env python3
"""
Benchmark harness comparing dnc_tree_multi's two non-core-taxa partition
methods (quartet vs tree) against each other and against the existing
core-tree (dnc_tree_k) baseline, on both correctness (RF distance against a
known model tree) and speed - same methodology as bench.py, extended with
dnc_tree_multi's extra knobs (core_size, partition_method).

Usage:
    python benchmarks/multi_partition.py
    python benchmarks/multi_partition.py -o benchmarks/results/multi_partition.json
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

from dnctree.stupid_tree import caterpillar, even_steven
from dnctree.algtesting import TestingPartialDistanceMatrix, testing_divide_n_conquer
from dnctree.version import __version__ as dnctree_version


# (label, tree_builder, n_taxa) -- tree_builder() returns a newick string for
# a fresh model tree at fixed branch length 0.1. Both shapes are
# deterministic given their size/depth, so the only randomness is in the
# distance perturbation added below (seeded per rep) and each algorithm's
# own core-taxa sampling. even_steven's depth->n_taxa is 3*2^(depth-2), so
# depths 5/6/7 (24/48/96 taxa) closely match the caterpillar sizes below for
# a fair side-by-side across tree shapes.
SCENARIOS = [
    ("caterpillar-20", lambda: caterpillar(20, branchlength=0.1), 20),
    ("caterpillar-50", lambda: caterpillar(50, branchlength=0.1), 50),
    ("caterpillar-100", lambda: caterpillar(100, branchlength=0.1), 100),
    ("even_steven-24", lambda: even_steven(5, branchlength=0.1), 24),
    ("even_steven-48", lambda: even_steven(6, branchlength=0.1), 48),
    ("even_steven-96", lambda: even_steven(7, branchlength=0.1), 96),
]

# (label, algorithm, base_case_size, core_size, partition_method)
CONFIGS = [
    ("core-tree", "core-tree", 10, 10, None),
    ("multi-quartet", "multi", 10, 10, "quartet"),
    ("multi-tree", "multi", 10, 10, "tree"),
]

EPSILON = 0.05
ERROR_DISTRIBUTION = "uniform"
N_REPS = 30


def git_commit():
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True, cwd="."
        )
        return out.stdout.strip()
    except Exception:
        return None


def run_one(model_tree_newick, algorithm, base_case_size, core_size, partition_method, seed):
    random.seed(seed)
    dm = TestingPartialDistanceMatrix(model_tree_newick, EPSILON, ERROR_DISTRIBUTION)

    start = time.perf_counter()
    t = testing_divide_n_conquer(
        dm,
        algorithm=algorithm,
        base_case_size=base_case_size,
        core_size=core_size,
        partition_method=partition_method or "quartet",  # ignored unless algorithm == "multi"
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
    for scenario_label, tree_builder, n_taxa in SCENARIOS:
        model_tree_newick = tree_builder()
        for config_label, algorithm, base_case_size, core_size, partition_method in CONFIGS:
            times, rfs, tms = [], [], []
            for rep in range(N_REPS):
                elapsed, rf_normalized, tm = run_one(
                    model_tree_newick, algorithm, base_case_size, core_size, partition_method, seed=rep
                )
                times.append(elapsed)
                rfs.append(rf_normalized)
                tms.append(tm)
            results.append({
                "scenario": scenario_label,
                "n_taxa": n_taxa,
                "config": config_label,
                "algorithm": algorithm,
                "base_case_size": base_case_size,
                "core_size": core_size,
                "partition_method": partition_method,
                "n_reps": N_REPS,
                "time_seconds": summarize(times),
                "rf_normalized": summarize(rfs),
                "tree_matching_distance": summarize(tms),
            })
    return results


def print_summary(results):
    header = f"{'scenario':16} {'config':14} {'time (s)':>10} {'RF (norm)':>10} {'TM dist':>10}"
    print(header)
    print("-" * len(header))
    for r in results:
        print(f"{r['scenario']:16} {r['config']:14} "
              f"{r['time_seconds']['mean']:>10.4f} "
              f"{r['rf_normalized']['mean']:>10.3f} "
              f"{r['tree_matching_distance']['mean']:>10.2f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--output", metavar="FILE",
                     help="Also write results as JSON to this file.")
    ap.add_argument("--label", default="dnctree",
                     help="Label for this run. Default: %(default)s")
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
