#!/usr/bin/env python3
"""
Compare two benchmark result files produced by bench.py (e.g. a dnctree
baseline vs. a fastphylo-backed run) on both axes: speed and correctness.

Usage:
    python benchmarks/compare.py baseline.json candidate.json
"""
import argparse
import json


def load(path):
    with open(path) as fh:
        return json.load(fh)


def index_by_key(results):
    return {(r["scenario"], r["algorithm"]): r for r in results}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("baseline")
    ap.add_argument("candidate")
    args = ap.parse_args()

    baseline = load(args.baseline)
    candidate = load(args.candidate)

    print(f"Baseline:  {baseline['label']} ({baseline.get('git_commit')})")
    print(f"Candidate: {candidate['label']} ({candidate.get('git_commit')})")
    print()

    base_by_key = index_by_key(baseline["results"])
    cand_by_key = index_by_key(candidate["results"])

    header = (f"{'scenario':16} {'algorithm':10} {'speedup':>9} "
              f"{'base RF':>9} {'cand RF':>9} {'base TM':>9} {'cand TM':>9}")
    print(header)
    print("-" * len(header))

    for key, base in base_by_key.items():
        cand = cand_by_key.get(key)
        scenario, algorithm = key
        if cand is None:
            print(f"{scenario:16} {algorithm:10} {'(missing in candidate)':>9}")
            continue
        base_t = base["time_seconds"]["mean"]
        cand_t = cand["time_seconds"]["mean"]
        speedup = base_t / cand_t if cand_t > 0 else float("inf")
        print(f"{scenario:16} {algorithm:10} {speedup:>8.2f}x "
              f"{base['rf_normalized']['mean']:>9.3f} {cand['rf_normalized']['mean']:>9.3f} "
              f"{base['tree_matching_distance']['mean']:>9.2f} {cand['tree_matching_distance']['mean']:>9.2f}")

    missing_in_base = set(cand_by_key) - set(base_by_key)
    for key in missing_in_base:
        scenario, algorithm = key
        print(f"{scenario:16} {algorithm:10} (missing in baseline)")


if __name__ == "__main__":
    main()
