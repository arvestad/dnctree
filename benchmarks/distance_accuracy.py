#!/usr/bin/env python3
"""
Compares distance-estimation accuracy: dnctree's own estimator vs fastphylo's,
against a known true distance from simulated sequence data. This is a
different question from bench.py (which tests tree-reconstruction accuracy
given an already-noisy distance matrix) - this tests the distance
*estimators* themselves.

Protein (JTT, WAG, LG): dnctree's ml_distance_estimate (modelmatcher +
scipy ML optimization) vs fastphylo.distance_matrix(). Sequence pairs are
simulated via modelmatcher's own sample_count_matrix, so both estimators see
exactly the same simulated data.

DNA (JC, K2P, TN93): dnctree has no independent DNA distance estimator of
its own (it already delegates 100% to fastphylo), so there's no second
"alternative" to compare against fastphylo the way there is for protein.
Instead, this uses the standard closed-form analytical distance formulas
(Jukes-Cantor 1969, Kimura 1980, Tamura-Nei 1993) as an independent
reference, computed from sequences simulated via each model's own
instantiated rate matrix (matrix-exponential substitution simulation, the
same general method modelmatcher uses for protein). Both the simulator and
the formulas are self-validated at the top of this script (large-n_sites
recovery of a known distance) before being used to judge fastphylo.

Usage:
    python benchmarks/distance_accuracy.py
    python benchmarks/distance_accuracy.py -o benchmarks/results/distance_accuracy.json
"""
import argparse
import json
import math
import random
import statistics
import sys

import fastphylo
import numpy as np
from modelmatcher import RateMatrix
from scipy.linalg import expm

from dnctree.distances import ml_distance_estimate

AA_ALPHABET = "ARNDCQEGHILKMFPSTWYV"
BASES = "ACGT"
PURINES = {0, 2}
PYRIMIDINES = {1, 3}

N_SITES = 500
N_REPS = 60
DIST_RANGE = (0.01, 2.0)


# ---------------------------------------------------------------------------
# Protein: dnctree (modelmatcher + ml_distance_estimate) vs fastphylo
# ---------------------------------------------------------------------------

def count_matrix_to_sequences(N, alphabet, rng):
    pairs = []
    for i, a in enumerate(alphabet):
        for j, b in enumerate(alphabet):
            pairs.extend([(a, b)] * int(N[i, j]))
    rng.shuffle(pairs)
    s1 = "".join(p[0] for p in pairs)
    s2 = "".join(p[1] for p in pairs)
    return s1, s2


def run_protein_model(model_name, rng):
    model = RateMatrix.instantiate(model_name)
    rows = []
    for _ in range(N_REPS):
        true_d = rng.uniform(*DIST_RANGE)
        N = model.sample_count_matrix(N_SITES, [true_d])
        s1, s2 = count_matrix_to_sequences(N, AA_ALPHABET, rng)

        d_dnctree = ml_distance_estimate(model, N)

        aln = fastphylo.Alignment([
            fastphylo.Sequence("a", s1, seq_type=fastphylo.SeqType.Protein),
            fastphylo.Sequence("b", s2, seq_type=fastphylo.SeqType.Protein),
        ])
        d_fastphylo = fastphylo.distance_matrix(aln, model=model_name)[0, 1]
        d_dnctree = d_dnctree.real if hasattr(d_dnctree, "real") else d_dnctree

        rows.append((true_d, float(d_dnctree), float(d_fastphylo)))
    return rows


# ---------------------------------------------------------------------------
# DNA: closed-form analytical formulas (independent reference) vs fastphylo
# ---------------------------------------------------------------------------

def is_transition(i, j):
    return (i in PURINES and j in PURINES) or (i in PYRIMIDINES and j in PYRIMIDINES)


def jc_Q():
    Q = np.full((4, 4), 1.0 / 3.0)
    np.fill_diagonal(Q, -1.0)
    return Q, np.array([0.25, 0.25, 0.25, 0.25])


def k2p_Q(kappa=2.0):
    beta = 1.0 / (kappa + 2.0)
    alpha = kappa * beta
    Q = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            if i != j:
                Q[i, j] = alpha if is_transition(i, j) else beta
    for i in range(4):
        Q[i, i] = -Q[i].sum()
    return Q, np.array([0.25, 0.25, 0.25, 0.25])


def tn93_Q(freq, kappa_r=2.5, kappa_y=1.5):
    Q = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            if i in PURINES and j in PURINES:
                Q[i, j] = kappa_r * freq[j]
            elif i in PYRIMIDINES and j in PYRIMIDINES:
                Q[i, j] = kappa_y * freq[j]
            else:
                Q[i, j] = freq[j]
    for i in range(4):
        Q[i, i] = -Q[i].sum()
    mean_rate = -sum(freq[i] * Q[i, i] for i in range(4))
    return Q / mean_rate, np.array(freq)


def simulate_pair(Q, freq, true_d, n_sites, rng):
    P = expm(Q * true_d)
    anc = rng.choices(range(4), weights=freq, k=n_sites)
    s1, s2 = [], []
    for a in anc:
        b = rng.choices(range(4), weights=P[a])[0]
        s1.append(BASES[a])
        s2.append(BASES[b])
    return "".join(s1), "".join(s2)


def jc_distance(s1, s2):
    n = len(s1)
    p = sum(1 for a, b in zip(s1, s2) if a != b) / n
    if p >= 0.75:
        return float("nan")
    return -0.75 * math.log(1 - 4.0 / 3.0 * p)


def k2p_distance(s1, s2):
    idx = {c: i for i, c in enumerate(BASES)}
    n = len(s1)
    n_trans = n_transv = 0
    for a, b in zip(s1, s2):
        if a == b:
            continue
        if is_transition(idx[a], idx[b]):
            n_trans += 1
        else:
            n_transv += 1
    P, Qv = n_trans / n, n_transv / n
    t1, t2 = 1 - 2 * P - Qv, 1 - 2 * Qv
    if t1 <= 0 or t2 <= 0:
        return float("nan")
    return -0.5 * math.log(t1) - 0.25 * math.log(t2)


def tn93_distance(s1, s2, freq):
    gA, gC, gG, gT = freq
    gR, gY = gA + gG, gC + gT
    idx = {c: i for i, c in enumerate(BASES)}
    n = len(s1)
    n_pur = n_pyr = n_transv = 0
    for a, b in zip(s1, s2):
        if a == b:
            continue
        i, j = idx[a], idx[b]
        if i in PURINES and j in PURINES:
            n_pur += 1
        elif i in PYRIMIDINES and j in PYRIMIDINES:
            n_pyr += 1
        else:
            n_transv += 1
    P1, P2, Qv = n_pur / n, n_pyr / n, n_transv / n
    a1_arg = 1 - gR * P1 / (2 * gA * gG) - Qv / (2 * gR)
    a2_arg = 1 - gY * P2 / (2 * gT * gC) - Qv / (2 * gY)
    a3_arg = 1 - Qv / (2 * gR * gY)
    if a1_arg <= 0 or a2_arg <= 0 or a3_arg <= 0:
        return float("nan")
    a1, a2, a3 = -math.log(a1_arg), -math.log(a2_arg), -math.log(a3_arg)
    return (2 * gA * gG / gR * a1 + 2 * gT * gC / gY * a2 +
            2 * (gR * gY - gA * gG * gY / gR - gT * gC * gR / gY) * a3)


DNA_TN93_FREQ = [0.3, 0.2, 0.2, 0.3]

DNA_MODELS = {
    "jc": (jc_Q(), jc_distance, None),
    "k2p": (k2p_Q(kappa=2.0), k2p_distance, None),
    "tn93": (
        tn93_Q(DNA_TN93_FREQ),
        lambda s1, s2: tn93_distance(s1, s2, DNA_TN93_FREQ),
        DNA_TN93_FREQ,
    ),
}


def self_validate_dna(rng):
    """Confirm the simulator + closed-form formulas agree with each other
    (large n_sites, single replicate per distance) before trusting either
    as ground truth to judge fastphylo against."""
    print("=== DNA self-validation (large n_sites; recovering the known simulated distance) ===")
    ok = True
    for name, ((Q, freq), dist_fn, _) in DNA_MODELS.items():
        for true_d in (0.05, 0.3, 0.8, 1.5):
            s1, s2 = simulate_pair(Q, freq, true_d, 200_000, rng)
            est = dist_fn(s1, s2)
            err = est - true_d
            flag = "" if abs(err) < 0.05 else "  <-- suspicious"
            print(f"  {name:5} true_d={true_d:.2f}  estimated={est:.4f}  error={err:+.4f}{flag}")
            if abs(err) >= 0.05:
                ok = False
    if not ok:
        print("WARNING: self-validation error exceeds 0.05 - investigate before trusting results")
    print()
    return ok


def run_dna_model(model_key, rng):
    (Q, freq), dist_fn, freq_for_fp = DNA_MODELS[model_key]
    fp_model_name = model_key  # fastphylo model names are exactly "jc"/"k2p"/"tn93"
    fp_seq_type = fastphylo.SeqType.DNA
    rows = []
    for _ in range(N_REPS):
        true_d = rng.uniform(*DIST_RANGE)
        s1, s2 = simulate_pair(Q, freq, true_d, N_SITES, rng)

        d_formula = dist_fn(s1, s2)

        aln = fastphylo.Alignment([
            fastphylo.Sequence("a", s1, seq_type=fp_seq_type),
            fastphylo.Sequence("b", s2, seq_type=fp_seq_type),
        ])
        d_fastphylo = fastphylo.distance_matrix(aln, model=fp_model_name)[0, 1]

        rows.append((true_d, d_formula, float(d_fastphylo)))
    return rows


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def summarize(rows, label_a, label_b):
    valid = [r for r in rows if not (math.isnan(r[1]) or math.isnan(r[2]))]
    dropped = len(rows) - len(valid)
    errors_a = [r[1] - r[0] for r in valid]
    errors_b = [r[2] - r[0] for r in valid]

    def bias_rmse(errors):
        bias = statistics.mean(errors)
        rmse = math.sqrt(statistics.mean(e * e for e in errors))
        return bias, rmse

    bias_a, rmse_a = bias_rmse(errors_a)
    bias_b, rmse_b = bias_rmse(errors_b)
    return {
        "n": len(valid),
        "dropped_nan": dropped,
        label_a: {"bias": bias_a, "rmse": rmse_a},
        label_b: {"bias": bias_b, "rmse": rmse_b},
    }


def print_summary(model, summary, label_a, label_b):
    a, b = summary[label_a], summary[label_b]
    print(f"{model:6} n={summary['n']:3} (dropped {summary['dropped_nan']:2})  "
          f"{label_a}: bias={a['bias']:+.4f} rmse={a['rmse']:.4f}   "
          f"{label_b}: bias={b['bias']:+.4f} rmse={b['rmse']:.4f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--output", metavar="FILE", help="Write results as JSON to this file")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rng = random.Random(args.seed)

    if not self_validate_dna(rng):
        sys.exit("DNA self-validation failed - aborting rather than reporting bad results")

    results = {}

    print("=== Protein: dnctree (ml_distance_estimate) vs fastphylo ===")
    for model_name in ("JTT", "WAG", "LG"):
        rows = run_protein_model(model_name, rng)
        summary = summarize(rows, "dnctree", "fastphylo")
        results[f"protein_{model_name}"] = {"rows": rows, "summary": summary}
        print_summary(model_name, summary, "dnctree", "fastphylo")

    print("\n=== DNA: closed-form formula (independent reference) vs fastphylo ===")
    for model_key in ("jc", "k2p", "tn93"):
        rows = run_dna_model(model_key, rng)
        summary = summarize(rows, "formula", "fastphylo")
        results[f"dna_{model_key}"] = {"rows": rows, "summary": summary}
        print_summary(model_key.upper(), summary, "formula", "fastphylo")

    if args.output:
        with open(args.output, "w") as fh:
            json.dump({
                "n_sites": N_SITES,
                "n_reps": N_REPS,
                "dist_range": DIST_RANGE,
                "seed": args.seed,
                "results": results,
            }, fh, indent=2)
        print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
