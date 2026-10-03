"""Predeclared synthetic shared-deletion search ablation; no biological claims.

The cases and order below are fixed in source. All methods use the same branch
order, node budget, and witness_trials=0: only the node bound changes. Ground
truth independently enumerates common deletion sets with Fraction arithmetic.
Run: PYTHONPATH=src python benchmarks/search_ablation.py
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import platform
import random
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction
from functools import cmp_to_key
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from extremarank import audit_topk  # noqa: E402


def reference_rank(rows, ids, keep, k, direction="up"):
    """Independent exact signed-ratio comparison; no implementation imports."""
    moments = []
    for j in range(len(ids)):
        values = [Fraction.from_float(float(rows[i][j])) for i in keep]
        moments.append((sum(values, Fraction()), sum((v*v for v in values), Fraction())))

    def compare(a, b):
        sa, qa = moments[a]
        sb, qb = moments[b]
        sign_a, sign_b = (sa > 0) - (sa < 0), (sb > 0) - (sb < 0)
        if direction == "absolute":
            ratio_a = sa*sa/qa if qa else Fraction()
            ratio_b = sb*sb/qb if qb else Fraction()
            c = (ratio_a > ratio_b) - (ratio_a < ratio_b)
        elif sign_a != sign_b:
            c = (sign_a > sign_b) - (sign_a < sign_b)
        elif sign_a == 0:
            c = 0
        else:
            left, right = sa*sa*qb, sb*sb*qa
            c = sign_a * ((left > right) - (left < right))
        if direction == "down":
            c = -c
        return -c if c else ((ids[a] > ids[b]) - (ids[a] < ids[b]))

    return tuple(ids[j] for j in sorted(range(len(ids)), key=cmp_to_key(compare))[:k])


def exhaustive_truth(rows, ids, k, budget, direction="up"):
    d = len(rows)
    baseline = reference_rank(rows, ids, tuple(range(d)), k, direction)
    changes, examined, first = 0, 0, None
    for count in range(budget+1):
        for removed in itertools.combinations(range(d), count):
            removed_set = set(removed)
            keep = tuple(i for i in range(d) if i not in removed_set)
            ranking = reference_rank(rows, ids, keep, k, direction)
            examined += 1
            if set(ranking) != set(baseline):
                changes += 1
                if first is None:
                    first = {"deleted_indices": list(removed), "topk": list(ranking)}
    return {"status": "REFUTED" if changes else "CERTIFIED",
            "baseline_topk": list(baseline), "subsets_examined": examined,
            "changed_topk_subsets": changes, "first_witness": first}


def witness_is_valid(rows, ids, k, budget, result, direction="up"):
    if result.status != "REFUTED":
        return result.witness_topk is None and result.deleted_indices is None
    removed = result.deleted_indices
    if removed is None or not 1 <= len(removed) <= budget or len(set(removed)) != len(removed):
        return False
    if any(not 0 <= i < len(rows) for i in removed):
        return False
    keep = tuple(i for i in range(len(rows)) if i not in set(removed))
    actual = reference_rank(rows, ids, keep, k, direction)
    return actual == result.witness_topk and set(actual) != set(result.baseline_topk)


def cases():
    """Fixed synthetic construction, disclosed rather than tuned to results."""
    out = []
    for d, g, budget in ((12, 8, 2), (16, 12, 3)):
        scale_pattern = [float(2**(i//2)) for i in range(d)]
        rows = []
        for i in range(d):
            # Two finite-variance near-constant winners, variable outsiders.
            row = [9.0 + (-0.5, 0.0, 0.5)[i % 3],
                   7.0 + (-0.25, 0.25)[i % 2]]
            row += [scale_pattern[(i+2*j) % d] for j in range(g-2)]
            rows.append(row)
        out.append({"name": f"stable_wide_moment_box_D{d}_G{g}", "rows": rows,
                    "ids": [f"g{j:02d}" for j in range(g)], "k": 2, "budget": budget,
                    "construction": "two low-variance positive winners; cyclic powers-of-two outsiders"})
    for d, g, seed in ((12, 8, 1208), (16, 12, 1612)):
        rng = random.Random(seed)
        rows = [[rng.randint(-6, 9) / 4 for _ in range(g)] for _ in range(d)]
        out.append({"name": f"mixed_sign_D{d}_G{g}", "rows": rows,
                    "ids": [f"m{j:02d}" for j in range(g)], "k": 2, "budget": 2,
                    "construction": f"random.Random({seed}), integers [-6,9]/4; no case selection"})
    a, h = (2.0, 2.1, 1.9, 10.0), (1.0, 1.1, 0.9, 2.0)
    out.append({"name": "variance_reversal_D4_G2", "rows": list(map(list, zip(a, h))),
                "ids": ["a", "h"], "k": 1, "budget": 1,
                "construction": "a>h in every donor; deleting donor index 3 reverses signed paired-t ranking"})
    out.append({"name": "separate_adversaries_overlap_D4_G2",
                "rows": [[2.0*x, x] for x in (1.0, 2.0, 3.0, 4.0)],
                "ids": ["a", "b"], "k": 1, "budget": 2,
                "construction": "proportional columns: every shared subset ties, lexical a always wins"})
    source = out[0]
    out.append({**source, "name": "stable_down_D12_G8", "direction": "down",
                "rows": [[-x for x in row] for row in source["rows"]],
                "construction": "sign reversal of fixed stable D12 construction; down-ranking contract"})
    out.append({**source, "name": "stable_absolute_D12_G8", "direction": "absolute",
                "rows": [[-x if j < 2 else x for j, x in enumerate(row)] for row in source["rows"]],
                "construction": "negative near-constant winners and positive variable outsiders; absolute-ranking contract"})
    return out


def display_t(values):
    xs = [Fraction.from_float(float(x)) for x in values]
    n, s, q = len(xs), sum(xs, Fraction()), sum((x*x for x in xs), Fraction())
    return math.copysign(math.sqrt(float((n-1)*s*s/(n*q-s*s))), s)


def main():
    records = []
    for case in cases():
        rows, ids, k, budget = (case[name] for name in ("rows", "ids", "k", "budget"))
        direction = case.get("direction", "up")
        truth_start = time.perf_counter()
        truth = exhaustive_truth(rows, ids, k, budget, direction)
        truth_seconds = time.perf_counter() - truth_start
        methods = {}
        for method in ("two_family", "moment_box", "none"):
            start = time.perf_counter()
            result = audit_topk(rows, ids, k, budget, max_nodes=100000,
                                bound_method=method, witness_trials=0, direction=direction)
            seconds = time.perf_counter() - start
            valid = (result.status == truth["status"] and
                     list(result.baseline_topk) == truth["baseline_topk"] and
                     witness_is_valid(rows, ids, k, budget, result, direction))
            methods[method] = {**result.as_dict(), "wall_seconds": seconds,
                               "matches_exhaustive_truth": valid}
            if not valid:
                raise AssertionError((case["name"], method, truth, result.as_dict()))
        encoded = json.dumps(rows, separators=(",", ":")).encode()
        record = {"name": case["name"], "kind": "synthetic", "construction": case["construction"],
                  "shape": [len(rows), len(ids)], "gene_ids": ids, "k": k,
                  "max_deletions": budget, "direction": direction,
                  "effects_sha256_json": hashlib.sha256(encoded).hexdigest(),
                  "effects": rows, "truth": {**truth, "wall_seconds": truth_seconds},
                  "methods": methods}
        record["two_family_node_reduction_vs_moment_box"] = methods["moment_box"]["nodes"] / methods["two_family"]["nodes"]
        record["two_family_node_reduction_vs_none"] = methods["none"]["nodes"] / methods["two_family"]["nodes"]
        if case["name"].startswith("variance_reversal"):
            record["explanation"] = {
                "all_donors_a_greater_than_h": all(row[0] > row[1] for row in rows),
                "baseline_t": {ids[j]: display_t([row[j] for row in rows]) for j in range(2)},
                "delete_index_3_t": {ids[j]: display_t([row[j] for row in rows[:3]]) for j in range(2)},
                "all_delete_one_variances_positive": all(
                    len(set(rows[i][j] for i in range(4) if i != removed)) > 1
                    for removed in range(4) for j in range(2))}
        records.append(record)
        print(case["name"], truth["status"],
              {name: {"nodes": value["nodes"], "seconds": round(value["wall_seconds"], 6)}
               for name, value in methods.items()})
    report = {"run_at_utc": datetime.now(timezone.utc).isoformat(),
              "environment": {"python": sys.version, "python_executable": sys.executable,
                              "platform": platform.platform()},
              "contract": {"case_selection": "fixed deterministic source constructions; no biological data",
                           "max_nodes": 100000, "witness_trials": 0,
                           "only_ablation_change": "bound_method", "timing_repeats": 1,
                           "timing_limit": "single local run; illustrative, not a hardware-normalized performance claim"},
              "all_methods_match_exhaustive_truth": True, "cases": records}
    destination = ROOT / "results" / "search_ablation.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print("Wrote", destination)


if __name__ == "__main__":
    main()
