"""Declared synthetic scaling of the exact scalar oracle, not biological data."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import platform
import random
import time
from extremarank import PreparedValues, compare_scores


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/scaling.json"))
    args = parser.parse_args()
    rng = random.Random(314159)
    matched = []
    for d in (8, 12, 16, 20):
        # Four independent generated vectors per D; include sorted-mixed and
        # all-positive cases. Same prepared moments and exact comparison cost.
        genes = [[rng.gauss(0 if j % 2 else 1.5, 1) for _ in range(d)] for j in range(4)]
        elapsed = {}
        answers = {}
        for method in ("two_family", "exhaustive"):
            start = time.perf_counter()
            prepared = [PreparedValues(x) for x in genes]
            answers[method] = [p.extrema(d//2, method) for p in prepared]
            elapsed[method] = time.perf_counter() - start
        mismatches = sum(compare_scores(a.minimum, b.minimum) != 0 or
                         compare_scores(a.maximum, b.maximum) != 0
                         for a,b in zip(answers["two_family"], answers["exhaustive"]))
        assert not mismatches
        matched.append({"donors": d, "genes": 4, "retained": d//2,
                        "exact_mismatches": mismatches,
                        "two_family_candidates_per_gene": answers["two_family"][0].candidate_count,
                        "exhaustive_candidates_per_gene": math.comb(d,d//2),
                        "seconds": elapsed,
                        "speedup": elapsed["exhaustive"]/elapsed["two_family"]})
    large = []
    for d in (32, 64, 128):
        genes = [[rng.gauss(1.5 if j % 2 else 0, 1) for _ in range(d)] for j in range(512)]
        start = time.perf_counter()
        answers = [PreparedValues(x).extrema(d//2) for x in genes]
        elapsed = time.perf_counter() - start
        large.append({"donors": d, "genes": len(genes), "retained": d//2,
                      "two_family_seconds": elapsed,
                      "two_family_candidates_per_gene": answers[0].candidate_count,
                      "exhaustive_candidates_per_gene": str(math.comb(d,d//2)),
                      "exhaustive_executed": False})
    payload = {"seed": 314159, "kind": "synthetic scalar-oracle benchmark",
               "python": platform.python_version(), "platform": platform.platform(),
               "matched": matched, "large": large,
               "claim": "Same exact-arithmetic scalar target; speed versus exhaustive enumeration, not approximate DE methods."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2)+"\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
