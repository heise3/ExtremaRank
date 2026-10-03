#!/usr/bin/env python3
"""Independent exhaustive exact-dyadic verification of the candidate theorem.

No implementation import, NumPy, or approximate numeric comparison is used.
Binary64 values are lifted to integer units with a common power-of-two scale;
the scale cancels from signed S/sqrt(Q). All feasible subsets are enumerated.
"""

import argparse
from datetime import datetime, timezone
import itertools
import json
import math
from pathlib import Path
import platform
import random
import sys
import time


def compare(a, b):
    sa, qa = a
    sb, qb = b
    sign_a, sign_b = (sa > 0) - (sa < 0), (sb > 0) - (sb < 0)
    if sign_a != sign_b:
        return (sign_a > sign_b) - (sign_a < sign_b)
    if sign_a == 0:
        return 0
    ca, cb = sa * sa * qb, sb * sb * qa
    relation = (ca > cb) - (ca < cb)
    return relation if sign_a > 0 else -relation


def lift_to_integer_units(values):
    ratios = [float(value).as_integer_ratio() for value in values]
    denominator = max((denominator for _, denominator in ratios), default=1)
    return [numerator * (denominator // scale) for numerator, scale in ratios]


def range_for_subsets(values, a, b, subsets):
    low = high = None
    minimum_sum = maximum_sum = None
    for selected in subsets:
        score = (a + sum(values[j] for j in selected),
                 b + sum(values[j] * values[j] for j in selected))
        if low is None or compare(score, low) < 0:
            low = score
        if high is None or compare(score, high) > 0:
            high = score
        minimum_sum = score[0] if minimum_sum is None else min(minimum_sum, score[0])
        maximum_sum = score[0] if maximum_sum is None else max(maximum_sum, score[0])
    return low, high, minimum_sum, maximum_sum


def verify(cases, seed):
    rng = random.Random(seed)
    cache = {}
    coverage = {"integer": 0, "normal_binary64": 0, "nonempty_fixed": 0,
                "all_positive_sums": 0, "all_negative_sums": 0,
                "mixed_or_zero_sums": 0, "all_zero_substitutions": 0,
                "positive_constant_substitutions": 0,
                "negative_constant_substitutions": 0,
                "disparate_scale_substitutions": 0}
    enumerated = scanned = 0
    failures = []
    started = time.perf_counter()
    for index in range(cases):
        total = 5 + index % 5
        fixed_count = 1 + (index // 5) % 3
        undecided_count = total - fixed_count
        choose = 1 + (index // 15) % undecided_count
        integer_family = bool(index % 2)
        if integer_family:
            fixed = [float(rng.randint(-8, 8)) for _ in range(fixed_count)]
            values = [float(rng.randint(-8, 8)) for _ in range(undecided_count)]
            coverage["integer"] += 1
        else:
            fixed = [rng.normalvariate(0, 5) for _ in range(fixed_count)]
            values = [rng.normalvariate(0, 5) for _ in range(undecided_count)]
            coverage["normal_binary64"] += 1
        if index % 7 == 0:
            fixed = [value + 25 for value in fixed]
        elif index % 7 == 1:
            fixed = [value - 25 for value in fixed]
        if index % 997 == 0:
            fixed, values = [0.0] * fixed_count, [0.0] * undecided_count
            coverage["all_zero_substitutions"] += 1
        elif index % 991 == 0:
            fixed, values = [3.0] * fixed_count, [3.0] * undecided_count
            coverage["positive_constant_substitutions"] += 1
        elif index % 983 == 0:
            fixed, values = [-3.0] * fixed_count, [-3.0] * undecided_count
            coverage["negative_constant_substitutions"] += 1
        elif index % 977 == 0:
            exponents = (-500, -100, 0, 100, 500)
            fixed = [math.ldexp(value, exponents[j % len(exponents)])
                     for j, value in enumerate(fixed)]
            values = [math.ldexp(value, exponents[(j + 2) % len(exponents)])
                      for j, value in enumerate(values)]
            coverage["disparate_scale_substitutions"] += 1
        values.sort()
        lifted = lift_to_integer_units(fixed + values)
        fi, ui = lifted[:fixed_count], lifted[fixed_count:]
        a, b = sum(fi), sum(value * value for value in fi)
        key = (undecided_count, choose)
        if key not in cache:
            all_subsets = list(itertools.combinations(range(undecided_count), choose))
            family_subsets = set()
            for left in range(undecided_count - choose + 1):
                family_subsets.add(tuple(range(left, left + choose)))
            for prefix in range(choose + 1):
                family_subsets.add(tuple(list(range(prefix)) + list(range(
                    undecided_count - (choose - prefix), undecided_count))))
            cache[key] = all_subsets, sorted(family_subsets)
        all_subsets, family_subsets = cache[key]
        exact_min, exact_max, minimum_sum, maximum_sum = range_for_subsets(ui, a, b, all_subsets)
        family_min, family_max, _, _ = range_for_subsets(ui, a, b, family_subsets)
        enumerated += len(all_subsets)
        scanned += len(family_subsets)
        coverage["nonempty_fixed"] += 1
        label = ("all_positive_sums" if minimum_sum > 0 else
                 "all_negative_sums" if maximum_sum < 0 else "mixed_or_zero_sums")
        coverage[label] += 1
        if compare(exact_min, family_min) != 0 or compare(exact_max, family_max) != 0:
            failures.append({"case": index, "included": [value.hex() for value in fixed],
                             "values": [value.hex() for value in values], "choose": choose})
            break
    return {
        "verification_passed": not failures,
        "seed": seed, "requested_cases": cases,
        "completed_cases": sum(coverage[key] for key in ("integer", "normal_binary64")),
        "coverage": coverage, "failure_count": len(failures), "failures": failures,
        "exhaustive_subset_evaluations": enumerated,
        "candidate_subset_evaluations": scanned,
        "arithmetic": "exact binary64 dyadics lifted to common integer units; exact cross products",
        "implementation_imported": False,
        "sampling_note": "stdlib Random; independently rerunnable 200k check, not the earlier NumPy stream",
        "retained_count_minimum": 2,
        "total_donor_counts": [5, 6, 7, 8, 9], "fixed_donor_counts": [1, 2, 3],
        "elapsed_seconds": time.perf_counter() - started,
        "python": sys.version, "platform": platform.platform(),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "claim_limit": "finite exact verification supports theorem review; does not replace proof or establish priority",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=int, default=200000)
    parser.add_argument("--seed", type=int, default=20261004)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1]
                        / "results" / "independent_theorem_check.json")
    args = parser.parse_args()
    if args.cases < 1:
        parser.error("cases must be positive")
    result = verify(args.cases, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in
                      ("verification_passed", "completed_cases", "failure_count",
                       "coverage", "elapsed_seconds")}, indent=2))
    if result["failure_count"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
