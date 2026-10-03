"""Independent common-subset ground truth for the complete ranking audit."""
import itertools
import math
import random
import unittest
from fractions import Fraction
from functools import cmp_to_key

from extremarank import audit_topk, extrema, score, compare_scores


def reference_rank(rows, ids, keep, k, direction="up"):
    moments = []
    for j in range(len(ids)):
        xs = [Fraction.from_float(float(rows[i][j])) for i in keep]
        moments.append((sum(xs, Fraction()), sum((x*x for x in xs), Fraction())))

    def compare(a, b):
        sa, qa = moments[a]
        sb, qb = moments[b]
        sga, sgb = (sa > 0) - (sa < 0), (sb > 0) - (sb < 0)
        if direction == "absolute":
            ratio_a = sa*sa/qa if qa else Fraction()
            ratio_b = sb*sb/qb if qb else Fraction()
            relation = (ratio_a > ratio_b) - (ratio_a < ratio_b)
        elif sga != sgb:
            relation = (sga > sgb) - (sga < sgb)
        elif not sga:
            relation = 0
        else:
            left, right = sa*sa*qb, sb*sb*qa
            relation = sga * ((left > right) - (left < right))
        if direction == "down":
            relation = -relation
        return -relation if relation else ((ids[a] > ids[b]) - (ids[a] < ids[b]))

    return tuple(ids[j] for j in sorted(range(len(ids)), key=cmp_to_key(compare))[:k])


def exhaustive_truth(rows, ids, k, budget, direction="up"):
    baseline = reference_rank(rows, ids, tuple(range(len(rows))), k, direction)
    changed = []
    for r in range(1, budget+1):
        for deleted in itertools.combinations(range(len(rows)), r):
            keep = tuple(i for i in range(len(rows)) if i not in set(deleted))
            actual = reference_rank(rows, ids, keep, k, direction)
            if set(actual) != set(baseline):
                changed.append((deleted, actual))
    return baseline, changed


class AuditTests(unittest.TestCase):
    def check_against_truth(self, rows, ids, k, budget, witness_trials=0, direction="up"):
        baseline, changed = exhaustive_truth(rows, ids, k, budget, direction)
        for method in ("two_family", "moment_box", "none"):
            with self.subTest(method=method, witness_trials=witness_trials):
                result = audit_topk(rows, ids, k, budget, max_nodes=100000,
                                    bound_method=method, witness_trials=witness_trials, direction=direction)
                self.assertEqual(result.baseline_topk, baseline)
                self.assertEqual(result.status, "REFUTED" if changed else "CERTIFIED")
                if changed:
                    removed = result.deleted_indices
                    self.assertIsNotNone(removed)
                    self.assertTrue(1 <= len(removed) <= budget)
                    self.assertEqual(tuple(sorted(set(removed))), removed)
                    self.assertTrue(all(0 <= i < len(rows) for i in removed))
                    keep = tuple(i for i in range(len(rows)) if i not in set(removed))
                    actual = reference_rank(rows, ids, keep, k, direction)
                    self.assertEqual(actual, result.witness_topk)
                    self.assertNotEqual(set(actual), set(baseline))
                else:
                    self.assertIsNone(result.deleted_indices)
                    self.assertIsNone(result.witness_topk)

    def test_small_matrices_exhaustive_common_subset_truth(self):
        for seed in range(48):
            rng = random.Random(701+seed)
            d, g = 4+seed % 4, 2+seed % 4
            rows = [[rng.randint(-8, 8)/4 for _ in range(g)] for _ in range(d)]
            ids = [f"gene_{g-j}" for j in range(g)]
            k, budget = 1+seed % (g-1), min(3, d-2)
            with self.subTest(seed=seed, shape=(d, g), k=k, budget=budget):
                self.check_against_truth(rows, ids, k, budget)

    def test_witness_proposals_require_actual_common_subset(self):
        for seed in range(12):
            rng = random.Random(1300+seed)
            rows = [[rng.normalvariate(0, 2) for _ in range(4)] for _ in range(6)]
            with self.subTest(seed=seed):
                self.check_against_truth(rows, ["z", "a", "c", "b"], 2, 2, witness_trials=8)

    def test_down_and_absolute_match_common_subset_exhaustive_truth(self):
        for direction in ("down", "absolute"):
            for seed in range(36):
                rng = random.Random(1900+seed)
                d, g = 4+seed % 4, 2+seed % 4
                rows = [[rng.randint(-12, 12)/8 for _ in range(g)] for _ in range(d)]
                with self.subTest(direction=direction, seed=seed):
                    self.check_against_truth(rows, [f"g{g-j}" for j in range(g)],
                                             1+seed % (g-1), min(3, d-2),
                                             witness_trials=8 if seed % 2 else 0, direction=direction)

    def test_absolute_opposite_sign_ties_and_crossing_zero(self):
        cases = [
            [[x, -x, 0] for x in (1, 2, 3, 4)],
            [[-2, 2, 0], [1, -1, 2], [1, -1, -2], [0, 0, 0]],
            [[1, -1, 0], [1, -1, 0], [1, -1, 0], [1, -1, 0]],
            [[0, 0, 0] for _ in range(4)],
        ]
        for direction in ("down", "absolute"):
            for index, rows in enumerate(cases):
                with self.subTest(direction=direction, case=index):
                    self.check_against_truth(rows, ["z", "a", "b"], 1, 2,
                                             witness_trials=8, direction=direction)

    def test_all_zero_exact_gene_id_ties(self):
        rows, ids = [[0.0]*4 for _ in range(5)], ["z", "a", "aa", "b"]
        self.check_against_truth(rows, ids, 2, 3)
        for method in ("two_family", "moment_box", "none"):
            result = audit_topk(rows, ids, 2, 3, bound_method=method, witness_trials=8)
            self.assertEqual(result.baseline_topk, ("a", "aa"))
            self.assertEqual(result.status, "CERTIFIED")

    def test_proportional_columns_overlap_is_not_refutation(self):
        rows = [[2*x, x] for x in (1.0, 2.0, 3.0, 4.0)]
        env_a, env_b = extrema([r[0] for r in rows], 3), extrema([r[1] for r in rows], 3)
        self.assertLess(compare_scores(env_a.minimum, env_b.maximum), 0)
        self.assertNotEqual(env_a.minimum.indices, env_b.maximum.indices)
        self.check_against_truth(rows, ["a", "b"], 1, 2, witness_trials=8)

    def test_search_limit_reports_unresolved_without_witness(self):
        rows = [[2*x, x] for x in (1.0, 2.0, 3.0, 4.0)]
        for method in ("two_family", "moment_box", "none"):
            with self.subTest(method=method):
                result = audit_topk(rows, ["a", "b"], 1, 1, max_nodes=1,
                                    bound_method=method, witness_trials=0)
                self.assertEqual(result.status, "UNRESOLVED")
                self.assertEqual(result.nodes, 1)
                self.assertEqual(result.retained_counts_completed, (4,))
                self.assertIsNone(result.deleted_indices)
                self.assertIsNone(result.witness_topk)

    def test_offset_cancellation_subnormal_and_high_dynamic_range(self):
        tiny, huge = math.ulp(0.0), float(2**500)
        matrices = [
            [[100+x, -100+x, x] for x in (-9, -5, -1, 2, 3, 8)],
            [[2**54, 1, 4], [1, 1.25, 3], [-2**54, -1, 2], [0, 0.5, 1], [0.125, 2, 0]],
            [[tiny, 0, -tiny], [2*tiny, tiny, 0], [-tiny, -2*tiny, tiny], [0, 0, 0]],
            [[huge, huge, 1], [-huge, huge+math.ulp(huge), -1],
             [math.ulp(huge), huge, 0.5], [1, -huge, 2]],
        ]
        for direction in ("up", "down", "absolute"):
            for index, rows in enumerate(matrices):
                with self.subTest(case=index, direction=direction):
                    self.check_against_truth(rows, ["a", "b", "c"], 1,
                                             min(2, len(rows)-2), direction=direction)

    def test_variance_reversal_with_all_rowwise_mean_margins_positive(self):
        rows = [[2, 1], [2.1, 1.1], [1.9, 0.9], [10, 2]]
        self.assertTrue(all(a > h for a, h in rows))
        self.check_against_truth(rows, ["a", "h"], 1, 1, witness_trials=8)
        baseline_t = [score([row[j] for row in rows]).t for j in range(2)]
        after_t = [score([row[j] for row in rows[:3]]).t for j in range(2)]
        self.assertLess(baseline_t[0], baseline_t[1])
        self.assertGreater(after_t[0], after_t[1])

    def test_zero_budget_and_all_features_selected(self):
        rows = [[1, -1], [2, 2], [3, 1]]
        for k, budget in ((1, 0), (2, 1)):
            with self.subTest(k=k, budget=budget):
                result = audit_topk(rows, ["a", "b"], k, budget, max_nodes=1)
                self.assertEqual(result.status, "CERTIFIED")
                self.assertEqual(result.nodes, 0)

    def test_invalid_inputs(self):
        valid = dict(effects=[[1, 2], [2, 3], [3, 4]], gene_ids=["a", "b"], k=1, max_deletions=1)
        cases = [{"effects": []}, {"effects": [[1, 2], [3]]},
                 {"gene_ids": ["a", "a"]}, {"k": True}, {"k": 0},
                 {"max_deletions": True}, {"max_deletions": 2},
                 {"max_nodes": True}, {"max_nodes": 0}, {"witness_trials": -1},
                 {"bound_method": "bad"}, {"direction": "bad"},
                 {"effects": [[1, 2], [float("nan"), 3], [3, 4]]},
                 {"effects": [[1, 2], [float("inf"), 3], [3, 4]]}]
        for patch in cases:
            with self.subTest(patch=patch):
                with self.assertRaises((TypeError, ValueError)):
                    audit_topk(**(valid | patch))


if __name__ == "__main__":
    unittest.main()
