"""Independent rational oracle for the new Welch target and shared deletions."""
import itertools
import math
import random
import unittest
from fractions import Fraction as F
from functools import cmp_to_key

from extremarank import PreparedWelch, audit_welch, compare_welch
from extremarank.unpaired import rank_welch


def oracle_scores(rows, groups, deleted=()):
    out = []
    for j in range(len(rows[0])):
        means, variances = [], []
        for group in ("T", "R"):
            xs = [F(rows[i][j]) for i, g in enumerate(groups) if g == group and i not in deleted]
            mean = sum(xs)/len(xs)
            means.append(mean)
            variances.append(sum((x-mean)**2 for x in xs)/(len(xs)-1)/len(xs))
        delta = means[0]-means[1]
        out.append(((delta > 0)-(delta < 0), delta**2/sum(variances) if sum(variances) else None, delta))
    return out


def oracle_rank(rows, groups, ids, direction, deleted=()):
    scores = oracle_scores(rows, groups, deleted)
    def compare(a, b):
        sa, qa, _ = scores[a]
        sb, qb, _ = scores[b]
        if direction == "absolute":
            sa, sb = abs(sa), abs(sb)
        if sa != sb:
            c = (sa > sb)-(sa < sb)
        elif not sa:
            c = 0
        elif qa is None or qb is None:
            c = sa*((qa is None)-(qb is None))
        else:
            c = sa*((qa > qb)-(qa < qb))
        if direction == "down":
            c = -c
        return -c if c else (ids[a] > ids[b])-(ids[a] < ids[b])
    return tuple(sorted(range(len(ids)), key=cmp_to_key(compare)))


def feasible_deletions(groups, budget):
    for r in range(1, budget+1):
        for deleted in itertools.combinations(range(len(groups)), r):
            if all(sum(g == group and i not in deleted for i, g in enumerate(groups)) >= 2 for group in ("T", "R")):
                yield deleted


class WelchTests(unittest.TestCase):
    def test_random_exact_ranks_complete_status_and_witness(self):
        rng = random.Random(82177)
        for trial in range(80):
            n, m = rng.randint(2, 5), rng.randint(2, 5)
            groups = ["T"]*n+["R"]*m
            rows = [[rng.randint(-16, 16)/8 for _ in range(5)] for _ in groups]
            ids = ("10", "2", "a", "b", "z")
            budget = min(3, n+m-4)
            deletions = list(feasible_deletions(groups, budget))
            prepared = PreparedWelch(rows, groups, "T", "R")
            self.assertEqual(set(deletions), set(prepared.deletions(budget)))
            self.assertEqual(len(deletions), prepared.deletion_count(budget))
            for direction in ("up", "down", "absolute"):
                baseline = oracle_rank(rows, groups, ids, direction)
                for removed in [(), *deletions]:
                    self.assertEqual(rank_welch(prepared.evaluate(removed), ids, direction),
                                     oracle_rank(rows, groups, ids, direction, removed))
                selected = set(baseline[:2])
                changed = any(set(oracle_rank(rows, groups, ids, direction, d)[:2]) != selected for d in deletions)
                result = audit_welch(rows, groups, ids, "T", "R", 2, budget, 100000, direction)
                self.assertEqual(result.status, "REFUTED" if changed else "CERTIFIED")
                if changed:
                    self.assertIn(result.deleted_indices, deletions)
                    witness = oracle_rank(rows, groups, ids, direction, result.deleted_indices)[:2]
                    self.assertEqual(result.witness_topk, tuple(ids[j] for j in witness))
                    self.assertNotEqual(set(witness), selected)

    def test_extreme_scales_cancellation_zero_variance(self):
        tiny = math.ulp(0.0)
        rows = [[1e300, tiny, 4, 0, 0], [1e300, 2*tiny, 4, 0, 1], [1e300, 3*tiny, 4, 0, -1],
                [-1e300, -tiny, -4, 0, 0], [-1e300, -2*tiny, -4, 0, 1], [-1e300, -3*tiny, -4, 0, -1]]
        groups = ["T"]*3+["R"]*3
        prepared = PreparedWelch(rows, groups, "T", "R")
        for direction in ("up", "down", "absolute"):
            for deleted in [(), *feasible_deletions(groups, 2)]:
                self.assertEqual(rank_welch(prepared.evaluate(deleted), tuple("abcde"), direction),
                                 oracle_rank(rows, groups, tuple("abcde"), direction, deleted))
        self.assertEqual(prepared.evaluate()[3].t, 0)
        self.assertEqual(prepared.evaluate()[2].t, math.inf)

    def test_unresolved_does_not_claim_stable(self):
        rows = [[2, 0], [2, 0], [2, 0], [0, 0], [0, 0], [0, 0]]
        result = audit_welch(rows, ["T"]*3+["R"]*3, ["a", "b"], "T", "R", 1, 2, max_subsets=1)
        self.assertEqual(result.status, "UNRESOLVED")
        self.assertEqual(result.subsets_checked, 1)
        self.assertIsNone(result.deleted_indices)
        self.assertEqual(audit_welch(rows, ["T"]*3+["R"]*3, ["a", "b"], "T", "R", 1, 2).status, "CERTIFIED")

    def test_invalid_designs_and_budget(self):
        with self.assertRaises(ValueError):
            PreparedWelch([[1], [2], [3]], ["T", "R", "R"], "T", "R")
        with self.assertRaises(ValueError):
            PreparedWelch([[1], [2], [3], [math.nan]], ["T", "T", "R", "R"], "T", "R")
        p = PreparedWelch([[1], [2], [3], [4]], ["T", "T", "R", "R"], "T", "R")
        with self.assertRaises(ValueError):
            p.evaluate((0,))
        with self.assertRaises(ValueError):
            audit_welch(p.rows, p.groups, ["a"], "T", "R", 1, 1)

    def test_display_handles_mixed_extreme_scales_within_one_feature(self):
        tiny = math.ulp(0.0)
        rows = [[1e300], [tiny], [2], [-1e300], [-tiny], [-2]]
        groups = ["T"]*3+["R"]*3
        prepared = PreparedWelch(rows, groups, "T", "R")
        self.assertTrue(math.isfinite(prepared.evaluate()[0].t))
        from extremarank import score
        paired = score([1e300, tiny, 2])
        self.assertTrue(math.isfinite(paired.t))
        self.assertTrue(math.isfinite(paired.value))


if __name__ == "__main__":
    unittest.main()
