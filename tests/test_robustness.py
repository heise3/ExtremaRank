"""Independent complete subsets verify certificates, bounds and minimum changes."""
import itertools
import random
import unittest
from fractions import Fraction as F

from extremarank.diagnostics import PreparedPaired
from extremarank.robustness import SearchProblem, audit_feature_robustness
from extremarank.unpaired import PreparedWelch, compare_welch, audit_welch
from test_unpaired import oracle_rank, feasible_deletions


def paired_rank(rows, ids, direction, deleted=()):
    # Ordinary paired t equals Welch t against an all-zero reference group,
    # whose variance contribution is zero. This independently reuses centered
    # rational variance rather than the production dyadic moments/oracle.
    remaining = [row for i, row in enumerate(rows) if i not in deleted]
    matrix = remaining+[[0]*len(ids)]*2
    return oracle_rank(matrix, ["T"]*len(remaining)+["R"]*2, ids, direction)


def verify_problem(rows, groups, ids, direction, design, k, budget):
    prepared = PreparedPaired(rows) if design == "paired" else PreparedWelch(rows, groups, "T", "R")
    ranker = (lambda deleted: paired_rank(rows, ids, direction, deleted)) if design == "paired" else (
        lambda deleted: oracle_rank(rows, groups, ids, direction, deleted))
    baseline = ranker(())
    original = set(baseline[:k])
    baseline_sign = []
    for j in range(len(ids)):
        if design == "paired":
            value = sum(F(row[j]) for row in rows)
        else:
            value = sum(F(row[j]) for row, g in zip(rows, groups) if g == "T")/groups.count("T")-sum(
                F(row[j]) for row, g in zip(rows, groups) if g == "R")/groups.count("R")
        baseline_sign.append((value > 0)-(value < 0))
    deleted_sets = list(feasible_deletions(groups, budget)) if design == "welch" else [d for r in range(1, budget+1) for d in itertools.combinations(range(len(rows)), r)]
    minima = {(j, prop): None for j in range(len(ids)) for prop in ("membership", "direction")}
    for deleted in deleted_sets:
        selected = set(ranker(deleted)[:k])
        for j in range(len(ids)):
            if design == "paired":
                value = sum(F(row[j]) for i, row in enumerate(rows) if i not in deleted)
            else:
                value = sum(F(row[j]) for i, (row, g) in enumerate(zip(rows, groups)) if g == "T" and i not in deleted)/sum(
                    g == "T" and i not in deleted for i, g in enumerate(groups))-sum(F(row[j]) for i, (row, g) in enumerate(zip(rows, groups)) if g == "R" and i not in deleted)/sum(
                    g == "R" and i not in deleted for i, g in enumerate(groups))
            sign = (value > 0)-(value < 0)
            for prop, changed in (("membership", (j in selected) != (j in original)), ("direction", sign != baseline_sign[j])):
                if changed:
                    old = minima[j, prop]
                    minima[j, prop] = min(old, len(deleted)) if old else len(deleted)
    report = audit_feature_robustness(prepared, ids, k, budget, design, direction, ids, 100000, 100000)
    checks = 0
    for row in report["features"]:
        j = ids.index(row["feature_id"])
        for prop in ("membership", "direction"):
            result, minimum = row[prop], minima[j, prop]
            assert result["status"] == ("REFUTED" if minimum else "CERTIFIED"), (row, minimum)
            assert result["exact_minimum_change"] == minimum, (row, minimum)
            if minimum:
                deleted = tuple(result["witness"]["deleted_indices"])
                assert len(deleted) == minimum
                if prop == "membership":
                    assert (j in set(ranker(deleted)[:k])) != (j in original)
            else:
                assert result["minimum_change_lower_bound"] == budget+1
            checks += 1
    if design == "welch":
        result = audit_welch(rows, groups, ids, "T", "R", k, budget, 100000, direction, "branch-bound", 100000)
        assert result.status == report["status"]
    # Check every root enclosure against all feasible actual scalar scores.
    problem = SearchProblem(prepared, ids, design, direction, k)
    for r in range(1, budget+1):
        for counts in problem.allocations(r):
            bounds = problem.bounds((), tuple(range(len(rows))), counts)
            children = []
            for fixed_deleted in ((),(1,)):
                optional=tuple(i for i in range(1,len(rows)) if i not in fixed_deleted)
                needed=(counts[0]-1,*counts[1:])
                if problem.canonical((0,),optional,needed) is not None:
                    children.append((fixed_deleted,problem.bounds((0,),optional,needed)))
            for deleted in deleted_sets:
                if len(deleted) != r or any(sum(problem.labels[i] == group and i not in deleted for i in range(len(rows))) != count for group, count in enumerate(counts)):
                    continue
                for score, (lo, hi) in zip(problem.evaluate(deleted), bounds):
                    assert compare_welch(lo, score) <= 0 and compare_welch(score, hi) <= 0
                    checks += 1
                if 0 not in deleted:
                    for fixed_deleted,child_bounds in children:
                        if not set(fixed_deleted)<=set(deleted):continue
                        for score,(lo,hi) in zip(problem.evaluate(deleted),child_bounds):
                            assert compare_welch(lo,score)<=0 and compare_welch(score,hi)<=0
                            checks+=1
    # A limited search may refute early or certify via a bound. Its interval
    # must enclose the independently enumerated minimum in every case.
    for limit in (1,3):
        capped=audit_feature_robustness(prepared,ids,k,budget,design,direction,ids,limit,limit,witness_trials=1)
        for row in capped['features']:
            j=ids.index(row['feature_id'])
            for prop in ('membership','direction'):
                actual=minima[j,prop];result=row[prop]
                assert result['minimum_change_lower_bound'] <= (actual if actual else budget+1)
                if result['status']=='CERTIFIED':assert actual is None
                if result['minimum_change_upper_bound'] is not None:
                    assert actual is not None and actual<=result['minimum_change_upper_bound']
                if result['exact_minimum_change'] is not None:assert result['exact_minimum_change']==actual
                checks+=1
    return checks


class RobustnessTests(unittest.TestCase):
    def test_extreme_scales_cancellation_and_subnormals(self):
        xs = [1e300, -1e300, 5e-324, -5e-324, 0.0, 1e-200, -1e-200, 1.0]
        rows = [[x, float(2**40+i), (i-4)*5e-324, 0.0] for i,x in enumerate(xs)]
        groups = ["T"]*4+["R"]*4
        for design in ("paired", "welch"):
            for direction in ("up", "down", "absolute"):
                verify_problem(rows, groups, tuple("abcd"), direction, design, 2, 2)

    def test_independent_complete_small_problems(self):
        rng = random.Random(60319)
        for _ in range(35):
            groups = ["T"]*rng.randint(2, 4)+["R"]*rng.randint(2, 4)
            rows = [[rng.randint(-20, 20)/8 for j in range(4)] for g in groups]
            for design in ("paired", "welch"):
                for direction in ("up", "down", "absolute"):
                    verify_problem(rows, groups, tuple("abcd"), direction, design, 2, min(2, len(rows)-4))

    def test_incomplete_search_bounds_and_zero_budget(self):
        rows = [[1, 0], [2, 1], [3, -1], [4, 0], [-1, 2], [3, 8]]
        p = PreparedPaired(rows)
        report = audit_feature_robustness(p, ("a", "b"), 1, 3, max_nodes=1, max_scenarios=1, witness_trials=0)
        for row in report["features"]:
            for prop in ("membership", "direction"):
                self.assertIn(row[prop]["status"], ("UNRESOLVED", "CERTIFIED"))
                self.assertIsNone(row[prop]["exact_minimum_change"])
        self.assertEqual(audit_feature_robustness(p, ("a", "b"), 1, 0)["status"], "CERTIFIED")

    def test_stable_welch_uses_bounds(self):
        groups = ["T"]*10+["R"]*10
        rows = [[100+i/10, 1+(i%3)/10, 0] for i in range(10)]+[[i/10, 1+(i%3)/10, 0] for i in range(10)]
        result = audit_welch(rows, groups, ("a", "b", "c"), "T", "R", 1, 3,
                             10000, "up", "branch-bound", 10000)
        self.assertEqual(result.status, "CERTIFIED")
        self.assertLess(result.subsets_checked, result.feasible_deletion_sets)
        self.assertGreater(result.query_bounds_pruned, 0)

    def test_indexed_many_feature_bounds_match_direct_counts(self):
        rng=random.Random(9001)
        rows=[[rng.randint(-10,10)/8 for j in range(60)] for i in range(8)]
        ids=tuple('g'+str(i) for i in range(60));groups=['T']*4+['R']*4
        for design in ('paired','welch'):
            prepared=PreparedPaired(rows) if design=='paired' else PreparedWelch(rows,groups,'T','R')
            for direction in ('up','down','absolute'):
                problem=SearchProblem(prepared,ids,design,direction,12)
                for retained in problem.allocations(2):
                    for included,optional,needed in [((),tuple(range(8)),retained),((0,),tuple(range(1,8)),(retained[0]-1,*retained[1:]))]:
                        bounds=problem.bounds(included,optional,needed)
                        queries=[(j,prop) for j in range(60) for prop in ('membership','direction')]
                        indexed=problem.prover(bounds,queries)
                        self.assertEqual([indexed(q) for q in queries],[problem.proves(q,bounds) for q in queries])


if __name__ == "__main__":
    unittest.main()
