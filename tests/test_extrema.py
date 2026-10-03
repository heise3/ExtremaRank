"""Independent exact-dyadic exhaustive tests of the two-family extrema API."""

import itertools
import math
import random
import unittest
from fractions import Fraction

from extremarank import compare_scores, extrema


def moment(values):
    dyadics = [Fraction.from_float(float(value)) for value in values]
    return sum(dyadics, Fraction()), sum((value * value for value in dyadics), Fraction())


def reference_compare(a, b):
    """Independent exact comparison of signed S/sqrt(Q), with zero Q neutral."""
    sa, qa = a
    sb, qb = b
    sign_a = (sa > 0) - (sa < 0)
    sign_b = (sb > 0) - (sb < 0)
    if sign_a != sign_b:
        return (sign_a > sign_b) - (sign_a < sign_b)
    if sign_a == 0:
        return 0
    cross_a, cross_b = sa * sa * qb, sb * sb * qa
    relation = (cross_a > cross_b) - (cross_a < cross_b)
    return relation if sign_a > 0 else -relation


def exhaustive_reference(values, choose, included):
    candidates = []
    for selected in itertools.combinations(range(len(values)), choose):
        candidates.append((selected, moment(list(included) + [values[j] for j in selected])))
    minimum = maximum = candidates[0]
    for candidate in candidates[1:]:
        if reference_compare(candidate[1], minimum[1]) < 0:
            minimum = candidate
        if reference_compare(candidate[1], maximum[1]) > 0:
            maximum = candidate
    return minimum, maximum


class ExtremaTests(unittest.TestCase):
    def assert_matches_exhaustive(self, values, choose, included=()):
        values, included = list(values), list(included)
        expected_min, expected_max = exhaustive_reference(values, choose, included)
        actual = extrema(values, choose, included=included, method="two_family")
        for score, reference in ((actual.minimum, expected_min), (actual.maximum, expected_max)):
            indices = list(score.indices)
            self.assertEqual(len(indices), choose)
            self.assertEqual(len(indices), len(set(indices)))
            self.assertTrue(all(isinstance(j, int) and 0 <= j < len(values) for j in indices))
            exact = moment(included + [values[j] for j in indices])
            self.assertEqual(score.sum, exact[0])
            self.assertEqual(score.sum_squares, exact[1])
            self.assertEqual(reference_compare(exact, reference[1]), 0)
        self.assertLessEqual(compare_scores(actual.minimum, actual.maximum), 0)
        return actual

    def test_small_integer_cases_all_feasible_cardinalities(self):
        for seed in range(30):
            rng = random.Random(seed)
            values = [rng.randint(-5, 5) for _ in range(4 + seed % 4)]
            included = [rng.randint(-8, 8) for _ in range(seed % 3)]
            for choose in range(max(0, 2 - len(included)), len(values) + 1):
                with self.subTest(seed=seed, choose=choose, included=included):
                    self.assert_matches_exhaustive(values, choose, included)

    def test_exact_binary64_normal_cases_with_forced_inclusion(self):
        for seed in range(18):
            rng = random.Random(3000 + seed)
            values = [rng.normalvariate(0, 3) for _ in range(5 + seed % 3)]
            included = [rng.normalvariate(2, 4) for _ in range(1 + seed % 3)]
            for choose in range(max(0, 2 - len(included)), len(values) + 1):
                with self.subTest(seed=seed, choose=choose):
                    self.assert_matches_exhaustive(values, choose, included)

    def test_all_zero_and_nonzero_constant_vectors(self):
        for value in (0, 3, -3):
            for choose in (0, 1, 2, 4):
                with self.subTest(value=value, choose=choose):
                    result = self.assert_matches_exhaustive([value] * 4, choose, [value, value])
                    self.assertEqual(compare_scores(result.minimum, result.maximum), 0)

    def test_mixed_sign_zero_sum_and_duplicate_ties(self):
        cases = [
            ([-2, -2, 0, 0, 2, 2], []),
            ([-4, -4, -1, 1, 4, 4], [2, -2]),
            ([1, 1, 1, 2, 2, 2], [10]),
            ([-2, -2, -2, -1, -1, -1], [-10]),
            ([1, -1, 2, -2, 0], [0, 0]),
        ]
        for values, included in cases:
            for choose in range(max(0, 2 - len(included)), len(values) + 1):
                with self.subTest(values=values, included=included, choose=choose):
                    self.assert_matches_exhaustive(values, choose, included)

    def test_arbitrary_positive_and_negative_forced_offsets(self):
        for included in ([30], [-30], [40, -10], [-40, 10], [100, -99], [0, 0]):
            for choose in range(max(0, 2 - len(included)), 7):
                with self.subTest(included=included, choose=choose):
                    self.assert_matches_exhaustive([-9, -5, -1, 2, 3, 8], choose, included)

    def test_cancellation_subnormals_and_disparate_scales(self):
        tiny = math.ulp(0.0)
        huge = float(2 ** 500)
        cases = [
            ([float(2 ** 54), 1, -float(2 ** 54), 0, 0.125], [0]),
            ([huge, -huge, math.ulp(huge), 1, -1], [3]),
            ([tiny, 0, -tiny, 2 * tiny, -2 * tiny], [0]),
            ([1, math.nextafter(1.0, 2.0), 1, 2, -1], [1]),
            ([huge, huge + math.ulp(huge), huge, -huge], [huge]),
            ([float(2 ** -500), -float(2 ** -500), 1, -1], [float(2 ** -100)]),
        ]
        for values, included in cases:
            for choose in range(1, len(values) + 1):
                with self.subTest(values=values, included=included, choose=choose):
                    self.assert_matches_exhaustive(values, choose, included)

    def test_score_comparison_matches_independent_cross_products(self):
        vectors = ([0, 0, 0], [1, -1, 0], [1, 1, 1], [2, 2, 2],
                   [-1, -1, -1], [-2, -2, -2], [1, 2, 3], [-1, -2, -3],
                   [float(2 ** 54), 1, -float(2 ** 54)])
        scores = [extrema(values, len(values)).minimum for values in vectors]
        for i, j in itertools.product(range(len(scores)), repeat=2):
            with self.subTest(i=i, j=j):
                self.assertEqual(compare_scores(scores[i], scores[j]),
                                 reference_compare(moment(vectors[i]), moment(vectors[j])))

    def test_invalid_inputs_are_rejected(self):
        cases = [
            ([1, 2], -1, {}), ([1, 2], 3, {}), ([1, 2], 1, {}),
            ([1, 2], True, {}), ([1, 2], 1.5, {}),
            ([float("nan"), 1], 2, {}), ([float("inf"), 1], 2, {}),
            ([1, 2], 2, {"included": [float("nan")]}),
            ([1, 2], 2, {"method": "unknown_method"}),
        ]
        for values, choose, options in cases:
            with self.subTest(values=values, choose=choose, options=options):
                with self.assertRaises((ValueError, TypeError)):
                    extrema(values, choose, **options)


if __name__ == "__main__":
    unittest.main()
