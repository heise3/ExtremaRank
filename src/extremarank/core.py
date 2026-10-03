"""Exact subset extrema of signed self-normalized paired effects.

Inputs are converted to finite IEEE-754 binary64 values. Arithmetic used to
select extrema and compare scores is exact for those values; float display
values are never used in a decision. See docs/THEOREM.md for the two-family
characterization. This is a descriptive ranking audit, not a DE test.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
import math
from typing import Iterable


@dataclass(frozen=True)
class Score:
    """Exact moments in a common dyadic scale, plus an attaining subset."""

    _sum: int
    _sum_squares: int
    _denominator: int
    count: int
    indices: tuple[int, ...] = ()

    @property
    def sum(self) -> Fraction:
        return Fraction(self._sum, self._denominator)

    @property
    def sum_squares(self) -> Fraction:
        return Fraction(self._sum_squares, self._denominator ** 2)

    @property
    def value(self) -> float:
        if not self._sum_squares:
            return 0.0
        squared = float(Fraction(self._sum ** 2, self._sum_squares))
        return math.copysign(math.sqrt(squared), self._sum) if self._sum else 0.0

    @property
    def t(self) -> float:
        """Extended paired t display: all-zero is descriptive zero, not a test."""
        if self.count < 2:
            raise ValueError("paired t requires at least two retained observations")
        if not self._sum_squares or not self._sum:
            return 0.0
        variance_numerator = self.count * self._sum_squares - self._sum ** 2
        if variance_numerator == 0:
            return math.copysign(math.inf, self._sum)
        try:
            squared = float(Fraction((self.count - 1) * self._sum ** 2,
                                     variance_numerator))
        except OverflowError:
            squared = math.inf
        return math.copysign(math.sqrt(squared), self._sum)

    def as_dict(self) -> dict:
        t = self.t
        return {"sum": str(self.sum), "sum_squares": str(self.sum_squares),
                "count": self.count, "indices": list(self.indices),
                "score": self.value, "t": t if math.isfinite(t) else str(t),
                "all_zero": self._sum_squares == 0,
                "zero_variance": self.count * self._sum_squares == self._sum ** 2}


def compare_scores(left: Score, right: Score) -> int:
    """Compare signed S/sqrt(Q) exactly, including neutral all-zero scores."""
    ls = (left._sum > 0) - (left._sum < 0)
    rs = (right._sum > 0) - (right._sum < 0)
    if ls != rs:
        return (ls > rs) - (ls < rs)
    if not ls:
        return 0
    a = left._sum ** 2 * right._sum_squares
    b = right._sum ** 2 * left._sum_squares
    return ls * ((a > b) - (a < b))


def compare_t(left: Score, right: Score) -> int:
    """Compare extended paired t exactly, even for different retained counts."""
    if left.count < 2 or right.count < 2:
        raise ValueError("paired t requires at least two retained observations")
    ls = (left._sum > 0) - (left._sum < 0)
    rs = (right._sum > 0) - (right._sum < 0)
    if ls != rs:
        return (ls > rs) - (ls < rs)
    if not ls:
        return 0
    ld = left.count * left._sum_squares - left._sum ** 2
    rd = right.count * right._sum_squares - right._sum ** 2
    if ld == 0 or rd == 0:
        return ls * ((rd == 0) - (ld == 0)) * -1
    a = (left.count-1) * left._sum ** 2 * rd
    b = (right.count-1) * right._sum ** 2 * ld
    return ls * ((a > b) - (a < b))


def compare_abs_scores(left: Score, right: Score) -> int:
    """Exact comparison of |S|/sqrt(Q), with all-zero neutral zero."""
    if not left._sum or not right._sum:
        return ((left._sum != 0) > (right._sum != 0)) - ((left._sum != 0) < (right._sum != 0))
    a = left._sum ** 2 * right._sum_squares
    b = right._sum ** 2 * left._sum_squares
    return (a > b) - (a < b)


def _binary64(values: Iterable[float]) -> tuple[float, ...]:
    out = tuple(float(x) for x in values)
    if any(not math.isfinite(x) for x in out):
        raise ValueError("all effects must be finite binary64 values")
    return out


def _dyadic(values: tuple[float, ...]) -> tuple[tuple[int, ...], int]:
    ratios = tuple(x.as_integer_ratio() for x in values)
    denominator = max((q for _, q in ratios), default=1)
    # All binary64 denominators are powers of two, so divisibility is exact.
    return tuple(p * (denominator // q) for p, q in ratios), denominator


def score(values: Iterable[float]) -> Score:
    vals = _binary64(values)
    if len(vals) < 2:
        raise ValueError("at least two retained observations are required")
    ints, denominator = _dyadic(vals)
    return Score(sum(ints), sum(x*x for x in ints), denominator, len(vals),
                 tuple(range(len(vals))))


@dataclass(frozen=True)
class ExtremaResult:
    minimum: Score
    maximum: Score
    candidate_count: int
    method: str
    total_subsets: int

    def as_dict(self) -> dict:
        return {"minimum": self.minimum.as_dict(), "maximum": self.maximum.as_dict(),
                "candidate_count": self.candidate_count, "method": self.method,
                "total_subsets": self.total_subsets}


class PreparedValues:
    """Sort and encode one gene once; reuse for multiple retained counts."""

    def __init__(self, values: Iterable[float], included: Iterable[float] = ()):
        self.values = _binary64(values)
        self.included = _binary64(included)
        all_ints, self.denominator = _dyadic(self.included + self.values)
        self.fixed_sum = sum(all_ints[:len(self.included)])
        self.fixed_squares = sum(x*x for x in all_ints[:len(self.included)])
        self.ints = all_ints[len(self.included):]
        self.order = tuple(sorted(range(len(self.values)), key=lambda i: (self.values[i], i)))
        self.sorted_ints = tuple(self.ints[i] for i in self.order)
        self.prefix_sum = [0]
        self.prefix_squares = [0]
        for x in self.sorted_ints:
            self.prefix_sum.append(self.prefix_sum[-1] + x)
            self.prefix_squares.append(self.prefix_squares[-1] + x*x)

    def _moment(self, start: int, stop: int) -> tuple[int, int]:
        return (self.prefix_sum[stop] - self.prefix_sum[start],
                self.prefix_squares[stop] - self.prefix_squares[start])

    def extrema(self, choose: int, method: str = "two_family") -> ExtremaResult:
        u = len(self.values)
        if isinstance(choose, bool) or not isinstance(choose, int) or not 0 <= choose <= u:
            raise ValueError("choose must be an integer between zero and optional donor count")
        if len(self.included) + choose < 2:
            raise ValueError("at least two retained observations are required")
        if method not in {"two_family", "exhaustive", "contiguous_only", "ends_only"}:
            raise ValueError("unknown extrema method")

        # Keep descriptors instead of materializing all subsets. Moment queries
        # are O(1) arithmetic operations after sorting; witnesses are made only
        # for the winning descriptors. Absent forced donors, the union has u+2
        # raw candidates and fewer after removing the two duplicate end blocks.
        def candidates():
            if choose in (0, u):
                s, q = self._moment(0, choose)
                yield s, q, ("block", 0)
            elif method == "exhaustive":
                for idx in combinations(range(u), choose):
                    vals = (self.ints[i] for i in idx)
                    s, q = 0, 0
                    for x in vals:
                        s += x
                        q += x*x
                    yield s, q, ("indices", idx)
            else:
                if method in {"two_family", "contiguous_only"}:
                    for start in range(u - choose + 1):
                        s, q = self._moment(start, start + choose)
                        yield s, q, ("block", start)
                if method in {"two_family", "ends_only"}:
                    # j=0 and j=m duplicate the end blocks in the union.
                    js = range(1, choose) if method == "two_family" else range(choose + 1)
                    for j in js:
                        s1, q1 = self._moment(0, j)
                        s2, q2 = self._moment(u - (choose - j), u)
                        yield s1+s2, q1+q2, ("ends", j)

        lo = hi = None
        lo_desc = hi_desc = None
        count = 0
        n = len(self.included) + choose
        for s, q, desc in candidates():
            item = Score(self.fixed_sum + s, self.fixed_squares + q, self.denominator, n)
            if lo is None or compare_scores(item, lo) < 0:
                lo, lo_desc = item, desc
            if hi is None or compare_scores(item, hi) > 0:
                hi, hi_desc = item, desc
            count += 1

        def witness(item: Score, desc: tuple) -> Score:
            kind, v = desc
            if kind == "indices":
                idx = v
            elif kind == "block":
                idx = self.order[v:v+choose]
            else:
                idx = self.order[:v] + (self.order[u-(choose-v):] if choose-v else ())
            return Score(item._sum, item._sum_squares, item._denominator, n, tuple(sorted(idx)))

        assert lo is not None and hi is not None
        return ExtremaResult(witness(lo, lo_desc), witness(hi, hi_desc), count, method,
                             math.comb(u, choose))


def extrema(values: Iterable[float], choose: int, included: Iterable[float] = (),
            method: str = "two_family") -> ExtremaResult:
    """Exact extrema over all choose-sized subsets plus forced included values.

    minimum.indices and maximum.indices refer to positions in `values`, not
    `included`. The latter are retained in every candidate. The two ablation
    methods may miss extrema and must never be used as certificate bounds.
    """
    return PreparedValues(values, included).extrema(choose, method)
