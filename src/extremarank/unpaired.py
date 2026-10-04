"""Exact Welch rankings under shared deletions of independent samples.

This is bounded exhaustive search, not an extension of the paired two-family
theorem. No p values, Welch degrees of freedom, or covariates are estimated.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from fractions import Fraction
from functools import cmp_to_key
import heapq
from itertools import combinations
import math
from typing import Iterable

from .core import _binary64, _dyadic


@dataclass(frozen=True)
class WelchScore:
    delta: int
    squared_numerator: int
    squared_denominator: int
    effect: Fraction
    target_count: int
    reference_count: int

    @property
    def t(self) -> float:
        if not self.delta:
            return 0.0
        if not self.squared_denominator:
            return math.copysign(math.inf, -1.0 if self.delta < 0 else 1.0)
        try:
            squared = float(Fraction(self.squared_numerator, self.squared_denominator))
        except OverflowError:
            squared = math.inf
        return math.copysign(math.sqrt(squared), -1.0 if self.delta < 0 else 1.0)


def compare_welch(a: WelchScore, b: WelchScore, absolute: bool = False) -> int:
    """Order exact extended Welch statistics; floats are display values only."""
    sa = (a.delta > 0) - (a.delta < 0)
    sb = (b.delta > 0) - (b.delta < 0)
    if absolute:
        sa, sb = abs(sa), abs(sb)
    if sa != sb:
        return (sa > sb) - (sa < sb)
    if not sa:
        return 0
    if not a.squared_denominator or not b.squared_denominator:
        c = (a.squared_denominator == 0) - (b.squared_denominator == 0)
    else:
        x = a.squared_numerator * b.squared_denominator
        y = b.squared_numerator * a.squared_denominator
        c = (x > y) - (x < y)
    return sa * c


class PreparedWelch:
    """Encode a fixed sample x feature matrix once and subtract deleted moments."""

    def __init__(self, values: Iterable[Iterable[float]], groups: Iterable[str],
                 target: str, reference: str):
        self.rows = tuple(_binary64(row) for row in values)
        self.groups = tuple(groups)
        if not target or not reference or target == reference:
            raise ValueError("target and reference must be distinct nonempty group names")
        if not self.rows or not self.rows[0] or any(len(r) != len(self.rows[0]) for r in self.rows):
            raise ValueError("values must be a nonempty rectangular sample x feature matrix")
        if len(self.groups) != len(self.rows) or any(x not in {target, reference} for x in self.groups):
            raise ValueError("one target/reference group is required for each sample")
        self.target = tuple(i for i, g in enumerate(self.groups) if g == target)
        self.reference = tuple(i for i, g in enumerate(self.groups) if g == reference)
        if min(len(self.target), len(self.reference)) < 2:
            raise ValueError("Welch ranking requires at least two independent samples in each group")
        self.is_target = tuple(g == target for g in self.groups)
        self.encoded = []
        self.totals = []
        for j in range(len(self.rows[0])):
            ints, denominator = _dyadic(tuple(r[j] for r in self.rows))
            self.encoded.append((ints, denominator))
            self.totals.append((sum(ints[i] for i in self.target),
                                sum(ints[i] ** 2 for i in self.target),
                                sum(ints[i] for i in self.reference),
                                sum(ints[i] ** 2 for i in self.reference)))

    def evaluate(self, deleted: Iterable[int] = ()) -> list[WelchScore]:
        removed = tuple(deleted)
        if len(set(removed)) != len(removed) or any(isinstance(i, bool) or not isinstance(i, int)
                or not 0 <= i < len(self.rows) for i in removed):
            raise ValueError("deleted indices must be unique valid sample positions")
        dt = tuple(i for i in removed if self.is_target[i])
        dr = tuple(i for i in removed if not self.is_target[i])
        n, m = len(self.target) - len(dt), len(self.reference) - len(dr)
        if min(n, m) < 2:
            raise ValueError("each retained group must contain at least two samples")
        out = []
        for (ints, denominator), (st, qt, sr, qr) in zip(self.encoded, self.totals):
            st -= sum(ints[i] for i in dt)
            qt -= sum(ints[i] ** 2 for i in dt)
            sr -= sum(ints[i] for i in dr)
            qr -= sum(ints[i] ** 2 for i in dr)
            delta = st * m - sr * n
            vtarget, vreference = n * qt - st * st, m * qr - sr * sr
            squared_num = delta * delta * (n - 1) * (m - 1)
            squared_den = vtarget * m * m * (m - 1) + vreference * n * n * (n - 1)
            out.append(WelchScore(delta, squared_num, squared_den,
                       Fraction(delta, n * m * denominator), n, m))
        return out

    def deletions(self, budget: int):
        """Enumerate only shared deletions that retain >=2 per group, by size."""
        for r in range(1, budget + 1):
            for a in range(max(0, r - len(self.reference) + 2), min(r, len(self.target) - 2) + 1):
                for dt in combinations(self.target, a):
                    for dr in combinations(self.reference, r - a):
                        yield tuple(sorted(dt + dr))

    def deletion_count(self, budget: int) -> int:
        return sum(math.comb(len(self.target), a) * math.comb(len(self.reference), r-a)
                   for r in range(1, budget + 1)
                   for a in range(max(0, r-len(self.reference)+2), min(r, len(self.target)-2)+1))


def rank_welch(scores: list[WelchScore], ids: tuple[str, ...], direction: str,
               k: int | None = None) -> tuple[int, ...]:
    if direction not in {"up", "down", "absolute"}:
        raise ValueError("direction must be up, down, or absolute")
    def compare(a, b):
        c = compare_welch(scores[a], scores[b], direction == "absolute")
        if direction == "down":
            c = -c
        return -c if c else (ids[a] > ids[b]) - (ids[a] < ids[b])
    key = cmp_to_key(compare)
    return tuple(sorted(range(len(ids)), key=key) if k is None else
                 heapq.nsmallest(k, range(len(ids)), key=key))


@dataclass(frozen=True)
class WelchAuditResult:
    status: str
    baseline_topk: tuple[str, ...]
    witness_topk: tuple[str, ...] | None
    deleted_indices: tuple[int, ...] | None
    max_deletions: int
    subsets_checked: int
    feasible_deletion_sets: int
    max_subsets: int
    direction: str
    search_method: str = "bounded exhaustive shared-sample enumeration"
    arithmetic: str = "exact Welch comparisons of finite binary64 input values"
    target: str = "unordered Welch-statistic top-K; >=2 retained per group; lexical feature-ID ties"

    def as_dict(self) -> dict:
        return asdict(self)


def audit_welch(values: Iterable[Iterable[float]], groups: Iterable[str], feature_ids: Iterable[str],
                target: str, reference: str, k: int, max_deletions: int,
                max_subsets: int = 10000, direction: str = "up") -> WelchAuditResult:
    prepared = PreparedWelch(values, groups, target, reference)
    return audit_prepared_welch(prepared, feature_ids, k, max_deletions, max_subsets, direction)


def audit_prepared_welch(prepared: PreparedWelch, feature_ids: Iterable[str], k: int,
                         max_deletions: int, max_subsets: int = 10000,
                         direction: str = "up") -> WelchAuditResult:
    ids = tuple(feature_ids)
    if len(ids) != len(prepared.rows[0]) or any(not isinstance(x, str) or not x for x in ids) or len(set(ids)) != len(ids):
        raise ValueError("feature IDs must be nonempty unique strings matching the matrix")
    if isinstance(k, bool) or not isinstance(k, int) or not 1 <= k <= len(ids):
        raise ValueError("k must be an integer within the feature count")
    if isinstance(max_deletions, bool) or not isinstance(max_deletions, int) or not 0 <= max_deletions <= len(prepared.rows)-4:
        raise ValueError("budget must be between 0 and total sample count minus four")
    if isinstance(max_subsets, bool) or not isinstance(max_subsets, int) or max_subsets < 1:
        raise ValueError("max_subsets must be a positive integer")
    baseline = rank_welch(prepared.evaluate(), ids, direction, k)
    selected = set(baseline)
    outsiders = tuple(i for i in range(len(ids)) if i not in selected)
    count, feasible = 0, prepared.deletion_count(max_deletions)
    def result(status, deleted=None, witness=None):
        return WelchAuditResult(status, tuple(ids[i] for i in baseline),
            tuple(ids[i] for i in witness) if witness is not None else None,
            deleted, max_deletions, count, feasible, max_subsets, direction)
    if not outsiders or not max_deletions:
        return result("CERTIFIED")
    def better(a, b, scores):
        c = compare_welch(scores[a], scores[b], direction == "absolute")
        if direction == "down":
            c = -c
        return c > 0 or (c == 0 and ids[a] < ids[b])
    for deleted in prepared.deletions(max_deletions):
        if count >= max_subsets:
            return result("UNRESOLVED")
        scores = prepared.evaluate(deleted)
        count += 1
        worst, best = baseline[0], outsiders[0]
        for j in baseline[1:]:
            if better(worst, j, scores):
                worst = j
        for j in outsiders[1:]:
            if better(j, best, scores):
                best = j
        if better(best, worst, scores):
            return result("REFUTED", deleted, rank_welch(scores, ids, direction, k))
    return result("CERTIFIED")
