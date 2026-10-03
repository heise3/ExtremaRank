"""Shared-donor top-K audit with exact conditional extrema as safe bounds."""
from __future__ import annotations

from dataclasses import dataclass
from functools import cmp_to_key
import heapq
import math
from typing import Iterable

from .core import PreparedValues, Score, compare_scores, compare_abs_scores


@dataclass(frozen=True)
class AuditResult:
    status: str
    baseline_topk: tuple[str, ...]
    witness_topk: tuple[str, ...] | None
    deleted_indices: tuple[int, ...] | None
    max_deletions: int
    nodes: int
    bound_candidates: int
    retained_counts_completed: tuple[int, ...]
    bound_method: str
    direction: str

    def as_dict(self) -> dict:
        return {"status": self.status, "baseline_topk": list(self.baseline_topk),
                "witness_topk": list(self.witness_topk) if self.witness_topk is not None else None,
                "deleted_indices": list(self.deleted_indices) if self.deleted_indices is not None else None,
                "max_deletions": self.max_deletions, "nodes": self.nodes,
                "bound_candidates": self.bound_candidates,
                "retained_counts_completed": list(self.retained_counts_completed),
                "bound_method": self.bound_method,
                "direction": self.direction,
                "target": "unordered paired-t top-K set; lexical gene-ID ties",
                "arithmetic": "exact moments/comparisons of input binary64 values"}


def audit_topk(effects: Iterable[Iterable[float]], gene_ids: Iterable[str], k: int,
               max_deletions: int, max_nodes: int = 10000,
               bound_method: str = "two_family", witness_trials: int = 8,
               direction: str = "up") -> AuditResult:
    """Certify stability or return a shared deletion witness; otherwise unresolved.

    The matrix is donor x gene. Every deletion re-estimates mean and variance
    on the frozen input effects. A finite node limit can produce UNRESOLVED,
    never a false certificate. `moment_box` is a conservative ablation bound;
    `none` is exhaustive search, useful for checking the effect of the oracle.
    """
    rows = tuple(tuple(float(x) for x in row) for row in effects)
    ids = tuple(str(x) for x in gene_ids)
    d, g = len(rows), len(ids)
    if d < 2 or not g or any(len(row) != g for row in rows):
        raise ValueError("effects must be a nonempty rectangular donor x gene matrix with >=2 donors")
    if len(set(ids)) != g:
        raise ValueError("gene IDs must be unique")
    if isinstance(k, bool) or not isinstance(k, int) or not 1 <= k <= g:
        raise ValueError("k must be an integer in [1, gene count]")
    if isinstance(max_deletions, bool) or not isinstance(max_deletions, int) or not 0 <= max_deletions <= d-2:
        raise ValueError("max_deletions must leave at least two donors")
    if isinstance(max_nodes, bool) or not isinstance(max_nodes, int) or max_nodes < 1:
        raise ValueError("max_nodes must be a positive integer")
    if not isinstance(witness_trials, int) or witness_trials < 0:
        raise ValueError("witness_trials must be a nonnegative integer")
    if bound_method not in {"two_family", "moment_box", "none"}:
        raise ValueError("bound_method must be two_family, moment_box, or none")
    if direction not in {"up", "down", "absolute"}:
        raise ValueError("direction must be up, down, or absolute")
    prepared = [PreparedValues(row[j] for row in rows) for j in range(g)]
    totals = [(sum(p.ints), sum(x*x for x in p.ints)) for p in prepared]

    def target_compare(a: Score, b: Score) -> int:
        if direction == "absolute":
            return compare_abs_scores(a,b)
        c = compare_scores(a,b)
        return -c if direction == "down" else c

    def target_bounds(lo: Score, hi: Score) -> tuple[Score, Score]:
        if direction == "down":
            return hi, lo
        if direction == "absolute":
            zero = Score(0, 0, 1, lo.count)
            # The signed extrema give an exact absolute maximum. Crossing zero
            # only supplies a conservative absolute minimum; an attainable
            # closest-to-zero subset is a different combinatorial problem.
            lower = zero if lo._sum <= 0 <= hi._sum else (lo if compare_abs_scores(lo,hi) <= 0 else hi)
            upper = hi if compare_abs_scores(lo,hi) <= 0 else lo
            return lower, upper
        return lo, hi

    def evaluate(keep: tuple[int, ...]) -> list[Score]:
        keep_set = set(keep)
        removed = tuple(i for i in range(d) if i not in keep_set)
        out = []
        for p, (s, q) in zip(prepared, totals):
            ds = sum(p.ints[i] for i in removed)
            dq = sum(p.ints[i]**2 for i in removed)
            out.append(Score(s-ds, q-dq, p.denominator, len(keep), keep))
        return out

    def compare_gene(a: int, b: int, scores: list[Score]) -> int:
        c = target_compare(scores[a], scores[b])
        return -c if c else ((ids[a] > ids[b]) - (ids[a] < ids[b]))

    def rank(scores: list[Score]) -> tuple[int, ...]:
        key = cmp_to_key(lambda a, b: compare_gene(a, b, scores))
        return tuple(heapq.nsmallest(k, range(g), key=key))

    baseline = rank(evaluate(tuple(range(d))))
    selected = set(baseline)
    winners, outsiders = tuple(baseline), tuple(i for i in range(g) if i not in selected)
    baseline_ids = tuple(ids[i] for i in baseline)
    nodes = bound_candidates = 0
    completed = [d]

    def result(status, witness=None, keep=None):
        return AuditResult(status, baseline_ids,
                           tuple(ids[i] for i in witness) if witness is not None else None,
                           tuple(i for i in range(d) if i not in set(keep)) if keep is not None else None,
                           max_deletions, nodes, bound_candidates, tuple(completed), bound_method, direction)

    if not outsiders or max_deletions == 0:
        return result("CERTIFIED")

    def better(a: int, sa: Score, b: int, sb: Score) -> bool:
        c = target_compare(sa, sb)
        return c > 0 or (c == 0 and ids[a] < ids[b])

    def box(p: PreparedValues, m: int, n: int) -> tuple[Score, Score]:
        # Independent S,Q intervals are wider than the exact coupled moment
        # oracle. Extremal quotients of the enclosing rectangle give safe
        # signed bounds; artificial zero-Q/nonzero-S corners use +/-sqrt(n).
        xs = p.sorted_ints
        s_lo = p.fixed_sum + sum(xs[:m])
        s_hi = p.fixed_sum + (sum(xs[-m:]) if m else 0)
        squares = sorted(x*x for x in xs)
        q_lo = p.fixed_squares + sum(squares[:m])
        q_hi = p.fixed_squares + (sum(squares[-m:]) if m else 0)
        def lower():
            if s_lo < 0 and q_lo == 0:
                return Score(-n, n, 1, n)
            return Score(s_lo, q_lo if s_lo < 0 else q_hi, p.denominator, n)
        def upper():
            if s_hi > 0 and q_lo == 0:
                return Score(n, n, 1, n)
            return Score(s_hi, q_hi if s_hi < 0 else q_lo, p.denominator, n)
        return lower(), upper()

    for removed_count in range(1, max_deletions + 1):
        n = d - removed_count
        stack = [((), tuple(range(d)), n)]
        while stack:
            included, optional, m = stack.pop()
            if nodes >= max_nodes:
                return result("UNRESOLVED")
            nodes += 1
            if m == 0 or m == len(optional):
                keep = tuple(sorted(included + (optional if m else ())))
                actual = rank(evaluate(keep))
                if set(actual) != selected:
                    return result("REFUTED", actual, keep)
                continue
            low, high, witness_keeps = {}, {}, []
            if bound_method != "none":
                for j in range(g):
                    p = PreparedValues((rows[i][j] for i in optional),
                                       (rows[i][j] for i in included))
                    if bound_method == "two_family":
                        env = p.extrema(m)
                        lo, hi = env.minimum, env.maximum
                        bound_candidates += env.candidate_count
                        if len(witness_keeps) < witness_trials:
                            chosen = lo if j in selected else hi
                            witness_keeps.append(tuple(sorted(included + tuple(optional[i] for i in chosen.indices))))
                    else:
                        lo, hi = box(p, m, n)
                        bound_candidates += 4
                    lo, hi = target_bounds(lo,hi)
                    if j in selected:
                        low[j] = lo
                    else:
                        high[j] = hi
                worst = winners[0]
                for j in winners[1:]:
                    if better(worst, low[worst], j, low[j]):
                        worst = j
                best = outsiders[0]
                for j in outsiders[1:]:
                    if better(j, high[j], best, high[best]):
                        best = j
                if better(worst, low[worst], best, high[best]):
                    continue
                # These are proposals, not independent-adversary refutations.
                # Recompute the entire ranking under each one common subset.
                if witness_trials and bound_method == "two_family":
                    targeted = [tuple(sorted(included + tuple(optional[i] for i in s.indices)))
                                for s in (low[worst], high[best]) if len(s.indices) == m]
                    for keep in list(dict.fromkeys(targeted + witness_keeps))[:witness_trials]:
                        actual = rank(evaluate(keep))
                        if set(actual) != selected:
                            return result("REFUTED", actual, keep)
            pivot, rest = optional[0], optional[1:]
            if m <= len(rest):
                stack.append((included, rest, m))
            if m > 0:
                stack.append((included + (pivot,), rest, m-1))
        completed.append(n)
    return result("CERTIFIED")
