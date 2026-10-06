"""Shared exact membership/direction searches and bounded minimum changes.

Bounds enclose statistics; their endpoints need not share an attaining subset.
Only an evaluated common deletion can refute a property. Smaller deletion
cardinalities must be excluded before a minimum is labelled exact.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations, islice
import hashlib
import json
import math
import struct

from .core import PreparedValues
from .diagnostics import PreparedPaired
from .unpaired import PreparedWelch, WelchScore, compare_welch, rank_welch


def _stat(delta, variance):
    delta, variance = F(delta), F(variance)
    square = delta * delta
    ratio = square / variance if variance else square
    return WelchScore((delta > 0)-(delta < 0), ratio.numerator,
                      ratio.denominator if variance else 0, delta, 0, 0)


def _paired_stat(s):
    return WelchScore((s._sum > 0)-(s._sum < 0), (s.count-1)*s._sum*s._sum,
                      s.count*s._sum_squares-s._sum*s._sum,
                      F(s._sum, s.count*s._denominator), s.count, 0)


def _group_box(ints, included, optional, m):
    fixed = [ints[i] for i in included]
    xs = sorted(ints[i] for i in optional)
    n = len(fixed)+m
    sf, qf = sum(fixed), sum(x*x for x in fixed)
    slo, shi = sf+sum(xs[:m]), sf+(sum(xs[-m:]) if m else 0)
    squares = sorted(x*x for x in xs)
    qlo = qf+sum(squares[:m])
    qhi = qf+(sum(squares[-m:]) if m else 0)
    # Relax forced inclusion: the minimum n-point variance among all allowed
    # values is attained by a sorted window. A relaxed minimum is a safe floor.
    allowed = sorted(fixed+xs)
    s, q = sum(allowed[:n]), sum(x*x for x in allowed[:n])
    window_floor = n*q-s*s
    for start in range(1, len(allowed)-n+1):
        a, b = allowed[start-1], allowed[start+n-1]
        s, q = s-a+b, q-a*a+b*b
        window_floor = min(window_floor, n*q-s*s)
    min_s2 = 0 if slo <= 0 <= shi else min(slo*slo, shi*shi)
    vlo = max(0, window_floor, n*qlo-max(slo*slo, shi*shi))
    vhi = max(0, n*qhi-min_s2)
    return F(slo, n), F(shi, n), F(vlo, n*n*(n-1)), F(vhi, n*n*(n-1))


class SearchProblem:
    def __init__(self, prepared, ids, design, direction, k):
        self.prepared, self.ids = prepared, tuple(ids)
        self.design, self.direction, self.k = design, direction, k
        self.n, self.g = len(prepared.rows), len(self.ids)
        if design not in {"paired", "welch"} or direction not in {"up", "down", "absolute"}:
            raise ValueError("invalid design or ranking direction")
        if self.g != len(prepared.rows[0]) or len(set(self.ids)) != self.g or any(not isinstance(x, str) or not x for x in self.ids):
            raise ValueError("feature IDs must be unique nonempty strings matching matrix columns")
        if isinstance(k, bool) or not isinstance(k, int) or not 1 <= k <= self.g:
            raise ValueError("top-K must be within the feature count")
        self.labels = tuple(int(not x) for x in prepared.is_target) if design == "welch" else (0,)*self.n
        self.counts = tuple(self.labels.count(g) for g in range(2 if design == "welch" else 1))
        self.maximum = self.n-2*len(self.counts)
        scores = self.evaluate(())
        self.baseline_order = rank_welch(scores, self.ids, direction)
        self.baseline_ranks = {j:r for r,j in enumerate(self.baseline_order,1)}
        self.baseline_topk = set(self.baseline_order[:k])
        self.baseline_signs = tuple(s.delta for s in scores)

    def evaluate(self, deleted):
        scores = self.prepared.evaluate(deleted)
        return ([_paired_stat(s) for s in scores] if self.design == "paired" else
                [WelchScore((s.delta > 0)-(s.delta < 0), s.squared_numerator,
                            s.squared_denominator, s.effect, s.target_count, s.reference_count) for s in scores])

    def allocations(self, r):
        if self.design == "paired":
            yield (self.n-r,)
        else:
            t, c = self.counts
            for a in range(max(0, r-c+2), min(r, t-2)+1):
                yield t-a, c-(r-a)

    def deletions(self, r):
        if self.design == "welch":
            for retained in self.allocations(r):
                for a in combinations(self.prepared.target, self.counts[0]-retained[0]):
                    for b in combinations(self.prepared.reference, self.counts[1]-retained[1]):
                        yield tuple(sorted(a+b))
        else:
            yield from combinations(range(self.n), r)

    def canonical(self, included, optional, needed):
        inc, opt, need = list(included), list(optional), list(needed)
        for group in range(len(need)):
            subset = [i for i in opt if self.labels[i] == group]
            if not 0 <= need[group] <= len(subset):
                return None
            if need[group] in {0, len(subset)}:
                if need[group]:
                    inc.extend(subset)
                opt = [i for i in opt if self.labels[i] != group]
                need[group] = 0
        return tuple(sorted(inc)), tuple(opt), tuple(need)

    def bounds(self, included, optional, needed):
        out = []
        if self.design == "paired":
            for j in range(self.g):
                p = PreparedValues((self.prepared.rows[i][j] for i in optional),
                                   (self.prepared.rows[i][j] for i in included))
                env = p.extrema(needed[0])
                out.append((_paired_stat(env.minimum), _paired_stat(env.maximum)))
        else:
            inc = [tuple(i for i in included if self.labels[i] == g) for g in (0, 1)]
            opt = [tuple(i for i in optional if self.labels[i] == g) for g in (0, 1)]
            for ints, _den in self.prepared.encoded:
                a = _group_box(ints, inc[0], opt[0], needed[0])
                b = _group_box(ints, inc[1], opt[1], needed[1])
                dlo, dhi = a[0]-b[1], a[1]-b[0]
                vlo, vhi = a[2]+b[2], a[3]+b[3]
                lo = _stat(dlo, vlo if dlo < 0 else vhi)
                hi = _stat(dhi, vlo if dhi > 0 else vhi)
                out.append((lo, hi))
        return out

    def better(self, j, a, l, b):
        c = compare_welch(a, b, self.direction == "absolute")
        if self.direction == "down":
            c = -c
        return c > 0 or (c == 0 and self.ids[j] < self.ids[l])

    def target_bounds(self, lo, hi):
        if self.direction == "down":
            return hi, lo
        if self.direction == "absolute":
            c = compare_welch(lo, hi, True)
            lower = _stat(0, 1) if lo.delta <= 0 <= hi.delta else (lo if c <= 0 else hi)
            return lower, hi if c <= 0 else lo
        return lo, hi

    def proves(self, query, bounds):
        j, prop = query
        if prop == "direction":
            lo, hi = bounds[j]
            sign = self.baseline_signs[j]
            return lo.delta > 0 if sign > 0 else hi.delta < 0 if sign < 0 else lo.delta == hi.delta == 0
        low, high = self.target_bounds(*bounds[j])
        possible = definite = 0
        for l, interval in enumerate(bounds):
            if l == j:
                continue
            ll, lh = self.target_bounds(*interval)
            possible += self.better(l, lh, j, low)
            definite += self.better(l, ll, j, high)
        return possible < self.k if j in self.baseline_topk else definite >= self.k

    def prover(self, bounds, queries):
        """Index interval endpoints for many queries; avoid G-by-G scans."""
        if sum(prop == "membership" for _,prop in queries) < 32:
            return lambda query: self.proves(query,bounds)
        target=[self.target_bounds(*interval) for interval in bounds]
        lows,highs=[list(x) for x in zip(*target)]
        low_order=rank_welch(lows,self.ids,self.direction)
        high_order=rank_welch(highs,self.ids,self.direction)
        def count_better(order,scores,j,point):
            lo,hi=0,len(order)
            while lo<hi:
                mid=(lo+hi)//2;other=order[mid]
                if self.better(other,scores[other],j,point):lo=mid+1
                else:hi=mid
            return lo
        def proves(query):
            j,prop=query
            if prop=="direction":return self.proves(query,bounds)
            possible=count_better(high_order,highs,j,lows[j])-int(self.better(j,highs[j],j,lows[j]))
            definite=count_better(low_order,lows,j,highs[j])
            return possible<self.k if j in self.baseline_topk else definite>=self.k
        return proves


def audit_feature_robustness(prepared, feature_ids, k, budget, design="paired", direction="up",
                             features=None, max_nodes=10000, max_scenarios=10000,
                             witness_trials=8, use_bounds=True, stop_at_first=False):
    """Certify each queried membership/direction and bound its minimum change.

    `features=None` queries original Top-K features. All comparisons use exact
    input moments. Limits apply globally across features and cardinalities.
    """
    p = SearchProblem(prepared, feature_ids, design, direction, k)
    digest=hashlib.sha256(json.dumps({"shape":[p.n,p.g],"feature_ids":p.ids,"design":design,"group_labels":p.labels},ensure_ascii=False,separators=(',',':')).encode('utf-8')+b'\0')
    for row in prepared.rows:digest.update(struct.pack('<'+'d'*p.g,*row))
    for value, name in ((budget, "budget"), (max_nodes, "max_nodes"), (max_scenarios, "max_scenarios"), (witness_trials, "witness_trials")):
        if isinstance(value, bool) or not isinstance(value, int) or value < (1 if name.startswith("max_") else 0):
            raise ValueError(f"{name} must be an integer within its allowed range")
    if budget > p.maximum:
        raise ValueError("budget must retain at least two observations per native group")
    chosen = tuple(p.ids[j] for j in p.baseline_order[:k]) if features is None else tuple(features)
    positions = {fid: j for j, fid in enumerate(p.ids)}
    if not chosen or len(set(chosen)) != len(chosen) or any(fid not in positions for fid in chosen):
        raise ValueError("queried features must be unique members of the fixed universe")
    queries = tuple((positions[fid], prop) for fid in chosen for prop in ("membership", "direction"))
    state = {q: {"lower": 1, "upper": None, "witness": None} for q in queries}
    nodes = evaluated_count = pruned = 0
    evaluated = set()
    first_witness = None
    stopped = False

    def inspect(deleted, r):
        nonlocal evaluated_count, first_witness
        if deleted in evaluated:
            return
        evaluated.add(deleted)
        evaluated_count += 1
        scores = p.evaluate(deleted)
        order = rank_welch(scores, p.ids, direction, k)
        selected = set(order)
        changed_set = selected != p.baseline_topk
        witness = {"deleted_indices": list(deleted), "topk": [p.ids[j] for j in order]}
        if changed_set and first_witness is None:
            first_witness = witness
        for q, item in state.items():
            j, prop = q
            changed = ((j in selected) != (j in p.baseline_topk)) if prop == "membership" else scores[j].delta != p.baseline_signs[j]
            if changed and (item["upper"] is None or r < item["upper"]):
                item.update(upper=r, witness={**witness, "effect_sign": scores[j].delta})

    for r in range(1, budget+1):
        active = tuple(q for q in queries if state[q]["upper"] is None and (not stop_at_first or q[1] == "membership"))
        if not active:
            break
        for deleted in islice(p.deletions(r), witness_trials):
            if evaluated_count >= max_scenarios:
                break
            inspect(deleted, r)
            if stop_at_first and first_witness:
                stopped = True
                break
        if stopped:
            break
        active = tuple(q for q in active if state[q]["upper"] is None)
        allocations = tuple(p.allocations(r))
        feasible_r = sum(math.prod(math.comb(n, count) for n, count in zip(p.counts, retained)) for retained in allocations)
        # Root bounds first; small residual cardinalities are cheaper to scan
        # directly than repeatedly construct conditional enclosures. This is
        # a cost choice only, never an approximation to the target.
        if feasible_r <= 2000:
            needs_search = set(active) if not use_bounds else set()
            for retained in allocations if use_bounds else ():
                if nodes >= max_nodes:
                    stopped = True
                    needs_search.update(active)
                    break
                nodes += 1
                bounds = p.bounds((), tuple(range(p.n)), retained)
                proves = p.prover(bounds,active)
                remaining = {q for q in active if not proves(q)}
                pruned += len(active)-len(remaining)
                needs_search.update(remaining)
            pending = set()
            if stopped:
                pending = {q for q in needs_search if state[q]["upper"] is None}
            elif needs_search:
                for deleted in p.deletions(r):
                    unresolved = {q for q in needs_search if state[q]["upper"] is None}
                    if not unresolved:
                        break
                    if deleted in evaluated:
                        continue
                    if evaluated_count >= max_scenarios or nodes >= max_nodes:
                        stopped = True
                        pending = unresolved
                        break
                    nodes += 1
                    inspect(deleted, r)
                    if stop_at_first and first_witness:
                        stopped = True
                        pending = {q for q in active if state[q]["upper"] is None}
                        break
            for q in active:
                if state[q]["upper"] is None and q not in pending:
                    state[q]["lower"] = r+1
            if stopped:
                break
            continue
        stack = [((), tuple(range(p.n)), retained, active) for retained in allocations]
        pending = set()
        while stack:
            if nodes >= max_nodes or evaluated_count >= max_scenarios:
                pending = {q for node in stack for q in node[3] if state[q]["upper"] is None}
                stopped = True
                break
            inc, opt, need, candidates = stack.pop()
            candidates = tuple(q for q in candidates if state[q]["upper"] is None)
            if not candidates:
                continue
            canonical = p.canonical(inc, opt, need)
            if canonical is None:
                continue
            inc, opt, need = canonical
            nodes += 1
            if not opt:
                deleted = tuple(i for i in range(p.n) if i not in set(inc))
                inspect(deleted, r)
                if stop_at_first and first_witness:
                    stopped = True
                    pending = {q for q in active if state[q]["upper"] is None}
                    break
                continue
            if use_bounds:
                bounds = p.bounds(inc, opt, need)
                proves = p.prover(bounds,candidates)
                remaining = tuple(q for q in candidates if not proves(q))
                pruned += len(candidates)-len(remaining)
                candidates = remaining
            if not candidates:
                continue
            pivot, rest = opt[0], opt[1:]
            group = p.labels[pivot]
            stack.append((inc, rest, need, candidates))
            if need[group] > 0:
                updated = list(need)
                updated[group] -= 1
                stack.append((inc+(pivot,), rest, tuple(updated), candidates))
        for q in active:
            if state[q]["upper"] is None and q not in pending:
                state[q]["lower"] = r+1
        if stopped:
            break
    rows = []
    for fid in chosen:
        j = positions[fid]
        row = {"feature_id": fid, "baseline_rank": p.baseline_ranks[j],
               "baseline_topk": j in p.baseline_topk, "baseline_effect_sign": p.baseline_signs[j]}
        for prop in ("membership", "direction"):
            item = state[(j, prop)]
            upper, lower = item["upper"], item["lower"]
            status = "REFUTED" if upper is not None else "CERTIFIED" if lower > budget else "UNRESOLVED"
            row[prop] = {"status": status, "minimum_change_lower_bound": lower,
                         "minimum_change_upper_bound": upper,
                         "exact_minimum_change": upper if upper is not None and lower == upper else None,
                         "certified_through": min(budget, lower-1), "witness": item["witness"]}
        rows.append(row)
    memberships = [row["membership"]["status"] for row in rows if row["baseline_topk"]]
    all_baseline_queried = p.baseline_topk <= {positions[fid] for fid in chosen}
    status = ("REFUTED" if first_witness else "CERTIFIED" if all_baseline_queried and all(x == "CERTIFIED" for x in memberships) else "UNRESOLVED")
    topk_lower = min(row["membership"]["minimum_change_lower_bound"] for row in rows if row["baseline_topk"]) if all_baseline_queried else 1
    topk_upper = len(first_witness["deleted_indices"]) if first_witness else None
    return {"status": status, "baseline_topk": [p.ids[j] for j in p.baseline_order[:k]],
            "budget": budget, "direction": direction, "design": design, "features": rows,
            "first_topk_witness": first_witness, "nodes": nodes,
            "topk_minimum_change": {"lower_bound": topk_lower, "upper_bound": topk_upper,
                                    "exact": topk_upper if topk_lower == topk_upper else None},
            "scenarios_checked": evaluated_count, "query_bounds_pruned": pruned,
            "limits": {"max_nodes": max_nodes, "max_scenarios": max_scenarios},
            "input_binary64_sha256":digest.hexdigest(),
            "hash_scope":"UTF-8 compact JSON {shape,feature_ids,design,group_labels} in input order, NUL, then row-major little-endian binary64 values; group labels 0=target/paired, 1=reference",
            "scope": "exact membership and effect-sign properties of frozen finite binary64 native inputs; lower/upper minimum deletion bounds; lexical feature-ID ties",
            "bound_method": "paired conditional two-family" if design == "paired" else "rational group mean/variance enclosures with relaxed minimum-variance windows"}


def write_feature_report(report,unit_ids,output):
    """Common exports for direct paired effects and prepared study inputs."""
    from .cli import _write_json, _write_tsv
    report['unit_ids'] = list(unit_ids)
    if report['first_topk_witness']:
        witness=report['first_topk_witness'];witness['deleted_units']=[unit_ids[i] for i in witness['deleted_indices']]
    flattened=[]
    for row in report['features']:
        item={key:value for key,value in row.items() if key not in {'membership','direction'}}
        for prop in ('membership','direction'):
            witness=row[prop]['witness']
            if witness:witness['deleted_units']=[unit_ids[i] for i in witness['deleted_indices']]
            item.update({prop+'_'+key:json.dumps(value,ensure_ascii=False) if key=='witness' else value for key,value in row[prop].items()})
        flattened.append(item)
    _write_json(output/'feature_robustness.json',report)
    _write_tsv(output/'feature_robustness.tsv',tuple(flattened[0]),flattened)
