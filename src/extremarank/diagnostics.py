"""Fully inspected delete-one diagnostics, independent of multi-deletion proofs."""
from __future__ import annotations

from functools import cmp_to_key
from html import escape
import json

from .core import PreparedValues, Score, compare_t, compare_abs_scores
from .unpaired import rank_welch


class PreparedPaired:
    def __init__(self, rows):
        rows = tuple(tuple(float(x) for x in row) for row in rows)
        if len(rows) < 2 or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
            raise ValueError("paired values require a rectangular matrix with >=2 units and >=1 feature")
        self.rows = rows
        self.prepared = [PreparedValues(row[j] for row in rows) for j in range(len(rows[0]))]
        self.totals = [(sum(p.ints), sum(x*x for x in p.ints)) for p in self.prepared]

    def evaluate(self, deleted=()):
        n = len(self.rows) - len(deleted)
        return [Score(s - sum(p.ints[i] for i in deleted),
                      q - sum(p.ints[i]**2 for i in deleted), p.denominator, n)
                for p, (s, q) in zip(self.prepared, self.totals)]


def rank_paired(scores, ids, direction):
    def compare(a, b):
        c = compare_abs_scores(scores[a], scores[b]) if direction == "absolute" else compare_t(scores[a], scores[b])
        if direction == "down":
            c = -c
        return -c if c else (ids[a] > ids[b]) - (ids[a] < ids[b])
    return tuple(sorted(range(len(ids)), key=cmp_to_key(compare)))


def delete_one_diagnostics(prepared, ids, unit_ids, k, direction, design, metadata=()):
    ranker = rank_welch if design == "welch" else rank_paired
    baseline_scores = prepared.evaluate()
    baseline = ranker(baseline_scores, ids, direction)
    selected = set(baseline[:k])
    base_rank = {j: r for r, j in enumerate(baseline, 1)}
    best_rank, worst_rank = dict(base_rank), dict(base_rank)
    count_topk = {j: 0 for j in baseline}
    sign_changed = {j: 0 for j in baseline}
    def sign(score):
        value = score.delta if design == "welch" else score._sum
        return (value > 0) - (value < 0)
    scenarios, skipped = [], []
    for i, uid in enumerate(unit_ids):
        valid = (len(prepared.target if prepared.is_target[i] else prepared.reference) > 2
                 if design == "welch" else len(unit_ids) > 2)
        if not valid:
            skipped.append(uid)
            continue
        scores = prepared.evaluate((i,))
        order = ranker(scores, ids, direction)
        after = set(order[:k])
        for r, j in enumerate(order, 1):
            best_rank[j] = min(best_rank[j], r)
            worst_rank[j] = max(worst_rank[j], r)
            count_topk[j] += j in after
            sign_changed[j] += sign(scores[j]) != sign(baseline_scores[j])
        scenarios.append({"deleted_unit": uid, "deleted_index": i,
            "metadata": metadata[i] if metadata else {},
            "topk_changed": after != selected,
            "jaccard": len(after & selected) / len(after | selected),
            "exited_features": [ids[j] for j in baseline[:k] if j not in after],
            "entered_features": [ids[j] for j in order[:k] if j not in selected]})
    features = [{"feature_id": ids[j], "baseline_rank": base_rank[j], "baseline_t": baseline_scores[j].t,
        "best_observed_rank": best_rank[j], "worst_observed_rank": worst_rank[j],
        "loo_selected_count": count_topk[j], "loo_scenarios": len(scenarios),
        "loo_selection_fraction": count_topk[j] / len(scenarios) if scenarios else None,
        "loo_sign_changes": sign_changed[j], "baseline_topk": j in selected} for j in baseline]
    return {"scope": "baseline plus every feasible one-unit deletion; no multi-deletion rank bounds",
        "selection_fraction_meaning": "fraction of inspected deletion scenarios, not a probability",
        "n_scenarios": len(scenarios), "n_changed": sum(row["topk_changed"] for row in scenarios),
        "skipped_units": skipped, "scenarios": scenarios, "features": features}


def html_report(report, diagnostics, granular=None):
    def table(headers, rows):
        return "<table><thead><tr>" + "".join("<th>"+escape(h)+"</th>" for h in headers) + "</tr></thead><tbody>" + "".join(
            "<tr>"+"".join("<td>"+escape(str(cell))+"</td>" for cell in row)+"</tr>" for row in rows) + "</tbody></table>"
    prep = report.get("preparation", {})
    audit = report["topk_audit"]
    status = escape(report["status"])
    sample_rows = [(s["deleted_unit"], s["metadata"].get("group", "paired"),
                   round(s["jaccard"], 4), ", ".join(s["exited_features"]), ", ".join(s["entered_features"]))
                  for s in diagnostics["scenarios"]]
    features = diagnostics["features"]
    # Show initial candidates and the most rank-sensitive outsiders, with scope
    # stated explicitly. Full feature diagnostics are retained in TSV/JSON.
    displayed = [f for f in features if f["baseline_topk"]]
    extras = sorted((f for f in features if not f["baseline_topk"]),
                    key=lambda f: (-(f["worst_observed_rank"]-f["best_observed_rank"]), f["feature_id"]))[:20]
    feature_rows = [(f["feature_id"], f["baseline_rank"], f["best_observed_rank"], f["worst_observed_rank"],
                     f["loo_selected_count"], f["loo_scenarios"], f["loo_sign_changes"]) for f in displayed+extras]
    witness = audit.get("deleted_donors") or audit.get("deleted_units") or []
    description = ("Every allowed deletion preserves the initial unordered top-K set." if status == "CERTIFIED" else
                   "The listed common deletion changes the initial top-K set." if status == "REFUTED" else
                   "The declared search budget was reached; stability remains unresolved." if status == "UNRESOLVED" else
                   "The shared multi-deletion search was not run.")
    extra = ""
    if granular:
        def radius(item):
            if item["exact_minimum_change"] is not None:
                return str(item["exact_minimum_change"])+" (exact)"
            lo, hi = item["minimum_change_lower_bound"], item["minimum_change_upper_bound"]
            return str(lo)+".."+str(hi)+" (bounded)" if hi is not None else ">="+str(lo)+"; no witness within inspected budget"
        rows = [(row["feature_id"], row["baseline_rank"], row["membership"]["status"],
                 row["membership"]["certified_through"], radius(row["membership"]),
                 row["direction"]["status"], radius(row["direction"])) for row in granular["features"][:200]]
        extra = '<div class="card"><h2>Individual candidate guarantees and minimum changes</h2><p>Membership means preserving each candidate\'s baseline inside/outside Top-K status. Direction means preserving the positive/zero/negative mean-effect sign. Exact minima exclude every smaller feasible deletion; incomplete searches retain bounds. All queried features and common deletion witnesses are in feature_robustness.json.</p><div class="scroll">'+table(
            ["Feature", "Baseline rank", "Membership", "Certified through", "Minimum membership change", "Direction", "Minimum direction change"], rows)+'</div></div>'
    return """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>ExtremaRank study audit</title><style>body{font-family:system-ui,sans-serif;max-width:1200px;margin:32px auto;padding:0 20px;color:#17283a;background:#f8fafc}h1,h2{color:#163f66}.card{background:white;border:1px solid #d5dfeb;border-radius:10px;padding:20px;margin:20px 0}table{border-collapse:collapse;width:100%;font-size:14px}th,td{border-bottom:1px solid #d5dfeb;text-align:left;padding:10px;overflow-wrap:anywhere}th{background:#eef4fa}code{overflow-wrap:anywhere}.scroll{overflow:auto}small{color:#49637c}</style><h1>ExtremaRank study audit</h1>""" + (
        '<div class="card"><h2>'+status+'</h2><p>'+escape(description)+'</p><p>Design: '+escape(prep.get("design", "paired"))+
        ' · Contrast: '+escape(prep.get("contrast", "prepared effects"))+'</p><p>Deletion witness: '+escape(", ".join(witness) or "none")+
        '</p><p>Target: '+escape(report.get("target", audit.get("target", "")))+'</p><small>Exactness concerns supplied binary64 values and this declared ranking. Preprocessing is frozen; no p values, FDR, biological truth or confounder adjustment is established.</small></div>'+
        '<div class="card"><h2>Every feasible one-unit deletion</h2><p>'+str(diagnostics["n_changed"])+ ' of '+str(diagnostics["n_scenarios"])+
        ' inspected deletions changed the set. Skipped units (retained count too small): '+escape(", ".join(diagnostics["skipped_units"]) or "none")+
        '.</p><div class="scroll">'+table(["Deleted unit", "Group", "Jaccard", "Exits", "Entries"], sample_rows)+'</div></div>'+
        '<div class="card"><h2>Candidate rank sensitivity</h2><p>Ranges include baseline and feasible delete-one scenarios only. Selection counts are scenario frequencies, not probabilities. All features are available in diagnostics.tsv.</p><div class="scroll">'+
        table(["Feature", "Baseline", "Best inspected rank", "Worst inspected rank", "LOO selected", "LOO scenarios", "Sign changes"], feature_rows)+
        '</div></div>'+extra+'<p>ExtremaRank '+escape(report.get("software", {}).get("version", ""))+' · Offline report; no external scripts or network requests.</p></html>')
