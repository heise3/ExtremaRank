"""Compare provided refit ranking tables, with no unobserved-set certificate."""
from __future__ import annotations

import argparse
from html import escape
import json
import math
from pathlib import Path

from . import __version__
from .cli import _write_json, _write_tsv
from .preparation import read_table


def compare_rankings(baseline: dict[str, float], scenarios: dict[str, dict[str, float]],
                     deleted_units: dict[str, list[str]], k: int, direction: str = "up") -> dict:
    """Report only the supplied scenarios; scores must share a fixed universe."""
    if not baseline or any(not isinstance(x, str) or not x for x in baseline):
        raise ValueError("baseline requires nonempty feature IDs")
    if isinstance(k, bool) or not isinstance(k, int) or not 1 <= k <= len(baseline):
        raise ValueError("top-K must be within the feature count")
    if direction not in {"up", "down", "absolute"}:
        raise ValueError("invalid direction")
    if not scenarios or set(scenarios) != set(deleted_units):
        raise ValueError("scenario IDs and deletion manifest must match and be nonempty")
    ids = tuple(baseline)
    def ranking(scores):
        if set(scores) != set(ids):
            raise ValueError("every scenario must retain the complete fixed feature universe")
        converted = {fid: float(value) for fid, value in scores.items()}
        if any(not math.isfinite(x) for x in converted.values()):
            raise ValueError("ranking scores must be finite numbers; no implicit missing-value removal")
        def key(fid):
            value = converted[fid]
            return (value if direction == "down" else -abs(value) if direction == "absolute" else -value, fid)
        return sorted(ids, key=key)
    original = ranking(baseline)
    selected = set(original[:k])
    base_ranks = {fid: i for i, fid in enumerate(original, 1)}
    min_rank, max_rank = dict(base_ranks), dict(base_ranks)
    selected_count = dict.fromkeys(ids, 0)
    rows = []
    for sid, scores in scenarios.items():
        units = deleted_units[sid]
        if not isinstance(sid, str) or not sid or not isinstance(units, list) or not units or any(
                not isinstance(x, str) or not x for x in units) or len(set(units)) != len(units):
            raise ValueError("each scenario needs a nonempty ID and a JSON list of unique nonempty deleted unit IDs")
        order = ranking(scores)
        after = set(order[:k])
        for i, fid in enumerate(order, 1):
            min_rank[fid] = min(min_rank[fid], i)
            max_rank[fid] = max(max_rank[fid], i)
            selected_count[fid] += fid in after
        rows.append({"scenario_id": sid, "deleted_units": units, "topk_changed": after != selected,
            "jaccard": len(after & selected) / len(after | selected),
            "exited_features": [fid for fid in original[:k] if fid not in after],
            "entered_features": [fid for fid in order[:k] if fid not in selected]})
    return {"status": "OBSERVED_CHANGED" if any(row["topk_changed"] for row in rows) else "OBSERVED_STABLE",
        "scope": "only the provided refit tables; no certificate for any unobserved deletion combination",
        "refit_provenance": "caller-supplied scores and deletion manifest; external model fits not verified by this comparator",
        "direction": direction, "top_k": k, "baseline_topk": original[:k], "n_scenarios": len(rows),
        "selection_fraction_meaning": "frequency across supplied scenarios, not a probability",
        "scenarios": rows, "features": [{"feature_id": fid, "baseline_rank": base_ranks[fid],
            "best_observed_rank": min_rank[fid], "worst_observed_rank": max_rank[fid],
            "observed_selection_count": selected_count[fid], "observed_scenarios": len(rows),
            "baseline_topk": fid in selected} for fid in original]}


def main(argv=None):
    parser = argparse.ArgumentParser(prog="extremarank compare", description=__doc__)
    parser.add_argument("baseline", type=Path, help="feature_id,score table")
    parser.add_argument("refits", type=Path, help="scenario_id,feature_id,score long table")
    parser.add_argument("--manifest", type=Path, required=True, help="scenario_id,deleted_units; deleted_units is a JSON array")
    parser.add_argument("--top-k", type=int, default=20)
    parser.add_argument("--direction", choices=("up", "down", "absolute"), default="up")
    parser.add_argument("--method", default="external refit", help="record the external method and statistic")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        bh, br, bsha = read_table(args.baseline)
        rh, rr, rsha = read_table(args.refits)
        mh, mr, msha = read_table(args.manifest)
        if not {"feature_id", "score"} <= set(bh) or not {"scenario_id", "feature_id", "score"} <= set(rh) or not {"scenario_id", "deleted_units"} <= set(mh):
            raise ValueError("required baseline/refit/manifest columns are missing")
        baseline, scenarios, deleted = {}, {}, {}
        for row in br:
            rec = dict(zip(bh, row))
            fid = rec["feature_id"].strip()
            if fid in baseline:
                raise ValueError("duplicate baseline feature ID")
            baseline[fid] = float(rec["score"])
        for row in rr:
            rec = dict(zip(rh, row))
            sid, fid = rec["scenario_id"].strip(), rec["feature_id"].strip()
            scores = scenarios.setdefault(sid, {})
            if fid in scores:
                raise ValueError("duplicate feature ID within a refit scenario")
            scores[fid] = float(rec["score"])
        for row in mr:
            rec = dict(zip(mh, row))
            sid = rec["scenario_id"].strip()
            if sid in deleted:
                raise ValueError("duplicate manifest scenario ID")
            deleted[sid] = json.loads(rec["deleted_units"])
        report = compare_rankings(baseline, scenarios, deleted, args.top_k, args.direction)
        report.update({"software": {"name": "ExtremaRank", "version": __version__}, "method": args.method,
                       "input_sha256": {"baseline": bsha, "refits": rsha, "manifest": msha}})
        args.output.mkdir(parents=True, exist_ok=True)
        _write_json(args.output / "comparison.json", report)
        _write_tsv(args.output / "feature_sensitivity.tsv", tuple(report["features"][0]), report["features"])
        rows = [{key: json.dumps(value, ensure_ascii=False) if isinstance(value, list) else value
                 for key, value in row.items()} for row in report["scenarios"]]
        _write_tsv(args.output / "scenario_changes.tsv", tuple(rows[0]), rows)
        content = '<!doctype html><html lang="en"><meta charset="utf-8"><title>ExtremaRank external refits</title><style>body{font:16px system-ui;max-width:1100px;margin:30px auto;padding:20px}td,th{border-bottom:1px solid #ccc;padding:12px;text-align:left}table{border-collapse:collapse;width:100%}</style><h1>'+escape(report["status"])+ '</h1><p>'+escape(report["method"])+ '</p><p>'+escape(report["scope"])+ '</p><p>'+escape(report["refit_provenance"])+ '</p><table><tr><th>Scenario</th><th>Deleted units</th><th>Jaccard</th><th>Exits</th><th>Entries</th></tr>'
        for row in report["scenarios"]:
            content += '<tr>'+''.join('<td>'+escape(str(value))+'</td>' for value in
                (row["scenario_id"], ', '.join(row["deleted_units"]), row["jaccard"], ', '.join(row["exited_features"]), ', '.join(row["entered_features"])))+'</tr>'
        (args.output / "report.html").write_text(content+'</table></html>', encoding="utf-8")
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(f"{report['status']}: {report['n_scenarios']} supplied refit scenarios; {args.output.resolve()}")
    return 0
