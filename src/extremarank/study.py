"""Study preparation, paired/Welch audits, and an offline diagnostic report."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import time

from . import __version__
from .cli import _write_json, _write_tsv, _display
from .diagnostics import PreparedPaired, delete_one_diagnostics, html_report, rank_paired
from .preparation import prepare_study
from .unpaired import PreparedWelch, audit_prepared_welch, rank_welch


def main(argv=None, prepare_only=False, feature_mode=False):
    parser = argparse.ArgumentParser(prog="extremarank prepare" if prepare_only else "extremarank study",
        description="Prepare a declared contrast from a matrix and strictly aligned sample metadata.")
    parser.add_argument("matrix", type=Path, help="feature_id x samples (default), or sample_id x features")
    parser.add_argument("--metadata", type=Path, required=True, help="CSV/TSV with sample_id, group, and donor_id for paired designs")
    parser.add_argument("--target", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--design", choices=("paired", "welch"), default="paired")
    parser.add_argument("--orientation", choices=("features", "samples"), default="features")
    parser.add_argument("--transform", choices=("none", "log2", "cpm-log2"), default="none")
    parser.add_argument("--missing", choices=("error", "drop-features"), default="error")
    parser.add_argument("--incomplete-pairs", choices=("error", "drop"), default="error")
    parser.add_argument("--min-cpm", type=float, default=0.0)
    parser.add_argument("--pseudocount", type=float, default=1.0)
    parser.add_argument("--output", type=Path, required=True)
    if not prepare_only:
        parser.add_argument("--budget", type=int, default=2)
        parser.add_argument("--top-k", type=int, default=None)
        parser.add_argument("--direction", choices=("up", "down", "absolute"), default="up")
        parser.add_argument("--max-nodes", type=int, default=10000, help="branch-and-bound node limit")
        parser.add_argument("--max-subsets", type=int, default=10000, help="Welch nonbaseline deletion-set limit")
        parser.add_argument("--welch-search", choices=("branch-bound", "enumeration"), default="branch-bound")
        parser.add_argument("--feature-audit", action="store_true", default=feature_mode,
                            help="audit individual membership/direction and minimum-change bounds")
        parser.add_argument("--features", choices=("baseline", "all"), default="baseline")
        parser.add_argument("--feature-id", action="append", help="query specified feature IDs instead")
        parser.add_argument("--no-diagnostics", action="store_true", help="omit delete-one diagnostics and HTML")
    args = parser.parse_args(argv)
    if not prepare_only and (args.feature_id or args.features == "all"):
        args.feature_audit = True
    started = time.perf_counter()
    try:
        study = prepare_study(args.matrix, args.metadata, args.target, args.reference,
            args.design, args.orientation, args.transform, args.missing, args.incomplete_pairs,
            args.min_cpm, args.pseudocount)
        n, g = len(study.unit_ids), len(study.feature_ids)
        if not prepare_only:
            maximum = n-2 if args.design == "paired" else n-4
            if not 0 <= args.budget <= maximum:
                raise ValueError(f"budget must be between 0 and {maximum}; retained groups/pairs require >=2 observations")
            k = min(20, g) if args.top_k is None else args.top_k
            if not 1 <= k <= g or args.max_nodes < 1 or args.max_subsets < 1:
                raise ValueError("invalid top-K or search limit")
        reserved = {args.output / name for name in ("effects.csv", "values.csv", "groups.tsv", "preparation.json",
            "baseline.tsv", "summary.tsv", "envelopes.tsv", "audit.json", "diagnostics.tsv", "diagnostics.json", "sample_influence.tsv", "report.html", "feature_robustness.json", "feature_robustness.tsv")}
        if args.matrix.resolve() in {p.resolve() for p in reserved} or args.metadata.resolve() in {p.resolve() for p in reserved}:
            raise ValueError("output filenames would overwrite a source matrix or metadata; choose a separate directory")
        args.output.mkdir(parents=True, exist_ok=True)
        prepared_path = args.output / ("effects.csv" if args.design == "paired" else "values.csv")
        with prepared_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["donor_id" if args.design == "paired" else "sample_id", *study.feature_ids])
            for uid, row in zip(study.unit_ids, study.values):
                writer.writerow([uid, *(format(x, ".17g") for x in row)])
        if args.design == "welch":
            _write_tsv(args.output / "groups.tsv", ("sample_id", "group", "donor_id"),
                       [{"sample_id": uid, "group": group, "donor_id": uid} for uid, group in zip(study.unit_ids, study.groups)])
        study.manifest["prepared_path"] = prepared_path.name
        _write_json(args.output / "preparation.json", study.manifest)
        if prepare_only:
            print(f"PREPARED: {args.design}, {n} independent units, {g} features; {prepared_path.resolve()}")
            return 0
        if args.design == "paired":
            from .cli import main as paired_main
            paired_main([str(prepared_path), "--budget", str(args.budget), "--top-k", str(k),
                         "--direction", args.direction, "--max-nodes", str(args.max_nodes), "--output", str(args.output)])
            report = json.loads((args.output / "audit.json").read_text())
            prepared = PreparedPaired(study.values)
        else:
            prepared = PreparedWelch(study.values, study.groups, args.target, args.reference)
            scores = prepared.evaluate()
            order = rank_welch(scores, study.feature_ids, args.direction)
            audit = audit_prepared_welch(prepared, study.feature_ids, k, args.budget,
                                          args.max_subsets, args.direction, args.welch_search, args.max_nodes).as_dict()
            audit["deleted_units"] = ([study.unit_ids[i] for i in audit["deleted_indices"]]
                                      if audit["deleted_indices"] is not None else None)
            baseline = [{"feature_id": study.feature_ids[j], "baseline_rank": r,
                         "baseline_t": _display(scores[j].t), "target_minus_reference": _display(scores[j].effect),
                         "baseline_topk": r <= k} for r, j in enumerate(order, 1)]
            _write_tsv(args.output / "baseline.tsv", tuple(baseline[0]), baseline)
            report = {"schema_version": 2, "software": {"name": "ExtremaRank", "version": __version__},
                "status": audit["status"], "target": audit["target"], "arithmetic": audit["arithmetic"],
                "parameters": {"budget": args.budget, "top_k": k, "direction": args.direction,
                               "max_subsets": args.max_subsets},
                "input": {"n_units": n, "n_features": g, "unit_ids": study.unit_ids,
                          "target_count": len(prepared.target), "reference_count": len(prepared.reference)},
                "score_conventions": {"equal_means": "t=0, including zero standard error",
                    "unequal_constant_groups": "signed infinity; descriptive extension, no inferential p value"},
                "topk_audit": audit, "outputs": {"baseline": "baseline.tsv"}}
        report["preparation"] = study.manifest
        report["preparation"]["command_inputs"] = {"matrix": str(args.matrix), "metadata": str(args.metadata)}
        granular = None
        if args.feature_audit:
            from .robustness import audit_feature_robustness, write_feature_report
            features = args.feature_id or (study.feature_ids if args.features == "all" else None)
            granular = audit_feature_robustness(prepared, study.feature_ids, k, args.budget,
                args.design, args.direction, features, args.max_nodes, args.max_subsets)
            write_feature_report(granular,study.unit_ids,args.output)
            report["feature_robustness"] = {key: granular[key] for key in ("status", "nodes", "scenarios_checked", "query_bounds_pruned", "scope")}
            report["outputs"].update(feature_robustness="feature_robustness.tsv", minimum_changes="feature_robustness.json")
        if not args.no_diagnostics:
            diagnostics = delete_one_diagnostics(prepared, study.feature_ids, study.unit_ids, k,
                                                  args.direction, args.design, study.unit_metadata)
            _write_json(args.output / "diagnostics.json", diagnostics)
            _write_tsv(args.output / "diagnostics.tsv", tuple(diagnostics["features"][0]),
                       [{key: _display(value) if key == "baseline_t" else value for key, value in row.items()}
                        for row in diagnostics["features"]])
            _write_tsv(args.output / "sample_influence.tsv",
                       ("deleted_unit", "group", "topk_changed", "jaccard", "exited_features", "entered_features"),
                       [{"deleted_unit": row["deleted_unit"], "group": row["metadata"].get("group", "paired"),
                         "topk_changed": row["topk_changed"], "jaccard": row["jaccard"],
                         "exited_features": json.dumps(row["exited_features"], ensure_ascii=False),
                         "entered_features": json.dumps(row["entered_features"], ensure_ascii=False)}
                        for row in diagnostics["scenarios"]])
            report["delete_one_summary"] = {key: diagnostics[key] for key in
                                             ("scope", "n_scenarios", "n_changed", "skipped_units")}
            report["outputs"].update({"diagnostics": "diagnostics.tsv", "sample_influence": "sample_influence.tsv", "report": "report.html"})
            (args.output / "report.html").write_text(html_report(report, diagnostics, granular), encoding="utf-8")
        report.setdefault("timing_seconds", {})["study_total_before_writing"] = time.perf_counter()-started
        _write_json(args.output / "audit.json", report)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(f"{report['status']}: {args.design}, {n} units, {g} features; {args.output.resolve()}")
    return 0
