"""Run declared R models under whole-donor deletion; retain failed scenarios."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
from html import escape
from itertools import combinations, islice
import json
import math
from pathlib import Path
import re
import shutil
import subprocess

from . import __version__
from .cli import _write_json, _write_tsv
from .external import compare_rankings
from .preparation import read_table

MODELS = ("limma", "limma-voom", "edgeR", "DESeq2")


def run_refits(matrix, metadata, target, reference, model, output, k=20, budget=1,
               direction="up", categorical=(), numeric=(), paired=False,
               min_total_count=10, max_refits=1000, rscript="Rscript",
               min_count=1, min_count_samples=0):
    matrix, metadata, output = Path(matrix), Path(metadata), Path(output)
    if model not in MODELS or target == reference or not target or not reference:
        raise ValueError("declare a supported model and distinct target/reference")
    if any(isinstance(x, bool) or not isinstance(x, int) for x in (k, budget, max_refits, min_total_count, min_count, min_count_samples)) or k < 1 or budget < 1 or max_refits < 1 or min_total_count < 0 or min_count < 1 or min_count_samples < 0:
        raise ValueError("top-K/max-refits positive, budget >=1, min-total-count >=0")
    if direction not in {"up", "down", "absolute"}:
        raise ValueError("invalid ranking direction")
    mh, mr, matrix_sha = read_table(matrix)
    sh, sr, metadata_sha = read_table(metadata)
    if mh[0] != "feature_id" or len(mh) < 2 or not {"sample_id", "group"} <= set(sh):
        raise ValueError("matrix needs feature_id x samples; metadata needs sample_id/group")
    covariates = tuple(categorical)+tuple(numeric)
    if len(set(covariates)) != len(covariates) or any(not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", x) or x in {"group", "sample_id", "donor_id"} or x not in sh for x in covariates):
        raise ValueError("covariates must be distinct declared metadata columns, excluding reserved IDs/group")
    meta = {}
    for row in sr:
        rec = dict(zip(sh, (x.strip() for x in row)))
        sid = rec["sample_id"]
        if not sid or sid in meta or not rec["group"]:
            raise ValueError("metadata IDs/groups must be unique/nonempty")
        rec.setdefault("donor_id", sid)
        if not rec["donor_id"] or any(not rec[x] for x in covariates):
            raise ValueError("missing donor/covariate values are not permitted")
        for name in numeric:
            if not math.isfinite(float(rec[name])):
                raise ValueError("numeric covariates must be finite")
        meta[sid] = rec
    if set(meta) != set(mh[1:]):
        raise ValueError("matrix and metadata sample ID sets must match exactly")
    samples = tuple(s for s in mh[1:] if meta[s]["group"] in {target, reference})
    if any(sum(meta[s]["group"] == g for s in samples) < 2 for g in (target, reference)):
        raise ValueError("at least two samples in each contrast group are required")
    donor_samples = {}
    for sid in samples:
        donor_samples.setdefault(meta[sid]["donor_id"], []).append(sid)
    if paired:
        if any(len(ss) != 2 or {meta[s]["group"] for s in ss} != {target, reference} for ss in donor_samples.values()):
            raise ValueError("paired refits require exactly one sample per group per donor")
        maximum = len(donor_samples)-2
    else:
        if any(len(ss) != 1 for ss in donor_samples.values()):
            raise ValueError("repeated donors require --paired; aggregate technical replicates first")
        maximum = len(samples)-4
    if budget > maximum:
        raise ValueError("deletion budget must retain >=2 samples per group / complete donor pairs")
    count_model = model != "limma"
    source_ids, values, excluded = set(), [], []
    columns = [mh.index(s) for s in samples]
    totals = [0.0]*len(samples)
    for row in mr:
        fid = row[0].strip()
        if not fid or fid in source_ids:
            raise ValueError("feature IDs must be unique/nonempty")
        source_ids.add(fid)
        if count_model:
            from .pseudobulk import _integer
            xs = [_integer(row[i]) for i in columns]
        else:
            xs = [float(row[i]) for i in columns]
        if any(not math.isfinite(x) for x in xs):
            raise ValueError("model inputs must be finite; no implicit imputation")
        if count_model and any(x > 2147483647 for x in xs):
            raise ValueError("count models require raw nonnegative integer counts <=2^31-1, not TPM/log counts")
        totals = [a+b for a, b in zip(totals, xs)]
        if count_model and (sum(xs) < min_total_count or sum(xs) == 0 or sum(x >= min_count for x in xs) < min_count_samples):
            excluded.append(fid)
        else:
            values.append((fid, xs))
    if not values or k > len(values) or (count_model and any(not math.isfinite(x) or x <= 0 for x in totals)):
        raise ValueError("empty filtered universe, invalid top-K, or nonpositive count library")
    donors = tuple(donor_samples)
    def scenarios():
        for r in range(1, budget+1):
            for deleted in combinations(donors, r):
                remaining = [s for s in samples if meta[s]["donor_id"] not in deleted]
                if all(sum(meta[s]["group"] == g for s in remaining) >= 2 for g in (target, reference)):
                    yield deleted
    if paired:
        feasible = sum(math.comb(len(donors), r) for r in range(1, budget+1))
    else:
        n, m = (sum(meta[s]["group"] == g for s in samples) for g in (target, reference))
        feasible = sum(math.comb(n,a)*math.comb(m,r-a) for r in range(1,budget+1)
                       for a in range(max(0,r-m+2),min(r,n-2)+1))
    chosen = list(islice(scenarios(), max_refits))
    names = ("model_input.tsv", "model_metadata.tsv", "model_scenarios.tsv", "baseline.tsv", "refits.tsv",
             "fit_status.tsv", "model_audit.json", "report.html", "R.log", "versions.tsv", "expected_topk.tsv", "model_masks.tsv", "feature_sensitivity.tsv", "nonfinite_scores.tsv")
    if matrix.resolve() in {(output/x).resolve() for x in names} or metadata.resolve() in {(output/x).resolve() for x in names}:
        raise ValueError("output filenames would overwrite model source inputs")
    output.mkdir(parents=True, exist_ok=True)
    _write_tsv(output/"model_input.tsv", ("feature_id", *samples),
               [{"feature_id": fid, **dict(zip(samples, (format(x,".17g") for x in xs)))} for fid, xs in values])
    model_meta = []
    for sid, total in zip(samples, totals):
        model_meta.append({"sample_id": sid, "group": "target" if meta[sid]["group"] == target else "reference",
            "donor_id": meta[sid]["donor_id"], "library_total": format(total,".17g"),
            **{"cv"+str(i): meta[sid][name] for i,name in enumerate(covariates)}})
    _write_tsv(output/"model_metadata.tsv", tuple(model_meta[0]), model_meta)
    manifest_rows = [{"scenario_id": "delete_"+str(i+1), "deleted_units": json.dumps(list(d)),
                      "retained_samples": json.dumps([s for s in samples if meta[s]["donor_id"] not in d])}
                     for i,d in enumerate(chosen)]
    _write_tsv(output/"model_scenarios.tsv", ("scenario_id", "deleted_units", "retained_samples"), manifest_rows)
    # R receives one explicit sample-mask column per scenario; no JSON library
    # or evaluation of user-provided R formula/code is needed.
    masks = [{"scenario_id": row["scenario_id"], **{s: int(meta[s]["donor_id"] not in chosen[i]) for s in samples}}
             for i,row in enumerate(manifest_rows)]
    _write_tsv(output/"model_masks.tsv", ("scenario_id",*samples), masks)
    r_path = shutil.which(rscript)
    if not r_path:
        raise ValueError("Rscript is required only for automatic model refits; native audits remain dependency-free")
    script = Path(__file__).with_name("r")/"refit.R"
    command = [r_path, str(script), str(output.resolve()), model,
               str(len(categorical)), str(len(numeric)), str(int(paired)), str(k)]
    run = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (output/"R.log").write_text(run.stdout, encoding="utf-8")
    contract = {"software": {"name":"ExtremaRank","version":__version__}, "model":model,
        "contrast":target+" - "+reference, "categorical_covariates":list(categorical), "numeric_covariates":list(numeric),
        "paired":paired, "budget":budget,"top_k":k,"direction":direction,
        "feasible_deletion_sets":feasible,"requested_refits":len(chosen),"all_feasible_scenarios_requested":len(chosen)==feasible,
        "fixed_feature_filter":{"min_total_count":min_total_count if count_model else None,
                                "min_count":min_count if count_model else None,"min_count_samples":min_count_samples if count_model else None,
                                "excluded":excluded,"retained":len(values)},
        "input_sha256":{"matrix":matrix_sha,"metadata":metadata_sha,"R_script":hashlib.sha256(script.read_bytes()).hexdigest()},
        "normalization":"Recomputed per retained scenario; edgeR/voom library totals include all source features before fixed filtering.",
        "scope":"Observed external-model refits only; no all-deletion native certificate, FDR guarantee or biological validation.",
        "identity_check":"Source IDs only; raw genotype/sample identity NOT_EVALUATED"}
    contract["excluded_other_group_samples"] = [s for s in mh[1:] if s not in samples]
    contract["ranking_statistic"]={"limma":"moderated t","limma-voom":"moderated t after voom","DESeq2":"Wald statistic",
                                    "edgeR":"sign(logFC) * sqrt(QL F); negative/nonfinite F is not clipped or silently excluded"}[model]
    if run.returncode or not (output/"fit_status.tsv").exists():
        report = {**contract,"status":"NOT_EVALUABLE","reason":"R adapter failed; inspect R.log","R_exit_code":run.returncode}
    else:
        fh, fr, _ = read_table(output/"fit_status.tsv")
        statuses = [dict(zip(fh,row)) for row in fr]
        valid = {row["scenario_id"] for row in statuses if row["status"]=="OK"}
        if "baseline" not in valid:
            report = {**contract,"status":"NOT_EVALUABLE","reason":"Baseline model/ranking is not estimable","fit_status":statuses}
        else:
            bh, br, _ = read_table(output/"baseline.tsv")
            baseline = {dict(zip(bh,row))["feature_id"]:float(dict(zip(bh,row))["score"]) for row in br}
            refits = {}
            rh, rr, _ = read_table(output/"refits.tsv")
            for row in rr:
                rec=dict(zip(rh,row)); refits.setdefault(rec["scenario_id"],{})[rec["feature_id"]]=float(rec["score"])
            deletions = {row["scenario_id"]:json.loads(row["deleted_units"]) for row in manifest_rows if row["scenario_id"] in valid}
            comparison = compare_rankings(baseline,refits,deletions,k,direction) if refits else {"status":"NOT_EVALUABLE","features":[],"scenarios":[]}
            report = {**contract,**comparison,"fit_status":statuses,"valid_refits":len(refits),
                      "failed_refits":sum(row["status"]!="OK" and row["scenario_id"]!="baseline" for row in statuses)}
            report["refit_provenance"]="R adapter executed the declared model per scenario; exported package versions, fit warnings and complete fixed-universe scores are retained."
            # Failures never disappear from the observed coverage statement.
            if report["failed_refits"] and report["status"]=="OBSERVED_STABLE":
                report["status"]="PARTIALLY_EVALUATED"
            for feature in report["features"]:
                fid=feature["feature_id"]
                baseline_member=fid in comparison["baseline_topk"]
                changed=[len(row["deleted_units"]) for row in comparison["scenarios"]
                         if fid in (row["exited_features"] if baseline_member else row["entered_features"])]
                feature["smallest_observed_membership_change"]=min(changed) if changed else None
            if report["features"]:
                _write_tsv(output/"feature_sensitivity.tsv",tuple(report["features"][0]),report["features"])
    _write_json(output/"model_audit.json",report)
    rows=report.get("fit_status",[])
    content='<!doctype html><html lang="en"><meta charset="utf-8"><title>ExtremaRank model refits</title><style>body{font:16px system-ui;max-width:1000px;margin:30px auto;padding:20px}td,th{padding:10px;border-bottom:1px solid #ddd;text-align:left}table{border-collapse:collapse;width:100%}</style><h1>'+escape(report["status"])+'</h1><p>'+escape(model+": "+contract["contrast"])+ '</p><p>'+escape(contract["scope"])+ '</p><p>Requested '+str(len(chosen))+' of '+str(feasible)+' feasible deletion sets. Failed/non-estimable fits are retained below.</p><table><tr><th>Scenario</th><th>Status</th><th>Detail</th></tr>'
    for row in rows:
        content+='<tr>'+''.join('<td>'+escape(row.get(key,""))+'</td>' for key in ("scenario_id","status","detail"))+'</tr>'
    content+='</table><h2>Observed candidate changes</h2><p>Smallest observed deletion counts are witnesses, not certified minimum changes. Non-estimable and unrequested scenarios remain outside this comparison.</p><table><tr><th>Deleted donors</th><th>Jaccard</th><th>Exits</th><th>Entries</th></tr>'
    for row in report.get("scenarios",[]):
        content+='<tr>'+''.join('<td>'+escape(str(value))+'</td>' for value in (', '.join(row['deleted_units']),row['jaccard'],', '.join(row['exited_features']),', '.join(row['entered_features'])))+'</tr>'
    (output/"report.html").write_text(content+'</table>',encoding="utf-8")
    return report


def main(argv=None):
    p=argparse.ArgumentParser(prog="extremarank refit",description="Refit declared R models after shared whole-donor deletion.")
    p.add_argument("matrix",type=Path); p.add_argument("--metadata",type=Path,required=True)
    p.add_argument("--target",required=True); p.add_argument("--reference",required=True)
    p.add_argument("--model",choices=MODELS,required=True); p.add_argument("--output",type=Path,required=True)
    p.add_argument("--top-k",type=int,default=20); p.add_argument("--budget",type=int,default=1)
    p.add_argument("--direction",choices=("up","down","absolute"),default="up")
    p.add_argument("--covariate",action="append",default=[]); p.add_argument("--numeric-covariate",action="append",default=[])
    p.add_argument("--paired",action="store_true"); p.add_argument("--min-total-count",type=int,default=10)
    p.add_argument("--max-refits",type=int,default=1000); p.add_argument("--rscript",default="Rscript")
    p.add_argument("--min-count",type=int,default=1); p.add_argument("--min-count-samples",type=int,default=0)
    a=p.parse_args(argv)
    try:
        report=run_refits(a.matrix,a.metadata,a.target,a.reference,a.model,a.output,a.top_k,a.budget,a.direction,
                          a.covariate,a.numeric_covariate,a.paired,a.min_total_count,a.max_refits,a.rscript,a.min_count,a.min_count_samples)
    except (OSError,ValueError) as e:
        p.error(str(e))
    print(report["status"]+": "+a.model+"; "+str(a.output.resolve()))
    return 1 if report["status"]=="NOT_EVALUABLE" else 0
