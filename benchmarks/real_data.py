#!/usr/bin/env python3
"""Prepare paired public RNA inputs and benchmark exact retained-subset extrema.

Default contract: budgets 0..4; all eligible genes for two_family; the first
100 eligible genes in original GeneID row order for exhaustive comparison and
family ablations. All normalization and filtering precede donor deletion.
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction
from functools import cmp_to_key
import gzip
import hashlib
import importlib.util
import io
import json
import math
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter

import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATASETS = ("GSE87290", "GSE50760")
PREDECLARED_MAX_DELETIONS = 4
PREDECLARED_MATCHED_GENES = 100


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def source_verifier():
    spec = importlib.util.spec_from_file_location("extremarank_source_fetch", REPO / "data" / "fetch_sources.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module.verify_sources


def soft_samples(path: Path) -> dict[str, dict]:
    records = {}
    for block in path.read_text().split("^SAMPLE = ")[1:]:
        lines = block.splitlines()
        sample = lines[0].strip()
        record = {"sample_id": sample, "characteristics": {}}
        for line in lines[1:]:
            if line.startswith("!Sample_title = "):
                record["title"] = line.split(" = ", 1)[1]
            elif line.startswith("!Sample_source_name_ch1 = "):
                record["source_name"] = line.split(" = ", 1)[1]
            elif line.startswith("!Sample_characteristics_ch1 = "):
                key, value = line.split(" = ", 1)[1].split(":", 1)
                record["characteristics"][key.strip().lower()] = value.strip()
        if sample in records:
            raise ValueError(f"Duplicate GEO sample {sample}")
        records[sample] = record
    return records


def read_counts(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    with gzip.open(path, "rt", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        if header[0] != "GeneID" or len(set(header[1:])) != len(header) - 1:
            raise ValueError("Unexpected count header or duplicate samples")
        rows, genes = [], []
        for row in reader:
            if len(row) != len(header):
                raise ValueError("Ragged count matrix")
            genes.append(row[0])
            rows.append([int(value) for value in row[1:]])
    counts = np.asarray(rows, dtype=np.int64)
    if len(set(genes)) != len(genes) or np.any(counts < 0):
        raise ValueError("Duplicate GeneID or negative counts")
    return genes, header[1:], counts


def pair_samples(accession: str, samples: dict, actual: list[str]) -> tuple[list[dict], list[dict]]:
    groups = {}
    for sample, metadata in samples.items():
        if accession == "GSE87290":
            match = re.fullmatch(r"PBMC_replicate_(\d+) (baseline|intravenous LPS) (high|low) inflammation", metadata["title"])
            if not match:
                raise ValueError(f"Unexpected source title: {metadata['title']}")
            donor, condition = "PBMC_replicate_" + match[1], match[2]
            if metadata["characteristics"].get("treatment") != condition:
                raise ValueError("Title/treatment disagreement")
        else:
            match = re.search(r"AMC_(\d+)-([123])$", metadata["title"])
            if not match:
                raise ValueError(f"Unexpected source title: {metadata['title']}")
            donor = "AMC_" + match[1]
            condition = {"1": "primary colorectal cancer", "2": "normal colon", "3": "metastasized cancer"}[match[2]]
            if metadata["source_name"] != condition:
                raise ValueError("Title/source tissue disagreement")
        group = groups.setdefault(donor, {})
        if condition in group:
            raise ValueError(f"Duplicate donor-condition: {donor}, {condition}")
        group[condition] = sample
    reference, target = ("baseline", "intravenous LPS") if accession == "GSE87290" else ("normal colon", "primary colorectal cancer")
    available = set(actual)
    if not available.issubset(samples):
        raise ValueError("Count samples absent from GEO metadata")
    pairs, excluded = [], []
    for donor in sorted(groups, key=lambda value: int(value.rsplit("_", 1)[1])):
        group = groups[donor]
        reference_id, target_id = group[reference], group[target]
        if reference_id not in available or target_id not in available:
            excluded.append({"donor_id": donor, "reason": "Incomplete pair in NCBI count matrix", "reference_sample_id": reference_id, "target_sample_id": target_id})
            continue
        if accession == "GSE87290":
            for key in ("sex", "race", "inflammatory response"):
                if samples[reference_id]["characteristics"][key] != samples[target_id]["characteristics"][key]:
                    raise ValueError(f"Within-pair metadata disagreement for {donor}: {key}")
        pairs.append({"donor_id": donor, "reference_sample_id": reference_id, "target_sample_id": target_id,
                      "reference_condition": reference, "target_condition": target})
    expected = 14 if accession == "GSE87290" else 18
    if len(pairs) != expected:
        raise ValueError(f"Pinned {accession} expected {expected} complete pairs, found {len(pairs)}")
    return pairs, excluded


def deterministic_gzip(path: Path, text: str) -> None:
    with path.open("wb") as stream:
        with gzip.GzipFile(filename="", mode="wb", fileobj=stream, mtime=0, compresslevel=9) as compressed:
            compressed.write(text.encode("utf-8"))


def prepare(accession: str, data_dir: Path) -> tuple[np.ndarray, list[str], list[str], dict]:
    raw = data_dir / "raw" / accession
    genes, actual_samples, counts = read_counts(raw / f"{accession}_raw_counts_GRCh38.p13_NCBI.tsv.gz")
    pairs, excluded = pair_samples(accession, soft_samples(raw / "samples_metadata.soft"), actual_samples)
    columns = {sample: column for column, sample in enumerate(actual_samples)}
    selected_samples = [pair[key] for pair in pairs for key in ("reference_sample_id", "target_sample_id")]
    selected = counts[:, [columns[sample] for sample in selected_samples]]
    library_sizes = selected.sum(axis=0, dtype=np.int64)
    if np.any(library_sizes <= 0):
        raise ValueError("Empty sample library")
    # Library denominators use all source genes; filtering never changes them.
    cpm = selected.astype(np.float64) * 1e6 / library_sizes
    minimum_samples = math.ceil(len(selected_samples) / 2)
    eligible = (cpm >= 1.0).sum(axis=1) >= minimum_samples
    feature_ids = [gene for gene, keep in zip(genes, eligible) if keep]
    logged = np.log2(cpm[eligible] + 1.0)
    effects = (logged[:, 1::2] - logged[:, ::2]).T.copy()
    if not np.isfinite(effects).all() or len(feature_ids) < PREDECLARED_MATCHED_GENES:
        raise ValueError("Unexpected prepared feature universe")
    donors = [pair["donor_id"] for pair in pairs]
    destination = data_dir / "prepared" / accession
    destination.mkdir(parents=True, exist_ok=True)
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(["donor_id", *feature_ids])
    for donor, row in zip(donors, effects):
        writer.writerow([donor, *[format(float(value), ".17g") for value in row]])
    effect_path = destination / "effects.csv.gz"
    deterministic_gzip(effect_path, buffer.getvalue())
    gene_path = destination / "gene_ids.txt"
    gene_path.write_text("\n".join(feature_ids) + "\n")
    write_tsv(destination / "donor_pairs.tsv", pairs)
    manifest = json.loads((data_dir / "sources.json").read_text())
    provenance = {
        "accession": accession, "n_source_genes": len(genes), "n_source_count_samples": len(actual_samples),
        "n_selected_samples": len(selected_samples), "n_donors": len(donors), "n_eligible_genes": len(feature_ids),
        "matrix_orientation": "donors_by_genes", "pairs": pairs, "incomplete_pairs_excluded": excluded,
        "source_metadata_samples_absent_from_counts": sorted(set(soft_samples(raw / "samples_metadata.soft")) - set(actual_samples)),
        "other_count_samples_outside_selected_contrast": sorted(set(actual_samples) - set(selected_samples)),
        "library_sizes_all_source_genes": dict(zip(selected_samples, map(int, library_sizes))),
        "filter": {"cpm_minimum": 1.0, "minimum_selected_samples": minimum_samples, "definition": "CPM >= 1 in at least half the selected complete-pair samples"},
        "transform": "log2(CPM + 1), target minus reference within each biological donor",
        "deletion_contract": "Normalized paired effects and feature universe are fixed. Signed paired t recomputes both mean and sample variance after every whole-donor deletion.",
        "eligible_order": "Original NCBI GeneID row order, without outcome-dependent reordering",
        "matched_gene_ids": feature_ids[:PREDECLARED_MATCHED_GENES],
        "effects_sha256": sha256(effect_path.read_bytes()), "gene_ids_sha256": sha256(gene_path.read_bytes()),
        "source_files": [record for record in manifest["files"] if accession in record["relative_path"]],
        "interpretation_limit": "Computational extrema and ranking robustness only; no biological truth, clinical utility, causal effect or differential-expression discovery claim.",
        "normalization_limit": "CPM with fixed pseudocount is a simple reproducible input transform; normalization uncertainty and count-model refitting are outside this benchmark.",
    }
    write_json(destination / "preparation.json", provenance)
    return effects, feature_ids, donors, provenance


def load_prepared(accession: str, data_dir: Path) -> tuple[np.ndarray, list[str], list[str], dict]:
    directory = data_dir / "prepared" / accession
    provenance = json.loads((directory / "preparation.json").read_text())
    payload = (directory / "effects.csv.gz").read_bytes()
    if sha256(payload) != provenance["effects_sha256"]:
        raise ValueError("Prepared effect checksum mismatch")
    gene_payload = (directory / "gene_ids.txt").read_bytes()
    if sha256(gene_payload) != provenance["gene_ids_sha256"]:
        raise ValueError("Prepared gene universe checksum mismatch")
    reader = csv.reader(io.StringIO(gzip.decompress(payload).decode()))
    features = next(reader)[1:]
    donors, rows = [], []
    for row in reader:
        donors.append(row[0])
        rows.append([float(value) for value in row[1:]])
    effects = np.asarray(rows, dtype=np.float64)
    if features != gene_payload.decode().splitlines() or effects.shape != (provenance["n_donors"], provenance["n_eligible_genes"]):
        raise ValueError("Prepared schema mismatch")
    return effects, features, donors, provenance


def compare_t(left, right) -> int:
    """Compare signed paired t across different retained counts, exactly."""
    ls = (left.sum > 0) - (left.sum < 0)
    rs = (right.sum > 0) - (right.sum < 0)
    if ls != rs:
        return (ls > rs) - (ls < rs)
    if not ls:
        return 0
    ld = left.count * left.sum_squares - left.sum ** 2
    rd = right.count * right.sum_squares - right.sum ** 2
    if ld == 0 or rd == 0:
        return ls * ((ld == 0) - (rd == 0))
    a = (left.count - 1) * left.sum ** 2 * rd
    b = (right.count - 1) * right.sum ** 2 * ld
    return ls * ((a > b) - (a < b))


def moment_summary(result) -> tuple[Fraction, Fraction]:
    mean = result.sum / result.count
    variance = (result.sum_squares - result.sum ** 2 / result.count) / (result.count - 1)
    return mean, variance


def benchmark(accession: str, effects: np.ndarray, features: list[str], donors: list[str], output: Path) -> dict:
    from extremarank import compare_scores, extrema

    full_rows, matched_rows, ablation_rows, gene_rows = [], [], [], []
    ndonors, ngenes = effects.shape
    values = [tuple(map(float, effects[:, column])) for column in range(ngenes)]
    full_seconds, full_candidates = 0.0, 0
    full_ablation_rows = []
    full_ablation_counts = {method: {removed: {"minimum_wrong": 0, "maximum_wrong": 0} for removed in range(PREDECLARED_MAX_DELETIONS + 1)} for method in ("contiguous_only", "ends_only")}
    baseline_scores = []
    # Full-spectrum fast algorithm; display summaries are outside the timed call.
    for column, (gene, gene_values) in enumerate(zip(features, values)):
        worst_minimum = worst_maximum = baseline = None
        for removed in range(PREDECLARED_MAX_DELETIONS + 1):
            start = perf_counter()
            result = extrema(gene_values, ndonors - removed, method="two_family")
            elapsed = perf_counter() - start
            full_seconds += elapsed
            full_candidates += result.candidate_count
            if removed == 0:
                baseline = result.minimum
                baseline_scores.append(baseline)
            if worst_minimum is None or compare_t(result.minimum, worst_minimum) < 0:
                worst_minimum = result.minimum
            if worst_maximum is None or compare_t(result.maximum, worst_maximum) > 0:
                worst_maximum = result.maximum
            full_rows.append({"feature_id": gene, "deletions": removed,
                              "retained": ndonors - removed, "minimum_t": result.minimum.t,
                              "maximum_t": result.maximum.t, "candidate_count": result.candidate_count,
                              "seconds": elapsed,
                              "minimum_retained_indices": ",".join(map(str, result.minimum.indices)),
                              "maximum_retained_indices": ",".join(map(str, result.maximum.indices))})
            for method, records in full_ablation_counts.items():
                ablated = extrema(gene_values, ndonors - removed, method=method)
                minimum_wrong = compare_scores(ablated.minimum, result.minimum) != 0
                maximum_wrong = compare_scores(ablated.maximum, result.maximum) != 0
                records[removed]["minimum_wrong"] += minimum_wrong
                records[removed]["maximum_wrong"] += maximum_wrong
                if minimum_wrong or maximum_wrong:
                    full_ablation_rows.append({"feature_id": gene, "deletions": removed, "method": method,
                                               "minimum_wrong": minimum_wrong, "maximum_wrong": maximum_wrong})
        baseline_mean, baseline_variance = moment_summary(baseline)
        worst_mean, worst_variance = moment_summary(worst_minimum)
        removed_donors = [donor for index, donor in enumerate(donors) if index not in worst_minimum.indices]
        gene_rows.append({"feature_id": gene, "baseline_t": baseline.t,
                          "baseline_mean": float(baseline_mean), "baseline_variance": float(baseline_variance),
                          "minimum_t_under_budget": worst_minimum.t, "maximum_t_under_budget": worst_maximum.t,
                          "worst_retained_count": worst_minimum.count, "worst_removed_donors": ",".join(removed_donors),
                          "worst_mean": float(worst_mean), "worst_variance": float(worst_variance),
                          "worst_mean_exact": str(worst_mean), "worst_variance_exact": str(worst_variance),
                          "robust_positive": worst_minimum.sum > 0, "robust_negative": worst_maximum.sum < 0})
        if (column + 1) % 2000 == 0:
            print(f"{accession}: full spectrum {column + 1}/{ngenes} genes", flush=True)
    write_tsv(output / f"{accession}_full_extrema.tsv", full_rows)
    write_tsv(output / f"{accession}_gene_audit.tsv", gene_rows)
    write_tsv(output / f"{accession}_full_ablation_failures.tsv", full_ablation_rows)
    def baseline_order(left: int, right: int) -> int:
        comparison = compare_scores(baseline_scores[left], baseline_scores[right])
        return -comparison if comparison else ((features[left] > features[right]) - (features[left] < features[right]))
    top20_indices = sorted(range(ngenes), key=cmp_to_key(baseline_order))[:20]
    write_tsv(output / f"{accession}_baseline_top20.tsv", [{"baseline_rank": rank, **gene_rows[index]} for rank, index in enumerate(top20_indices, 1)])
    method_seconds = {method: 0.0 for method in ("two_family", "exhaustive")}
    method_candidates = {method: 0 for method in method_seconds}
    ablation_errors = {method: {"minimum_wrong": 0, "maximum_wrong": 0, "both_equal": 0} for method in ("contiguous_only", "ends_only")}
    mismatches = 0
    # Pair by gene and budget. Alternate execution order deterministically to
    # balance first-call/cache effects; every timed call includes input sorting.
    for column, (gene, gene_values) in enumerate(zip(features[:PREDECLARED_MATCHED_GENES], values[:PREDECLARED_MATCHED_GENES])):
        for removed in range(PREDECLARED_MAX_DELETIONS + 1):
            results, timings = {}, {}
            order = ("two_family", "exhaustive") if (column + removed) % 2 == 0 else ("exhaustive", "two_family")
            for method in order:
                start = perf_counter()
                results[method] = extrema(gene_values, ndonors - removed, method=method)
                timings[method] = perf_counter() - start
                method_seconds[method] += timings[method]
                method_candidates[method] += results[method].candidate_count
            fast, oracle = results["two_family"], results["exhaustive"]
            minimum_equal = compare_scores(fast.minimum, oracle.minimum) == 0
            maximum_equal = compare_scores(fast.maximum, oracle.maximum) == 0
            mismatches += not (minimum_equal and maximum_equal)
            matched_rows.append({"feature_id": gene, "deletions": removed, "retained": ndonors - removed,
                                 "minimum_equal": minimum_equal, "maximum_equal": maximum_equal,
                                 "two_family_seconds": timings["two_family"], "exhaustive_seconds": timings["exhaustive"],
                                 "two_family_candidates": fast.candidate_count, "exhaustive_candidates": oracle.candidate_count,
                                 "minimum_witness_indices_equal": fast.minimum.indices == oracle.minimum.indices,
                                 "maximum_witness_indices_equal": fast.maximum.indices == oracle.maximum.indices})
            for method in ablation_errors:
                result = extrema(gene_values, ndonors - removed, method=method)
                min_wrong = compare_scores(result.minimum, oracle.minimum) != 0
                max_wrong = compare_scores(result.maximum, oracle.maximum) != 0
                ablation_errors[method]["minimum_wrong"] += min_wrong
                ablation_errors[method]["maximum_wrong"] += max_wrong
                ablation_errors[method]["both_equal"] += not (min_wrong or max_wrong)
                ablation_rows.append({"feature_id": gene, "deletions": removed, "method": method,
                                      "minimum_wrong": min_wrong, "maximum_wrong": max_wrong,
                                      "candidate_count": result.candidate_count})
    write_tsv(output / f"{accession}_matched_exact.tsv", matched_rows)
    write_tsv(output / f"{accession}_ablations.tsv", ablation_rows)
    summary = {
        "accession": accession, "n_donors": ndonors, "n_eligible_genes": ngenes,
        "deletions_tested": list(range(PREDECLARED_MAX_DELETIONS + 1)),
        "matched_gene_ids": features[:PREDECLARED_MATCHED_GENES], "matched_gene_order": "First100 eligible source rows",
        "full_two_family": {"gene_budget_calls": len(full_rows), "seconds": full_seconds, "candidate_count": full_candidates},
        "matched_exact": {"gene_budget_calls": len(matched_rows), "mismatches": int(mismatches),
                          "seconds": method_seconds, "candidate_count": method_candidates,
                          "exhaustive_over_two_family_runtime_ratio": method_seconds["exhaustive"] / method_seconds["two_family"]},
        "ablations": {method: {**record, "gene_budget_calls": len(matched_rows),
                                "minimum_wrong_rate": record["minimum_wrong"] / len(matched_rows),
                                "maximum_wrong_rate": record["maximum_wrong"] / len(matched_rows)} for method, record in ablation_errors.items()},
        "full_spectrum_ablations": {method: {str(removed): {**record, "n_genes": ngenes,
                                                            "minimum_wrong_rate": record["minimum_wrong"] / ngenes,
                                                            "maximum_wrong_rate": record["maximum_wrong"] / ngenes} for removed, record in records.items()} for method, records in full_ablation_counts.items()},
        "practical_gene_audit": {"baseline_top20_gene_ids": [features[index] for index in top20_indices],
                                 "robust_positive_genes": sum(row["robust_positive"] for row in gene_rows),
                                 "robust_negative_genes": sum(row["robust_negative"] for row in gene_rows),
                                 "baseline_top20_robust_positive": sum(gene_rows[index]["robust_positive"] for index in top20_indices),
                                 "output": f"{accession}_gene_audit.tsv",
                                 "definition": "A robust-positive gene has positive paired mean after every subset of at most4 donor deletions. This is a descriptive stability screen, not significance or biological truth.",
                                 "worst_case": "The smallest signed paired t across every retained count n,n-1,...,n-4; exact signed-rational comparison across retained counts. Mean, sample variance and actual removed-donor witness are supplied."},
        "timing_scope": "Entire extrema call including sorting and exact arithmetic; serialization, preprocessing and ablation calls excluded from matched runtime. One alternating-order matched run; ratios depend on hardware and Python.",
        "exactness_scope": "Score comparison on the exact binary64 paired-effect input; no floating tolerance used for equality. Equal extrema values may have different valid retained-index witnesses.",
        "validation_status": "passed" if mismatches == 0 else "FAILED",
    }
    write_json(output / f"{accession}_summary.json", summary)
    if mismatches:
        raise AssertionError(f"{accession}: {mismatches} exact extrema mismatches")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=REPO / "data")
    parser.add_argument("--output-dir", type=Path, default=REPO / "results" / "real")
    parser.add_argument("--dataset", choices=(*DATASETS, "all"), default="all")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--prepared-only", action="store_true", help="Read checksummed paired effects, without original source files")
    parser.add_argument("--download-missing", action="store_true", help="Restore missing pinned raw sources; default is fully offline")
    args = parser.parse_args()
    if args.prepare_only and args.prepared_only:
        parser.error("--prepare-only and --prepared-only cannot be combined")
    data_dir = args.data_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source_check = [] if args.prepared_only else source_verifier()(data_dir, offline=not args.download_missing)
    chosen = DATASETS if args.dataset == "all" else (args.dataset,)
    results, inputs = [], []
    for accession in chosen:
        effects, features, donors, provenance = load_prepared(accession, data_dir) if args.prepared_only else prepare(accession, data_dir)
        inputs.append({"accession": accession, "effects_sha256": provenance["effects_sha256"],
                       "gene_ids_sha256": provenance["gene_ids_sha256"], "shape": list(effects.shape), "donor_ids": donors})
        print(f"{accession}: prepared {len(donors)} real pairs, {len(features)} eligible genes", flush=True)
        if not args.prepare_only:
            results.append(benchmark(accession, effects, features, donors, output))
    environment = {"python": sys.version, "python_executable": sys.executable,
                   "numpy": np.__version__, "platform": platform.platform(), "machine": platform.machine(),
                   "processor": platform.processor(), "retrieved_sources_utc": "2026-10-03",
                   "run_at_utc": datetime.now(timezone.utc).isoformat()}
    write_json(output / "real_benchmark.json", {"environment": environment, "source_checks": source_check,
                                              "inputs": inputs, "summaries": results,
                                              "predeclared_contract": {"max_deletions": PREDECLARED_MAX_DELETIONS,
                                                                       "matched_first_eligible_genes": PREDECLARED_MATCHED_GENES,
                                                                       "cpm_threshold": 1.0, "pseudocount": 1.0},
                                              "interpretation_limit": "Exact computational input robustness; no biological truth or clinical claims."})
    print(json.dumps({"datasets": inputs, "matched_validation": [record["validation_status"] for record in results]}, indent=2), flush=True)


if __name__ == "__main__":
    main()
