"""CSV/TSV interface for exact paired-effect deletion audits (stdlib only)."""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from fractions import Fraction
from functools import cmp_to_key
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import tempfile
import time
from typing import Sequence
import zlib

from . import __version__
from .audit import audit_topk
from .core import PreparedValues, Score, compare_abs_scores, compare_t


@dataclass(frozen=True)
class InputTable:
    donor_ids: tuple[str, ...]
    gene_ids: tuple[str, ...]
    effects: tuple[tuple[float, ...], ...]
    sha256: str
    delimiter: str
    compression: str


def read_effects(path: Path) -> InputTable:
    """Read finite donor x gene effects, preserving the input donor order."""
    raw = path.read_bytes()
    compressed = path.suffix.lower() == ".gz"
    suffix = path.with_suffix("").suffix.lower() if compressed else path.suffix.lower()
    if compressed:
        try:
            decoded = gzip.decompress(raw)
        except (OSError, EOFError, zlib.error) as exc:
            raise ValueError("invalid or truncated gzip input") from exc
    else:
        decoded = raw
    try:
        source = decoded.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ValueError("input must be UTF-8 CSV or TSV") from exc
    if suffix == ".csv":
        delimiter = ","
    elif suffix == ".tsv":
        delimiter = "\t"
    else:
        try:
            delimiter = csv.Sniffer().sniff(source[:65536], delimiters=",\t").delimiter
        except csv.Error as exc:
            raise ValueError("cannot identify CSV/TSV delimiter; use a .csv or .tsv filename") from exc
    reader = csv.reader(io.StringIO(source, newline=""), delimiter=delimiter, strict=True)
    try:
        header = next(reader, None)
        if header is None or len(header) < 2 or header[0].strip() != "donor_id":
            raise ValueError("first header column must be donor_id, followed by at least one gene ID")
        gene_ids = tuple(x.strip() for x in header[1:])
        if any(not x for x in gene_ids) or len(set(gene_ids)) != len(gene_ids):
            raise ValueError("gene IDs must be nonempty and unique")
        donor_ids: list[str] = []
        effects: list[tuple[float, ...]] = []
        seen_donors: set[str] = set()
        for row in reader:
            if not row or all(not cell.strip() for cell in row):
                continue
            if len(row) != len(header):
                raise ValueError(f"line {reader.line_num}: expected {len(header)} columns, got {len(row)}")
            donor_id = row[0].strip()
            if not donor_id or donor_id in seen_donors:
                raise ValueError(f"line {reader.line_num}: donor IDs must be nonempty and unique")
            values: list[float] = []
            for gene_id, cell in zip(gene_ids, row[1:]):
                try:
                    value = float(cell)
                except ValueError as exc:
                    raise ValueError(f"line {reader.line_num}, gene {gene_id!r}: effect must be numeric") from exc
                if not math.isfinite(value):
                    raise ValueError(f"line {reader.line_num}, gene {gene_id!r}: effect must be finite")
                values.append(value)
            seen_donors.add(donor_id)
            donor_ids.append(donor_id)
            effects.append(tuple(values))
    except csv.Error as exc:
        raise ValueError(f"invalid CSV/TSV near line {reader.line_num}: {exc}") from exc
    if len(donor_ids) < 2:
        raise ValueError("input must contain at least two donors")
    return InputTable(tuple(donor_ids), gene_ids, tuple(effects),
                      hashlib.sha256(raw).hexdigest(), delimiter, "gzip" if compressed else "none")


def _display(value: float | Fraction) -> str:
    """Round-trip binary64 display; exact fractions are rounded for display only."""
    try:
        number = float(value)
    except OverflowError:
        return "inf" if value > 0 else "-inf"
    return format(number, ".17g") if math.isfinite(number) else str(number)


def _json_safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def _deleted(score: Score, donor_ids: tuple[str, ...]) -> list[str]:
    retained = set(score.indices)
    return [donor_id for i, donor_id in enumerate(donor_ids) if i not in retained]


def _ids_cell(ids: list[str]) -> str:
    return json.dumps(ids, ensure_ascii=False, separators=(",", ":"))


def _write_tsv(path: Path, columns: Sequence[str], rows: list[dict]) -> None:
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                     dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        try:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _write_json(path: Path, data: dict) -> None:
    text = json.dumps(_json_safe(data), indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="\n",
                                     dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        try:
            handle.write(text)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


SUMMARY_COLUMNS = (
    "gene_id", "baseline_rank", "baseline_t", "worst_t", "best_t",
    "worst_deleted_donors", "best_deleted_donors", "worst_retained_count",
    "best_retained_count", "worst_mean", "worst_sample_variance",
    "robust_positive", "robust_negative", "baseline_topk",
)
ENVELOPE_COLUMNS = (
    "gene_id", "deleted_count", "retained_count", "minimum_t", "maximum_t",
    "minimum_score", "maximum_score", "minimum_deleted_donors", "maximum_deleted_donors",
    "minimum_sum", "minimum_sum_squares", "maximum_sum", "maximum_sum_squares",
    "candidate_count", "total_subsets",
)


def _summarize(table: InputTable, budget: int, k: int, direction: str) -> tuple[list[dict], list[dict], list[Score], list[int]]:
    d, g = len(table.donor_ids), len(table.gene_ids)
    summaries, envelopes, baseline_scores = [], [], []
    for j, gene_id in enumerate(table.gene_ids):
        prepared = PreparedValues(row[j] for row in table.effects)
        worst = best = None
        robust_positive = True
        robust_negative = True
        for removed_count in range(budget + 1):
            envelope = prepared.extrema(d - removed_count)
            lo, hi = envelope.minimum, envelope.maximum
            if removed_count == 0:
                baseline_scores.append(lo)
            if worst is None or compare_t(lo, worst) < 0:
                worst = lo
            if best is None or compare_t(hi, best) > 0:
                best = hi
            robust_positive = robust_positive and lo.sum > 0
            robust_negative = robust_negative and hi.sum < 0
            envelopes.append({
                "gene_id": gene_id, "deleted_count": removed_count, "retained_count": d - removed_count,
                "minimum_t": _display(lo.t), "maximum_t": _display(hi.t),
                "minimum_score": _display(lo.value), "maximum_score": _display(hi.value),
                "minimum_deleted_donors": _ids_cell(_deleted(lo, table.donor_ids)),
                "maximum_deleted_donors": _ids_cell(_deleted(hi, table.donor_ids)),
                "minimum_sum": str(lo.sum), "minimum_sum_squares": str(lo.sum_squares),
                "maximum_sum": str(hi.sum), "maximum_sum_squares": str(hi.sum_squares),
                "candidate_count": envelope.candidate_count, "total_subsets": envelope.total_subsets,
            })
        assert worst is not None and best is not None
        mean = worst.sum / worst.count
        variance = (worst.sum_squares - worst.sum ** 2 / worst.count) / (worst.count - 1)
        summaries.append({
            "gene_id": gene_id, "baseline_rank": None, "baseline_t": _display(baseline_scores[-1].t),
            "worst_t": _display(worst.t), "best_t": _display(best.t),
            "worst_deleted_donors": _ids_cell(_deleted(worst, table.donor_ids)),
            "best_deleted_donors": _ids_cell(_deleted(best, table.donor_ids)),
            "worst_retained_count": worst.count, "best_retained_count": best.count,
            "worst_mean": _display(mean), "worst_sample_variance": _display(variance),
            "robust_positive": robust_positive, "robust_negative": robust_negative, "baseline_topk": False,
        })

    def compare_gene(a: int, b: int) -> int:
        if direction == "absolute":
            # Baseline retained counts agree, so |R| and |t| induce the same rank.
            comparison = compare_abs_scores(baseline_scores[a], baseline_scores[b])
        else:
            comparison = compare_t(baseline_scores[a], baseline_scores[b])
            if direction == "down":
                comparison = -comparison
        if comparison:
            return -comparison
        return (table.gene_ids[a] > table.gene_ids[b]) - (table.gene_ids[a] < table.gene_ids[b])

    order = sorted(range(g), key=cmp_to_key(compare_gene))
    for rank, j in enumerate(order, 1):
        summaries[j]["baseline_rank"] = rank
        summaries[j]["baseline_topk"] = rank <= k
    return [summaries[j] for j in order], envelopes, baseline_scores, order


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="extremarank", description="Exact deletion envelopes and shared-donor top-K audit for paired effects.")
    parser.add_argument("input", type=Path, help="UTF-8 donor x gene CSV/TSV, optionally gzip; first header must be donor_id")
    parser.add_argument("--budget", type=int, default=2, help="maximum deleted donors (default: 2); must leave >=2 donors")
    parser.add_argument("--top-k", type=int, default=None,
                        help="top-K set size (default: 20, reduced to gene count only when omitted)")
    parser.add_argument("--direction", choices=("up", "down", "absolute"), default="up",
                        help="rank by largest signed t, smallest signed t, or largest |t| (default: up)")
    parser.add_argument("--output", type=Path, required=True, help="output directory (existing result filenames are replaced)")
    parser.add_argument("--max-nodes", type=int, default=10000, help="shared top-K search node limit (default: 10000)")
    parser.add_argument("--skip-topk", action="store_true", help="compute every gene's exact envelopes; do not run shared top-K search")
    parser.add_argument("--version", action="version", version=f"ExtremaRank {__version__}")
    args = parser.parse_args(argv)
    start = time.perf_counter()
    try:
        table = read_effects(args.input)
        d, g = len(table.donor_ids), len(table.gene_ids)
        if not 0 <= args.budget <= d - 2:
            raise ValueError(f"--budget must be between 0 and {d - 2} for {d} donors; at least two must remain")
        k = min(20, g) if args.top_k is None else args.top_k
        if not 1 <= k <= g:
            raise ValueError(f"--top-k must be between 1 and {g}; only the omitted default is reduced")
        if args.max_nodes < 1:
            raise ValueError("--max-nodes must be positive")
        summaries, envelopes, baseline_scores, order = _summarize(table, args.budget, k, args.direction)
        envelopes_seconds = time.perf_counter() - start
        baseline_topk = [table.gene_ids[j] for j in order[:k]]
        if args.skip_topk:
            audit = {"status": "NOT_RUN", "baseline_topk": baseline_topk,
                     "reason": "--skip-topk was selected", "nodes": 0,
                     "deleted_indices": None, "deleted_donors": None, "witness_topk": None}
        else:
            result = audit_topk(table.effects, table.gene_ids, k, args.budget,
                               max_nodes=args.max_nodes, direction=args.direction)
            audit = result.as_dict()
            audit["deleted_donors"] = ([table.donor_ids[i] for i in result.deleted_indices]
                                       if result.deleted_indices is not None else None)
            if list(result.baseline_topk) != baseline_topk:
                raise RuntimeError("baseline ranking differs between envelope and shared top-K audit")
        completed = time.perf_counter() - start
        report = {
            "schema_version": 1, "software": {"name": "ExtremaRank", "version": __version__},
            "status": audit["status"],
            "input": {"path": str(args.input.resolve()), "sha256": table.sha256,
                      "sha256_scope": "original input file bytes, before gzip decompression",
                      "delimiter": "tab" if table.delimiter == "\t" else "comma",
                      "compression": table.compression,
                      "n_donors": d, "n_genes": g, "donor_ids": list(table.donor_ids)},
            "parameters": {"budget": args.budget, "top_k_requested": args.top_k,
                           "top_k": k, "top_k_default_reduced": args.top_k is None and g < 20,
                           "direction": args.direction, "max_nodes": args.max_nodes, "skip_topk": args.skip_topk},
            "tie_policy": {
                "gene_ranks": "exact ranking-score ties use ascending lexical gene ID",
                "donor_value_ties": "equal effects are sorted by original input row index",
                "within_cardinality_extrema": "first attaining two-family candidate in deterministic traversal order",
                "across_cardinality_extrema": "equal extended paired-t scores keep the fewest deletions",
            },
            "target": f"unordered baseline {args.direction} paired-t top-K gene set under any shared deletion of <=budget donors",
            "arithmetic": "exact moments and comparisons of parsed finite binary64 effects; decimal scores are display only",
            "score_conventions": {"all_zero": "descriptive t=0", "nonzero_constant": "signed infinity",
                                  "variance": "sample variance, denominator retained_count-1"},
            "baseline_topk_scores": [{"gene_id": table.gene_ids[j], "t": baseline_scores[j].t}
                                     for j in order[:k]],
            "topk_audit": audit,
            "outputs": {"summary": "summary.tsv", "exact_cardinality_envelopes": "envelopes.tsv"},
            "timing_seconds": {"input_and_envelopes": envelopes_seconds,
                               "shared_topk": completed - envelopes_seconds, "total_before_writing": completed},
        }
        args.output.mkdir(parents=True, exist_ok=True)
        _write_tsv(args.output / "summary.tsv", SUMMARY_COLUMNS, summaries)
        _write_tsv(args.output / "envelopes.tsv", ENVELOPE_COLUMNS, envelopes)
        _write_json(args.output / "audit.json", report)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(f"{audit['status']}: {g} genes, {d} donors, deletion budget {args.budget}, top-K {k}")
    print(f"Results: {args.output.resolve()}")
    return 0
