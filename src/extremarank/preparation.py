"""Strict matrix/metadata joins and explicit paired or independent contrasts."""
from __future__ import annotations

import csv
from dataclasses import dataclass
import gzip
import hashlib
import io
import math
from pathlib import Path
import zlib


def read_table(path: Path):
    raw = path.read_bytes()
    try:
        decoded = gzip.decompress(raw) if path.suffix.lower() == ".gz" else raw
        text = decoded.decode("utf-8-sig")
    except (OSError, EOFError, UnicodeDecodeError, zlib.error) as exc:
        raise ValueError("invalid UTF-8 or gzip table") from exc
    suffix = path.with_suffix("").suffix.lower() if path.suffix.lower() == ".gz" else path.suffix.lower()
    delimiter = "," if suffix == ".csv" else "\t" if suffix == ".tsv" else None
    if delimiter is None:
        try:
            delimiter = csv.Sniffer().sniff(text[:65536], delimiters=",\t").delimiter
        except csv.Error as exc:
            raise ValueError("use a CSV or TSV filename") from exc
    try:
        rows = list(csv.reader(io.StringIO(text, newline=""), delimiter=delimiter, strict=True))
    except csv.Error as exc:
        raise ValueError("invalid CSV/TSV") from exc
    rows = [row for row in rows if row and any(x.strip() for x in row)]
    if not rows:
        raise ValueError("empty table")
    header = tuple(x.strip() for x in rows[0])
    if any(not x for x in header) or len(set(header)) != len(header):
        raise ValueError("header IDs must be nonempty and unique")
    if any(len(row) != len(header) for row in rows[1:]):
        raise ValueError("table contains ragged rows")
    return header, rows[1:], hashlib.sha256(raw).hexdigest()


@dataclass(frozen=True)
class StudyInput:
    unit_ids: tuple[str, ...]
    feature_ids: tuple[str, ...]
    values: tuple[tuple[float, ...], ...]
    groups: tuple[str, ...]
    unit_metadata: tuple[dict, ...]
    manifest: dict


def prepare_study(matrix: Path, metadata: Path, target: str, reference: str,
                  design: str = "paired", orientation: str = "features",
                  transform: str = "none", missing: str = "error",
                  incomplete_pairs: str = "error", min_cpm: float = 0.0,
                  pseudocount: float = 1.0) -> StudyInput:
    if design not in {"paired", "welch"} or orientation not in {"features", "samples"}:
        raise ValueError("design must be paired/welch and orientation features/samples")
    if transform not in {"none", "log2", "cpm-log2"} or missing not in {"error", "drop-features"}:
        raise ValueError("invalid transformation or missing-value policy")
    if incomplete_pairs not in {"error", "drop"}:
        raise ValueError("incomplete-pairs must be error or drop")
    if not target or not reference or target == reference:
        raise ValueError("declare distinct target and reference groups")
    if not math.isfinite(pseudocount) or pseudocount <= 0 or not math.isfinite(min_cpm) or min_cpm < 0:
        raise ValueError("pseudocount must be positive finite; min-cpm nonnegative finite")
    if min_cpm and transform != "cpm-log2":
        raise ValueError("min-cpm requires the cpm-log2 transformation")
    h, rows, matrix_hash = read_table(matrix)
    expected = "feature_id" if orientation == "features" else "sample_id"
    if len(h) < 2 or h[0] != expected:
        raise ValueError(f"first matrix header must be {expected} for this orientation")
    row_ids = tuple(row[0].strip() for row in rows)
    if not row_ids or any(not x for x in row_ids) or len(set(row_ids)) != len(row_ids):
        raise ValueError("matrix row IDs must be nonempty and unique")
    feature_ids = row_ids if orientation == "features" else h[1:]
    sample_ids = h[1:] if orientation == "features" else row_ids
    parsed = []
    for row in rows:
        values = []
        for cell in row[1:]:
            cell = cell.strip()
            if cell.lower() in {"", "na", "nan", "null", "none"}:
                values.append(None)
                continue
            try:
                value = float(cell)
            except ValueError as exc:
                raise ValueError("matrix values must be numeric or explicit missing values") from exc
            if not math.isfinite(value):
                raise ValueError("infinite matrix values are not supported")
            values.append(value)
        parsed.append(values)
    values_by_sample = [list(col) for col in zip(*parsed)] if orientation == "features" else parsed
    mh, mrows, metadata_hash = read_table(metadata)
    if not {"sample_id", "group"} <= set(mh) or (design == "paired" and "donor_id" not in mh):
        raise ValueError("metadata requires sample_id and group; paired design also requires donor_id")
    meta = {}
    for row in mrows:
        record = dict(zip(mh, (x.strip() for x in row)))
        sid = record["sample_id"]
        if not sid or sid in meta or not record["group"]:
            raise ValueError("metadata sample IDs must be unique/nonempty and groups nonempty")
        record.setdefault("donor_id", sid)
        if not record["donor_id"]:
            raise ValueError("metadata donor IDs must be nonempty")
        meta[sid] = record
    if set(meta) != set(sample_ids):
        raise ValueError("matrix and metadata sample ID sets must match exactly; no silent intersection")
    chosen = [i for i, sid in enumerate(sample_ids) if meta[sid]["group"] in {target, reference}]
    if set(meta[sample_ids[i]]["group"] for i in chosen) != {target, reference}:
        raise ValueError("both declared groups must be present")
    chosen_set = set(chosen)
    excluded_groups = [sample_ids[i] for i in range(len(sample_ids)) if i not in chosen_set]
    pairs = {}
    dropped_pairs = []
    if design == "paired":
        for i in chosen:
            record = meta[sample_ids[i]]
            pair = pairs.setdefault(record["donor_id"], {})
            if record["group"] in pair:
                raise ValueError("each donor must have exactly one sample per contrast group; aggregate technical replicates first")
            pair[record["group"]] = i
        dropped_pairs = [donor for donor, pair in pairs.items() if len(pair) != 2]
        if dropped_pairs and incomplete_pairs == "error":
            raise ValueError("incomplete donor pairs; use --incomplete-pairs drop to explicitly exclude and record them")
        pairs = {donor: pair for donor, pair in pairs.items() if len(pair) == 2}
        keep_samples = {i for pair in pairs.values() for i in pair.values()}
        chosen = [i for i in chosen if i in keep_samples]
        if len(pairs) < 2:
            raise ValueError("at least two complete donor pairs are required")
    else:
        donors = [meta[sample_ids[i]]["donor_id"] for i in chosen]
        if len(set(donors)) != len(donors):
            raise ValueError("Welch design requires independent donor units; repeated donor IDs require a paired or other model")
        if any(sum(meta[sample_ids[i]]["group"] == group for i in chosen) < 2 for group in (target, reference)):
            raise ValueError("Welch design requires >=2 independent samples in each group")
    chosen_values = [values_by_sample[i][:] for i in chosen]
    missing_features = [j for j in range(len(feature_ids)) if any(row[j] is None for row in chosen_values)]
    if missing_features and (missing == "error" or transform == "cpm-log2"):
        raise ValueError("missing values present; drop-features is explicit and unavailable for CPM library totals")
    missing_set = set(missing_features)
    retained = [j for j in range(len(feature_ids)) if j not in missing_set]
    library_sizes = None
    if transform != "none":
        if any(row[j] < 0 for row in chosen_values for j in retained):
            raise ValueError("log2/CPM transformations require nonnegative input values")
        if transform == "cpm-log2":
            try:
                library_sizes = [math.fsum(row) for row in chosen_values]
            except OverflowError as exc:
                raise ValueError("library totals exceed finite binary64 range") from exc
            if any(not math.isfinite(total) or total <= 0 for total in library_sizes):
                raise ValueError("each sample must have a positive finite library total")
            cpm = [[x / total * 1e6 for x in row] for row, total in zip(chosen_values, library_sizes)]
            if min_cpm:
                retained = [j for j in retained if sum(row[j] >= min_cpm for row in cpm) >= math.ceil(len(chosen)/2)]
            chosen_values = cpm
        for row in chosen_values:
            for j in retained:
                row[j] = math.log2(row[j] + pseudocount)
    if not retained:
        raise ValueError("no features remain after declared exclusions")
    frozen = [[row[j] for j in retained] for row in chosen_values]
    if any(not math.isfinite(x) for row in frozen for x in row):
        raise ValueError("transformation produced nonfinite values")
    chosen_map = {source_i: k for k, source_i in enumerate(chosen)}
    if design == "paired":
        units = tuple(pairs)
        values = tuple(tuple(frozen[chosen_map[pair[target]]][j] - frozen[chosen_map[pair[reference]]][j]
                             for j in range(len(retained))) for pair in pairs.values())
        groups = ()
        unit_metadata = tuple({"donor_id": donor, "target_sample": sample_ids[pair[target]],
                               "reference_sample": sample_ids[pair[reference]]} for donor, pair in pairs.items())
    else:
        units = tuple(meta[sample_ids[i]]["donor_id"] for i in chosen)
        values = tuple(tuple(row) for row in frozen)
        groups = tuple(meta[sample_ids[i]]["group"] for i in chosen)
        unit_metadata = tuple(meta[sample_ids[i]] for i in chosen)
    if any(not math.isfinite(x) for row in values for x in row):
        raise ValueError("contrast subtraction produced nonfinite values")
    retained_set = set(retained)
    manifest = {"schema_version": 1, "design": design, "contrast": f"{target} - {reference}",
        "target": target, "reference": reference, "orientation": orientation,
        "transform": transform, "pseudocount": pseudocount if transform != "none" else None,
        "min_cpm": min_cpm, "cpm_filter": ">=min_cpm in >=half included samples; all source features define library totals",
        "matrix_sha256": matrix_hash, "metadata_sha256": metadata_hash,
        "matrix_source_features": len(feature_ids), "retained_features": len(retained),
        "included_samples": [sample_ids[i] for i in chosen], "independent_units": len(units),
        "excluded_other_group_samples": excluded_groups, "excluded_incomplete_donors": dropped_pairs,
        "missing_policy": missing, "excluded_missing_features": [feature_ids[j] for j in missing_features],
        "excluded_features": [x for j, x in enumerate(feature_ids) if j not in retained_set],
        "library_totals": dict(zip((sample_ids[i] for i in chosen), library_sizes)) if library_sizes else None,
        "unit_metadata": unit_metadata,
        "identity_check": "strict source IDs only; raw genotype/sample-swap verification NOT_EVALUATED",
        "guarantee_input": "finite binary64 prepared values; preprocessing frozen before deletion",
        "inference": "descriptive ranking sensitivity; no confounder adjustment or FDR inference"}
    return StudyInput(units, tuple(feature_ids[j] for j in retained), values, groups, unit_metadata, manifest)
