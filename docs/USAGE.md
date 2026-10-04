# Use ExtremaRank on paired omics effects

For v0.2 matrix/metadata preparation, independent groups, external models and
offline sample-influence reports, use the [study workflows](STUDY.md).
This page retains the original paired-effects interface.

ExtremaRank asks how donor deletion can change a gene's signed paired t score or a selected top-K gene set. Every retained subset recomputes both the mean and the sample variance. The scalar envelopes are exact over all allowed subsets of the parsed binary64 effects. The shared top-K search returns a certificate, a concrete common deletion witness, or an unresolved result at its node limit.

This is a descriptive robustness audit. It does not calculate differential-expression p values or control false discovery rates. Effects and preprocessing remain frozen during the audit: deletion does not rerun normalization, covariate adjustment, or expression filtering. Start with biologically meaningful paired donor effects, such as a normalized within-donor treatment-minus-control contrast; different preprocessing choices define different inputs.

## Install and run the included example

From the package directory, install the dependency-free Python package:

```sh
python -m pip install -e .
python -m extremarank examples/paired_effects.tsv --budget 2 --top-k 3 --output results/example_cli
```

The installed `extremarank` command accepts the same arguments:

```sh
extremarank examples/paired_effects.tsv --budget 2 --top-k 3 --output results/example_cli
```

An uninstalled source checkout also works:

```sh
PYTHONPATH=src python -m extremarank examples/paired_effects.tsv --budget 2 --top-k 3 --output results/example_cli
```

The example has six donors, stable positive and negative effects, one-donor and two-donor driven effects, a positive constant, and an all-zero feature. It is synthetic and illustrates the file interface; it does not establish biological performance.

## Input format

Use UTF-8 CSV or TSV, optionally gzip-compressed, with a header and one row per donor. The first column must be named `donor_id`. Every later column names a gene or other feature. IDs must be nonempty and unique; leading and trailing whitespace around IDs is removed. Preserve stable feature IDs in the header because exact ranking ties use their lexical order.

```text
donor_id	GENE_A	GENE_B
D01	0.75	-0.25
D02	0.5	-0.125
D03	1	0
D04	0.625	-0.5
```

In this block, the visible gaps represent tabs. Use comma separation for a `.csv` or `.csv.gz` file and tab separation for `.tsv` or `.tsv.gz`. Quoted fields, UTF-8 IDs, and a UTF-8 byte-order mark are accepted. For other filename extensions the delimiter is inferred from commas and tabs after any gzip decompression. Empty rows are ignored. Ragged rows, missing effects, duplicate IDs, NaN, infinity, and corrupt or truncated gzip streams are rejected. Decimal effects are converted once to finite binary64 numbers, and decisions use exact integer moments of those numbers. When generating an input file from Python floats, `format(value, ".17g")` preserves their binary64 value on re-reading.

At least two donors must remain after every deletion. Consequently, `--budget` must be between zero and the donor count minus two. Its default is two; with fewer than four donors, specify a smaller budget explicitly.

## Arguments and ranking direction

| Argument | Meaning |
|---|---|
| `input` | CSV/TSV file containing frozen paired effects, optionally `.gz` compressed. |
| `--output DIRECTORY` | Required output directory. Existing `summary.tsv`, `envelopes.tsv`, and `audit.json` are replaced. |
| `--budget B` | Audit every deletion count from zero through B; default 2. |
| `--top-k K` | Baseline top-K set size; default 20. When omitted, the default is reduced to the feature count for a smaller input. An explicit oversized K is rejected. |
| `--direction up` | Rank by largest signed paired t; default. |
| `--direction down` | Rank by smallest signed paired t. |
| `--direction absolute` | Rank by largest absolute paired t; opposite signs with equal magnitudes are tied. |
| `--max-nodes N` | Maximum shared-search nodes across deletion counts; default 10000. |
| `--skip-topk` | Compute all exact scalar envelopes and baseline ranks, but omit the shared top-K search. Its status is `NOT_RUN`. |

Examples for negative effects and for large inputs:

```sh
extremarank effects.csv --direction down --budget 3 --top-k 25 --output results/down
extremarank effects.tsv --direction absolute --budget 4 --top-k 25 --max-nodes 50000 --output results/absolute
extremarank effects.tsv --budget 4 --skip-topk --output results/envelopes_only
```

The bundled prepared matrices can be used directly, without manually decompressing them:

```sh
extremarank data/prepared/GSE87290/effects.csv.gz --budget 4 --top-k 20 --output results/GSE87290_cli
extremarank data/prepared/GSE50760/effects.csv.gz --budget 4 --top-k 20 --output results/GSE50760_cli
```

These commands audit the full eligible feature universe of each prepared matrix. Dataset pairing, exclusion rules, filtering, normalization, and source provenance are recorded in the corresponding `preparation.json` files and [data attribution](../data/DATA_LICENSE.md). Interpret these as reproducible computational examples under that fixed preprocessing contract.

One local sequential run, including gzip input and all output writes, took 2.85 seconds for GSE87290 (14 donors, 13,328 genes) and 3.94 seconds for GSE50760 (18 donors, 16,597 genes) at budget 4 and top-K 20. Both returned `REFUTED` with independently checked shared deletion witnesses. These are single-run measurements on the recorded machine; [the end-to-end record](../results/cli_end_to_end.json) includes provenance, output dimensions, and timing details.

All scalar minimum and maximum scores remain signed for every direction. In absolute mode they are signed interval endpoints, not exact minimum and maximum absolute scores. The shared absolute-rank search uses conservative bounds when a signed interval crosses zero; an overlap can require further search. The minimum absolute score is not claimed to be attained at a signed endpoint.

## Read the outputs

`summary.tsv` has one row per gene, sorted by baseline rank. Its columns are:

| Columns | Interpretation |
|---|---|
| `gene_id`, `baseline_rank`, `baseline_t` | Gene ID, one-based rank in the chosen direction, and signed t using all donors. |
| `worst_t`, `best_t` | Smallest and largest **signed** t over every allowed deletion count from zero through the budget. Exact comparisons account for different retained counts. |
| `worst_deleted_donors`, `best_deleted_donors` | Attaining deletion witnesses, encoded as JSON arrays within TSV cells. |
| `worst_retained_count`, `best_retained_count` | Retained counts for those witnesses. |
| `worst_mean`, `worst_sample_variance` | Mean and sample variance at the minimum signed-t witness, rounded for display. They do not separately minimize mean or variance. |
| `robust_positive` | True only when the retained mean is strictly positive for every allowed subset. Zero is not positive. This field is independent of ranking direction. |
| `robust_negative` | True only when the retained mean is strictly negative for every allowed subset; it supports the same sign audit for negative effects. |
| `baseline_topk` | Membership in the all-donor top-K set in the chosen direction. |

`envelopes.tsv` contains one row per gene and exact deletion count, including zero. `minimum_t` and `maximum_t` are signed extrema within that cardinality; `minimum_score` and `maximum_score` display the equivalent fixed-count ratio S/sqrt(Q). The corresponding deletion witnesses attain those extrema. Sum and sum-of-squares columns are exact rational strings for the binary64 inputs. `candidate_count` records the number of scanned two-family candidates, while `total_subsets` is the number of all possible retained subsets. The moment strings and witnesses allow independent reconstruction of the displayed scores.

`audit.json` records the original input file's SHA-256 hash, compression, dimensions, donor order, effective parameters, arithmetic and tie policies, timings, baseline top-K scores, and the shared top-K audit. For gzip input, the hash is over the original compressed bytes, matching the bundled preparation metadata. The top-level `status` matches `topk_audit.status`:

| Status | Meaning |
|---|---|
| `CERTIFIED` | The unordered baseline top-K set is unchanged for every common donor deletion of at most the budget. Ordering within that set can still change. |
| `REFUTED` | `topk_audit.deleted_donors` supplies a common deletion that changes the set; `witness_topk` gives the resulting set in rank order. |
| `UNRESOLVED` | The node limit was reached before the shared-search proof completed. The scalar envelopes are still exact. Increase the limit to pursue a certificate or witness. |
| `NOT_RUN` | `--skip-topk` was selected. The scalar envelopes and baseline top-K are available, without a shared-search conclusion. |

Every gene's extremal witness can use a different deletion subset. Overlapping gene envelopes do not themselves prove that a shared deletion changes the top-K set. A refutation uses one actual common deletion, and a certificate relies on sound bounds throughout the common-subset search.

## Numerical conventions and reproducibility

The paired t display is `mean / (sample_standard_deviation / sqrt(retained_count))`. Nonzero constant effects have signed infinite t; all-zero effects use a documented descriptive t of zero. These extensions support deterministic ranking and are not inferential t tests. Float displays use 17 significant digits. Exact fractions determine extrema and ranks even when a finite score rounds to an infinite display. JSON represents nonfinite displays as strings such as `"inf"` and `"-inf"`.

Gene ranking ties use ascending lexical gene ID. Equal donor values are sorted by input row index; within one cardinality the first attaining two-family candidate is retained. When different cardinalities give equal signed t, the summary keeps the witness with the fewest deletions. Witnesses are deterministic for a fixed input byte sequence and program version; tied extrema can have additional valid witnesses.

Run the independent interface and exact-arithmetic regressions from the source checkout:

```sh
PYTHONPATH=src python -m unittest discover -s tests -v
python benchmarks/independent_theorem_check.py --cases 200000 --output results/independent_theorem_check.json
```

The scalar two-family scan avoids enumerating all subsets. Its cost grows with genes, donors, and the requested deletion counts. The shared top-K search can still require exponentially many nodes in difficult cases; `--max-nodes` and `--skip-topk` expose that limit directly. Benchmark timings are machine- and input-dependent.
