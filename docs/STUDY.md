# Paired, independent-group, and external-model workflows (v0.3)

ExtremaRank offers native and external ways to inspect candidate-list sensitivity. Choose
the design and the ranking statistic deliberately; a shared file format does
not make the statistical models interchangeable.

| Workflow | Target | Guarantee and scope |
|---|---|---|
| `study --design paired` or the original effects CLI | Ordinary paired t on frozen within-donor differences | Exact scalar extrema via the two-family theorem; shared Top-K branch-and-bound may be unresolved at its limit. |
| `study --design welch` | Ordinary unequal-variance two-group Welch statistic on independent samples | Exact integer ranking comparisons; safe conditional pruning or reference enumeration. Does not inherit the paired extremum theorem. |
| `compare` | Scores exported from any external refit method | Only the supplied scenarios are inspected. `OBSERVED_STABLE` / `OBSERVED_CHANGED` never certify unobserved deletions or independently verify the external fits. |

Real validation includes paired bulk RNA, the complete Golub processed microarray,
Kang donor-level single-cell pseudobulk, technical CPTAC protein measurements,
and biological Nutrimouse lipid percentages. Automatic refits cover limma,
limma-voom, edgeR and DESeq2; some declared fits can be non-evaluable and remain
visible. See [v0.3 results and scope](V03_CN.md),
[individual certificates](ROBUSTNESS.md), and [refit/single-cell inputs](REFITS.md).

## Matrix and metadata

The default matrix orientation is **features in rows, samples in columns**:

```text
feature_id  S1  S2  S3  S4
GENE_A      4   5   1   2
GENE_B      1   2   4   5
```

Use actual tabs in TSV, or commas in CSV. Files may be gzip compressed. For
samples in rows, use `--orientation samples` and first header `sample_id`.
Feature and sample IDs must be unique and nonempty. Metadata has columns:

```text
sample_id  group    donor_id
S1         treated  D1
S2         treated  D2
S3         control  D3
S4         control  D4
```

Matrix and metadata sample ID sets must match exactly. Metadata is reordered
to the matrix; no fuzzy renaming or silent intersection is performed. Other
groups are excluded only by the explicit target/reference contrast and their
IDs are recorded. Contrast direction is always **target minus reference**.

Paired designs require `donor_id` and exactly one observation per donor per
contrast group. The whole donor pair is deleted together. Incomplete pairs
cause an error by default. `--incomplete-pairs drop` explicitly excludes them
and records their IDs. Aggregate technical replicates upstream with a suitable
method; duplicate donor/group rows are rejected.

Welch designs require one independent donor per sample across both groups.
`donor_id` is optional and defaults to `sample_id`; supplying it allows the
tool to reject repeated donors. Declaring unique IDs is not evidence that the
biological observations are independent or correctly labeled. Genotype checks,
cell-to-donor identity, batch/sex/age adjustment and repeated-measure modeling
remain upstream responsibilities. Additional metadata columns are retained as
context and are **not fitted as covariates**.

## Run a complete study

```bash
# Independent groups, processed expression or abundance matrix.
extremarank study data/prepared/Golub/matrix.csv.gz \
  --metadata data/prepared/Golub/metadata.tsv \
  --target AML --reference ALL --design welch \
  --budget 2 --top-k 20 --direction up --output golub-audit

# Genuine paired samples; normalized measurements are supplied unchanged.
extremarank study expression.tsv --metadata samples.tsv \
  --target treated --reference control --design paired \
  --budget 2 --top-k 20 --output paired-audit

# Explicit descriptive raw-count transformation and frozen filter.
extremarank study counts.tsv --metadata samples.tsv \
  --target treated --reference control --design paired \
  --transform cpm-log2 --min-cpm 1 --output count-audit
```

`--transform none` is the default for already prepared numeric measurements;
it permits negative processed values. `log2` computes log2(value+pseudocount).
`cpm-log2` uses each included sample's sum over **all source features**, filters
by CPM in at least half the included samples if `--min-cpm` is positive, and
computes log2(CPM+pseudocount). The default pseudocount is 1. Both log transforms
require nonnegative values; CPM requires positive finite library totals.
This option is a declared descriptive preparation, not a replacement for a
domain-specific normalization or RNA-seq significance pipeline. Preparation,
filtering and normalization are frozen before any deletion.

Missing values fail by default. For an already processed abundance matrix,
`--missing drop-features` explicitly removes every feature missing in any
included sample and lists the exclusions. No values are imputed and no
gene-specific donor subsets are used. This option is unavailable for CPM
preparation because missing counts would make library totals ambiguous.
Dropping features changes the universe and can introduce selection bias; choose
the policy before inspecting robustness. Infinity is always rejected.

`prepare` accepts the same preparation options and writes inputs and
`preparation.json` without auditing:

```bash
extremarank prepare expression.tsv --metadata samples.tsv \
  --target treated --reference control --design paired --output prepared
extremarank prepared/effects.csv --budget 2 --top-k 20 --output exact-envelopes
```

## Search limits and outputs

`--budget` is the **total** number of deleted biological units, not a separate
budget in each group. In Welch mode, only deletions retaining at least two
samples per group are feasible; the maximum budget is N-4. Paired mode requires
at least two retained pairs. These minima are mathematical requirements, not
recommended study sample sizes. `--max-nodes` bounds native branch-and-bound searches;
`--max-subsets` bounds the Welch nonbaseline deletion sets (default 10,000).
Neither is a wall-clock timeout. v0.3 defaults to safe conditional Welch pruning;
`--welch-search enumeration` retains the reference path. Search can still grow
combinatorially. See [individual certificates and minimum changes](ROBUSTNESS.md).
Automatic R refits and single-cell inputs are described in [REFITS.md](REFITS.md).

Welch uses

`t = (mean_target - mean_reference) / sqrt(s_target²/n_target + s_reference²/n_reference)`.

The [NIST definition](https://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm)
specifies this unequal-variance statistic. Decisions use integer dyadic
moments and sign-aware comparisons of squared rational statistics; no square
root or rounded display is used in a rank decision. If the means agree,
descriptive t is zero even with zero standard error. Unequal constant groups
have signed infinite t. These extensions do not supply inferential p values.

For n target and m reference integer-scaled observations, write sums S_t/S_r
and square sums Q_t/Q_r, Δ=m*S_t-n*S_r, V_t=n*Q_t-S_t²,
V_r=m*Q_r-S_r². The exact squared magnitude is

`Δ²*(n-1)*(m-1) / [V_t*m²*(m-1) + V_r*n²*(n-1)]`.

The common scale cancels within each feature. Cross products compare squared
magnitudes across features; sign is handled before squaring. This is ordinary
algebra and enumeration, not a newly claimed Welch extremum theorem.

`audit.json` records `CERTIFIED`, `REFUTED` with one actual common deletion,
or `UNRESOLVED` at the limit. Certificates concern the unordered original set,
not its internal order. A refuted shortlist is not biologically incorrect and
an influential donor is not automatically an outlier to remove.

Every `study` run also inspects **all feasible single-unit deletions**, unless
`--no-diagnostics` is selected:

- `report.html`: offline status, actual witness, changed candidates per sample,
  initial candidates and the 20 outsiders with the largest inspected rank span.
- `sample_influence.tsv`: each deleted unit, group, Jaccard, exits and entries.
- `diagnostics.tsv/json`: every feature's baseline, best/worst rank among
  baseline plus delete-one scenarios, selection count and sign changes.
- `effects.csv` (paired) or `values.csv` plus `groups.tsv` (Welch): round-trip
  prepared inputs; `preparation.json` records hashes, exclusions and metadata.
- `baseline.tsv` (Welch): signed t display, rank and target-reference mean
  difference. Paired runs retain `summary.tsv` and exact `envelopes.tsv`.

Diagnostic selection fractions are fractions of inspected deletion scenarios,
not statistical probabilities. They do not represent bounds for two or more
deletions, even if `audit.json` used a larger budget. With too few retained
observations, invalid delete-one units are listed as skipped. When no delete-one
scenario is feasible, selection fractions are null. Use a fresh output directory
per configuration; result filenames are overwritten on reruns.

## Connect external statistical models

For limma, DESeq2, covariate models or custom scores, refit the appropriate
method after deleting each biological unit and export **the full fixed feature
universe**, not only its significant hits. `compare` consumes three tables:

| File | Required columns |
|---|---|
| Baseline | `feature_id`, `score` |
| Refits (long format) | `scenario_id`, `feature_id`, `score` |
| Deletion manifest | `scenario_id`, `deleted_units` (quoted JSON array of unique unit IDs) |

All scores must be finite and every scenario must contain the same features
as baseline. The comparator rejects missing values and changed universes.
`up` chooses larger scores, `down` smaller scores, `absolute` larger magnitudes.
Ties use ascending lexical feature ID. Choose the direction according to the
actual exported statistic; `absolute` is appropriate for signed t/z statistics,
while raw p-value ranking would use `down` and is not an FDR audit.

```bash
extremarank compare data/prepared/Golub/limma/baseline.tsv \
  data/prepared/Golub/limma/refits.tsv.gz \
  --manifest data/prepared/Golub/limma/manifest.tsv \
  --method 'limma refitted moderated t; AML minus ALL' \
  --direction absolute --top-k 20 --output limma-audit
```

The result is `OBSERVED_STABLE` or `OBSERVED_CHANGED`, with inspected scenario
changes and feature rank spans. Supplied method names and deletion manifests
are provenance provided by the caller, not independently validated model fits.
This interface can ingest DESeq2/custom refit outputs, but the included executed
adapter is limma. It is **not** a built-in DESeq2/GLM/mixed-model fitter, and
does not apply the paired theorem to those models.

The included [limma exporter](../benchmarks/limma_refits.R) recomputes empirical
Bayes moderation at each of 38 single-sample deletions, with all 3,051 features
and declared AML-minus-ALL contrast. Regeneration requires R and limma;
comparison of its bundled tables uses only Python. The statistical workflow
follows the [limma author guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf).

## Reproduce the extension

```bash
python -m unittest discover -s tests -v
python benchmarks/independent_welch_check.py --cases 2000
python benchmarks/upgrade_real.py

# Optional: regenerate original public processed matrix and external fits.
Rscript benchmarks/prepare_golub.R data/source/Golub/golub.RData scratch-golub
python benchmarks/verify_golub_roundtrip.py scratch-golub
Rscript benchmarks/limma_refits.R scratch-golub/matrix.csv scratch-golub/metadata.tsv scratch-limma
```

The [frozen execution contract](../benchmarks/upgrade_contract.json) was written
before these real-data audits. Provider preprocessing is retained verbatim in
`data/source/Golub/provider_preprocessing.R`; it includes thresholds, filtering,
log10 transformation and column scaling. The prepared data are already
processed and do not support interpretation of their mean difference as raw
log fold change. Sample IDs are deterministic provider-column labels, not newly
verified patient identities. See [data attribution](../data/DATA_LICENSE.md).
