# Automatic model refits and donor-level single-cell inputs

## Refit supported R models

Install R and the relevant Bioconductor package (`limma`, `edgeR` or `DESeq2`)
separately. Core native audits do not require R. Count workflows accept exact
raw nonnegative integer counts, not TPM, normalized abundance or log counts.
The current adapters require each count ≤2^31−1 for consistent R portability.
`limma` accepts already prepared finite continuous expression or abundance.

```bash
extremarank refit data/prepared/Golub/matrix.csv.gz \
  --metadata data/prepared/Golub/metadata.tsv --target AML --reference ALL \
  --model limma --budget 1 --top-k 20 --output limma-audit

extremarank refit data/prepared/KangPB/celltype_1/matrix.tsv.gz \
  --metadata data/prepared/KangPB/celltype_1/metadata.tsv \
  --target stim --reference ctrl --paired --model DESeq2 \
  --min-count-samples 3 --budget 1 --top-k 20 --output deseq-audit

# Real biological lipid data; adjust for the recorded five-level diet factor.
extremarank refit data/prepared/Nutrimouse/matrix.tsv \
  --metadata data/prepared/Nutrimouse/metadata.tsv --target ppar --reference wt \
  --model limma --covariate diet --budget 1 --top-k 5 --output lipid-refits
```

`--covariate NAME` declares a categorical term; `--numeric-covariate NAME`
declares a finite numeric term. Repeat either option as needed. Every retained
scenario fits the full declared design. Paired fits include a donor factor and
remove both samples of every deleted donor. Unpaired fits require one sample
per donor. Repeated visits and technical replicates must be resolved upstream.
Rank-deficient or non-estimable designs remain failed scenarios; covariates,
contrast, genes and donors are never silently changed to obtain a fit.

Filters are fixed using pooled contrast inputs before any deletion: count
sum ≥`--min-total-count` (default 10), positive total, and count ≥`--min-count`
(default 1) in at least `--min-count-samples` samples (default 0). All source
genes define count library totals before filtering. TMM/voom or DESeq2 size
factors and dispersion/moderation are refitted on retained samples. Contrasts
are explicitly target minus reference. All external scores must be finite in
the same frozen feature universe.

| Model | Exported ranking score | Declared fitting behavior |
|---|---|---|
| limma | moderated t | lmFit and eBayes, trend=FALSE, robust=FALSE |
| limma-voom | moderated t | TMM, voom, lmFit, eBayes |
| edgeR | sign(logFC) × sqrt(QL F) | TMM, estimateDisp, glmQLFit robust=FALSE, glmQLFTest |
| DESeq2 | Wald statistic | replacement disabled; Cook's cutoff and independent filtering disabled for complete fixed-universe score export |

These choices are explicit adapter contracts, not universal recommended DE
pipelines. `--max-refits` defaults to 1,000; unrequested feasible subsets remain
outside observed coverage. The exporter retains coefficients, raw statistics,
p values, adjusted p values, fit warnings, design failures and package versions.
The software adds no FDR guarantee. `nonfinite_scores.tsv` preserves invalid
rows; negative QL F is not clipped and genes are not dropped after deletion.

Results are observed `OBSERVED_CHANGED`, `OBSERVED_STABLE`,
`PARTIALLY_EVALUATED`, or `NOT_EVALUABLE`. Even complete coverage of declared
external scenarios is labelled observed. `smallest_observed_membership_change`
is a witness size, not a native certified minimum. `report.html` shows fit
coverage and candidate changes. The existing `compare` command supports custom
externally supplied full-universe refit tables.

## Aggregate single-cell raw counts

Each cell's metadata must contain unique `cell_id` plus `sample_id`, `donor_id`,
`group` and `cell_type`. Count and metadata cell-ID sets must agree exactly.
The aggregation unit is **donor + condition + cell type**. Each deletion unit
is a donor, never an individual cell. `sample_id` identifies the original
library. Multiple libraries within a unit fail unless `--merge-technical`
explicitly declares that they are genuine technical replicates. Do not merge
biological visits or repeated experimental conditions with that option.

```bash
extremarank pseudobulk data/prepared/Kang10X \
  --format 10x --metadata data/prepared/Kang10X/cell_metadata.tsv.gz \
  --min-cells 10 --output my-pseudobulk

# Combined donor-paired count-model audit per usable cell type.
extremarank pseudobulk data/prepared/Kang10X \
  --format 10x --metadata data/prepared/Kang10X/cell_metadata.tsv.gz \
  --audit-design paired --audit-model limma-voom --target stim --reference ctrl \
  --min-count-samples 3 --budget 1 --top-k 20 --output my-sc-audit
```

Dense feature-by-cell CSV/TSV (including gzip) and 10X MTX directories work
without third-party Python libraries. MTX is streamed; duplicate coordinates
sum exactly and only Gene Expression features are retained. Low-cell units
are explicitly excluded and logged. Count sums remain exact integers ≤2^53−1.
Fractional, negative, missing and nonfinite inputs fail.

For H5AD, install `pip install '.[single-cell]'` and explicitly select the raw
count layer:

```bash
extremarank pseudobulk cells.h5ad --format h5ad --counts-layer counts \
  --metadata cells.tsv --min-cells 10 --output my-h5ad-pseudobulk
```

`--counts-layer X` is allowed only when X actually contains raw counts. There
is no fallback to normalized X or `.raw`. Dense, CSR and CSC arrays are read
in bounded cell chunks through backed HDF5 access. IDs and count layers are
checked; cell-level annotations/QC remain the user's responsibility.

Only explicitly declared categorical/numeric covariate columns are carried
forward; they must be constant within each aggregation unit. Other cell
metrics are ignored and recorded. Native audits preserve these fields as
context but do not adjust for covariates. Use automatic R refits to fit them.
Incomplete paired cell types are marked `NOT_EVALUABLE` in the combined
workflow. An explicit manual `study --incomplete-pairs drop` is available if
that policy was selected independently of outcomes.

Every output records input hashes, excluded units, original libraries, declared
raw-count layer, and aggregation choices. Raw genotype identity checks are
`NOT_EVALUATED`. The full Kang example checks 4,347,470 sums against independent
R Matrix row sums, not against a biological ground-truth expression label.
