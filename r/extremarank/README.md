# extremarank: native R package

An R package for checking whether omics candidates survive shared biological
sample deletion. The exact kernel is compiled C++ called directly from R.
**No Python runtime, reticulate or command-line bridge is used.**

```r
install.packages(c("Rcpp", "BH", "Matrix", "digest"))
install.packages(
  "https://github.com/heise3/ExtremaRank/releases/download/v0.4.0/extremarank_0.4.0.tar.gz",
  repos = NULL, type = "source")
library(extremarank)

effects <- matrix(c(5, 4, 6, 5, 1, 2, 1, 2, -2, -1, -3, -2), nrow = 4,
  dimnames = list(paste0("donor", 1:4), c("A", "B", "C")))
result <- extremarank_effects(effects, k = 1, budget = 1)
result
result$features[, c("feature_id", "membership_status", "direction_status")]
```

R >= 4.1, Rcpp, BH, Matrix and digest are required. Source installation needs
a C++17 compiler. Windows users need the Rtools version matching their R;
macOS users need the Xcode command line tools. A small standard R source archive
is also attached to the GitHub release; it does not bundle the large benchmark
datasets from the repository.
Alternatively, `remotes::install_github("heise3/ExtremaRank", subdir = "r/extremarank",
ref = "v0.4.0", upgrade = "never")` installs from the tagged repository and
resolves dependencies. The small release archive is preferable if you only need
the R package.

| Entry point | Input and use |
|---|---|
| `extremarank()` | Feature-by-sample matrix or SummarizedExperiment; paired or independent Welch audits |
| `extremarank_effects()` | Donor-by-feature matrix of already computed within-donor effects |
| `extremarank_extrema()` | Exact scalar extrema over retained optional values, with forced included values |
| `refit_extremarank()` | Direct R model refits, covariates and recorded failures |
| `pseudobulk_extremarank()` | Sparse/dense raw cell counts or SingleCellExperiment; donor-level aggregation |
| `pseudobulk_experiment()` | Convert one pseudobulk cell type to SummarizedExperiment |
| `write_extremarank()` | RDS, usable TSV tables, witnesses and model scores |

For matrix input, row names are unique feature IDs and column names are unique
sample IDs. Metadata must match the exact sample ID set. Declare `target` and
`reference`. Paired metadata also needs `donor_id` with one sample per group
and donor. Default preparation is `transform = "none"`; use already prepared
continuous measurements or explicitly select `log2`/`cpm-log2`. Each certificate
applies to the R-prepared finite binary64 values, frozen before deletion.

```r
result <- extremarank(expression, samples, "treated", "control",
  design = "paired", k = 20, budget = 2, transform = "cpm-log2")
result$features
result$audit$topk_minimum_change
result$diagnostics$features
write_extremarank(result, "audit-output")
```

Raw single cells are not independent donors. Supply `cell_id`, `sample_id`,
`donor_id`, `group` and `cell_type` in colData or a data.frame. Assays must be
selected explicitly. Aggregation rejects fractional counts and implicit merges
of separate sample IDs. Missing annotations are errors unless explicitly
excluded with `missing_cell_type = "drop"`; affected cells are recorded.

```r
pb <- pseudobulk_extremarank(sce, assay = "counts", min_cells = 10)
b <- pseudobulk_experiment(pb, "B cells")
result <- extremarank(b, target = "stim", reference = "ctrl", assay = "counts",
  design = "paired", transform = "cpm-log2", k = 20, budget = 2)
```

Install only the optional models or experiment packages you need:

```r
install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "DESeq2", "SingleCellExperiment"))
fits <- refit_extremarank(counts, samples, "treated", "control",
  model = "DESeq2", design = "paired", k = 20, budget = 1)
fits$fit_status
fits$scores$baseline
```

Native audits return `CERTIFIED`, `REFUTED` (actual recomputed witnesses) or
`UNRESOLVED` (valid conservative bounds). Witness indices are **one-based**.
Model refits return `OBSERVED_STABLE`, `OBSERVED_CHANGED`, `PARTIALLY_EVALUATED`
or `NOT_EVALUABLE`; they do not certify unseen deletions. Failed designs and
nonfinite ranking scores remain visible, without clipping or feature removal.
Top-K ranks and effect signs do not establish FDR, biological validity or
clinical usefulness. See `?extremarank` for the complete contract.

Tests run entirely in R, including frozen exact rational expected results:

```sh
R CMD build r/extremarank
R CMD check --no-manual extremarank_0.4.0.tar.gz
```

The independent development oracle is `benchmarks/validate_native_r.py`; Python
is used to verify the package, never to execute the installed R algorithm.
