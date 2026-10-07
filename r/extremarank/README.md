# extremarank: native R package

An R package for checking whether omics candidates survive shared biological
sample deletion. The exact kernel is compiled C++ called directly from R.
**No Python runtime, reticulate or command-line bridge is used.**

```r
install.packages(c("Rcpp", "BH", "Matrix", "digest"))
install.packages(
  "https://github.com/heise3/ExtremaRank/releases/download/v0.6.0/extremarank_0.6.0.tar.gz",
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
ref = "v0.6.0", upgrade = "never")` installs from the tagged repository and
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
R CMD check --no-manual extremarank_0.6.0.tar.gz
```

The independent development oracle is `benchmarks/validate_native_r.py`; Python
is used to verify the package, never to execute the installed R algorithm.


## v0.5 practical controls

`preflight_extremarank()` gives exact feasible deletion counts and approximate
storage estimates. Default LOO diagnostics are skipped at budget zero and
otherwise capped at 100; `max_diagnostics` is independent of search caps.
Use `summary(result)`, `as.data.frame(result)`, `plot(result)` or
`plot(result, type="influence")`. `robustness_curve(effects)` plots native
candidate certificate/refutation/unresolved counts over paired budgets.

Model refits expose per-candidate observed membership/sign changes and actual
witnesses. `coverage` reports planned/evaluated/valid/failed/not-run counts.
No executed deletion returns NOT_EVALUATED; an incomplete stable-looking run
returns PARTIALLY_EVALUATED. `store_scores="disk"` with `score_path` or
`store_scores="none"` avoids holding every score table in the returned object.
Checkpoint/resume applies to model refits and checks frozen input/parameters/
versions; max_refits can be extended. `workers=1:4` uses deterministic scenario
collection and seeded RNG initialization while preserving the caller RNG.

`delete_by="site"` removes whole metadata blocks. A declared named `contrast`
is never silently altered. `refit_extremarank_custom()` supports predeclared
repeated-measure/multi-group models with whole-unit deletion; user callbacks
must return the full frozen feature universe and refit retained samples.
DelayedArray/HDF5Array inputs are optional; cell aggregation is row-blocked and
only its smaller biological-unit matrix is materialized. Native audits also
materialize a sample matrix; max_dense_bytes is a storage guard, not a peak
RAM guarantee. See `vignette("native-workflows", package="extremarank")` and
[中文新功能](../../docs/V05_CN.md).

R 4.6 Windows x86_64 and both macOS architectures have version-specific release
binaries. The R-only installer at `r/install_extremarank.R` verifies SHA-256;
other R versions use the source archive. Restart R before updating a loaded
native package.

## v0.6 Mac pseudobulk

`backend="native"` (also selected by `"auto"`) scans borrowed dense/CSC/CSR/triplet
count views once. `backend="metal"` requests exact integer Metal aggregation on
macOS, with `max_gpu_bytes` bounding requested output and staging buffers.
`extremarank_backend_info()` reports device availability. Other platforms retain
the native CPU path; CUDA is not implemented. Apple unified memory is shared.
The cap excludes R objects, framework allocation and allocator rounding.

```r
extremarank_backend_info()
pb <- pseudobulk_extremarank(counts, cells, backend="native")
# On a Mac with Metal access:
# pb_gpu <- pseudobulk_extremarank(counts, cells, backend="metal", max_gpu_bytes=64*1024^2)
# stopifnot(identical(pb$data, pb_gpu$data))
```

See [the Mac/Windows guide](../../docs/MAC_CN.md) and
[measured validation](../../results/v06/README.md).
