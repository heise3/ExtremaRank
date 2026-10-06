\name{extremarank}
\alias{extremarank}
\alias{extremarank_effects}
\alias{extremarank_extrema}
\alias{refit_extremarank}
\alias{pseudobulk_extremarank}
\alias{pseudobulk_experiment}
\alias{write_extremarank}
\alias{read_extremarank}
\alias{preflight_extremarank}
\alias{refit_extremarank_custom}
\alias{robustness_curve}
\title{Native R Sample Deletion Robustness for Omics Rankings}
\description{Audit candidate Top-K membership and effect signs under shared biological sample deletion. Native R calls an arbitrary precision C++ kernel directly. No Python interpreter, subprocess bridge or reticulate installation is used. Model refits have a separate observed sensitivity scope.}
\usage{
extremarank(x, metadata = NULL, target, reference, design = c("paired", "welch"), k = NULL,
    budget = 2L, direction = c("up", "down", "absolute"), transform = c("none", "log2",
        "cpm-log2"), features = NULL, all_features = FALSE, max_nodes = 10000L, max_scenarios = 10000L,
    diagnostics = TRUE, assay = NULL, missing = c("error", "drop_features"), incomplete_pairs = c("error",
        "drop"), min_cpm = 0, pseudocount = 1, max_diagnostics = 100L, search_order = c("input",
        "influence"), max_dense_bytes = 512 * 1024^2, block_rows = 1000L)

extremarank_effects(effects, k = NULL, budget = 2L, direction = c("up", "down", "absolute"), features = NULL,
    all_features = FALSE, max_nodes = 10000L, max_scenarios = 10000L, diagnostics = TRUE,
    max_diagnostics = 100L, search_order = c("input", "influence"))

extremarank_extrema(x, retain, included = numeric())

refit_extremarank(x, metadata = NULL, target, reference, model = c("limma", "limma-voom", "edgeR",
    "DESeq2"), design = c("paired", "welch"), k = NULL, budget = 1L, direction = c("up",
    "down", "absolute"), categorical = character(), numeric = character(), min_total = 10,
    min_count = 1, min_count_samples = 0L, max_refits = 1000L, assay = NULL, incomplete_pairs = c("error",
        "drop"), features = NULL, all_features = FALSE, contrast = NULL, delete_by = NULL,
    store_scores = c("memory", "disk", "none"), score_path = NULL, checkpoint = NULL,
    resume = FALSE, progress = FALSE, workers = 1L, max_dense_bytes = 512 * 1024^2, block_rows = 1000L,
    seed = 51006L)

pseudobulk_extremarank(x, metadata = NULL, assay = NULL, min_cells = 10L, covariates = character(),
    merge_technical = FALSE, missing_cell_type = c("error", "drop"), max_dense_bytes = 512 *
        1024^2, block_rows = 1000L)

pseudobulk_experiment(x, cell_type)

write_extremarank(x, path)

read_extremarank(path)

preflight_extremarank(x, metadata = NULL, target, reference, design = c("paired", "welch"), budget = 2L,
    assay = NULL, incomplete_pairs = c("error", "drop"), max_diagnostics = 100L)

refit_extremarank_custom(x, metadata, fit_function, unit_column = "donor_id", k = NULL, budget = 1L,
    direction = c("up", "down", "absolute"), features = NULL, all_features = FALSE, max_refits = 1000L,
    store_scores = c("memory", "disk", "none"), score_path = NULL, checkpoint = NULL,
    resume = FALSE, progress = FALSE, workers = 1L, max_dense_bytes = 512 * 1024^2, block_rows = 1000L,
    seed = 51006L)

robustness_curve(x, budgets = 0:2, ...)
}
\arguments{
\item{x}{A numeric feature-by-sample matrix, numeric Matrix, DelayedMatrix or SummarizedExperiment. Declare its assay explicitly. Direct SingleCellExperiment audit input is rejected: aggregate raw cells first. For pseudobulk, columns are cells and raw nonnegative integer counts are required. For extrema, a finite numeric vector of optional values. For export or pseudobulk conversion, the corresponding result object.}
\item{metadata}{A data.frame aligned by exact ID sets. Audit/refit columns: sample_id, group, and donor_id (required for paired designs). Meaningful row names can supply sample_id. Cell columns: cell_id, sample_id, donor_id, group, cell_type. colData is used when omitted for experiment objects. Factors and leading-zero identifiers are preserved as strings.}
\item{target, reference}{Distinct, explicit group identifiers. Other groups are excluded and recorded.}
\item{design}{Paired within-donor target-minus-reference effects, or independent two-group Welch statistics. Repeated donor IDs are rejected in Welch designs.}
\item{k}{Top-K size; NULL selects the smaller of 20 and the frozen feature count. Ties use UTF-8 byte ordering of feature IDs.}
\item{budget}{Maximum deleted biological units; at least two observations per native group must remain. Paired deletions remove whole donors.}
\item{direction}{Descending signed statistics, ascending signed statistics, or descending absolute statistics. Effect signs always refer to target minus reference.}
\item{transform}{Frozen preparation: none, log2(x + pseudocount), or log2(CPM + pseudocount). CPM libraries use all source genes before the frozen feature filter. Preparation is never recomputed by the native audit.}
\item{features}{Specific feature IDs to audit for membership and effect signs. NULL audits the original Top-K. All source features still compete in every ranking.}
\item{all_features}{Audit every retained feature. Cannot be combined with features.}
\item{max_nodes, max_scenarios}{Positive global limits on search nodes and evaluated deletion scenarios. Hitting a cap preserves valid bounds and can return UNRESOLVED. Leave-one-out diagnostics run separately and do not consume these search caps.}
\item{diagnostics}{TRUE or "auto" computes bounded one-unit diagnostics when budget is positive. FALSE or "none" disables them; "all" explicitly enables them even at budget zero. max_diagnostics limits their count. Fractions concern only evaluated deletions and are descriptive.}
\item{assay}{Explicit assay name for an experiment object; NULL for an input matrix. Select raw counts for pseudobulk.}
\item{missing}{Reject missing values, or explicitly remove affected features before freezing preparation. CPM always requires complete counts. Infinite values are rejected.}
\item{incomplete_pairs}{Reject incomplete donor pairs, or explicitly exclude and record them. Duplicate donor/group samples are rejected.}
\item{min_cpm}{For CPM preparation, retain genes meeting this threshold in at least half the included samples.}
\item{pseudocount}{Finite positive value for log transforms.}
\item{effects}{Numeric donor-by-feature matrix of already computed within-donor effects. Row and column IDs are required. This orientation matches the original paired-effects interface.}
\item{retain}{Number of optional values to retain, in addition to every included value. At least two values must remain.}
\item{included}{Finite values forced into every candidate retained subset.}
\item{model}{Established R model to fit directly. limma accepts prepared continuous values; the other models accept raw integer counts.}
\item{categorical, numeric}{Explicit covariate column names. Covariates, filters and the declared design remain fixed after deletion. Normalization and variance or dispersion estimation are recomputed. Rank deficient designs, missing residual degrees of freedom, and nonfinite scores are recorded as NOT_EVALUABLE.}
\item{min_total, min_count, min_count_samples}{Count filter frozen on the baseline included libraries: total count threshold and a minimum number of samples meeting min_count. Not applied to continuous limma inputs.}
\item{max_refits}{Maximum deletion refits, excluding the baseline. Enumeration is streamed in ascending deletion cardinality. A cap is recorded even when all observed refits are stable.}
\item{min_cells}{Minimum cells per donor/group/cell-type pseudobulk unit. Excluded units are retained in the result.}
\item{covariates}{Cell metadata columns to carry into pseudobulk; each must be constant within a unit.}
\item{merge_technical}{Explicitly permit multiple sample IDs in one donor/group/cell-type unit. Use only for technical replicates; separate biological visits must not be merged implicitly.}
\item{missing_cell_type}{Reject missing annotations, or explicitly exclude and record the affected cell IDs. Raw counts are validated before exclusion.}
\item{cell_type}{An existing name in the pseudobulk result's data list.}
\item{max_diagnostics}{Nonnegative independent limit on LOO diagnostics; default 100. Unrun units are recorded.}
\item{search_order}{Input or influence-guided branch and witness order. Heuristics do not change exact proof comparisons.}
\item{max_dense_bytes}{Maximum bytes of a dense prepared matrix or pseudobulk output, default 512 MiB. Model worker copies are also bounded. Excludes model objects and exact-integer storage.}
\item{block_rows}{Maximum rows per input/aggregation block; automatically reduced to keep input blocks near 64 MiB. DelayedArray and HDF5Array are optional input backends. The smaller biological-unit output is materialized.}
\item{contrast}{Optional finite named coefficient vector. Model-matrix coefficient names are used (grouptarget, cv1, cv2, donor indicators). Treatment coding is explicit and independent of global contrasts options. Missing nonzero coefficients or a changed reference level for a contrast term make that fit unevaluable. NULL tests target minus reference.}
\item{delete_by}{Optional metadata column for whole-block deletion, such as site. A donor cannot span deletion blocks. Budget counts blocks; native two-observation constraints remain.}
\item{store_scores}{Keep scores in memory, save each to disk, or retain only summaries. Default memory.}
\item{score_path}{New/empty directory for disk-backed score RDS files. Its contract is checked on resume.}
\item{checkpoint}{Optional new/empty directory for atomic per-scenario fit records and an input/parameter/version fingerprint. Applies to model refits.}
\item{resume}{Reuse matching completed model fits from checkpoint; permits increasing max_refits. Any frozen input, parameter or built-in package-version change is rejected.}
\item{progress}{Emit model-refit progress messages.}
\item{seed}{Fixed integer seed used to give each scenario a reproducible L'Ecuyer-CMRG RNG initialization. The caller RNG state is restored. Included in checkpoint contracts.}
\item{workers}{One to four PSOCK processes; results are collected in deterministic scenario order. More workers can increase RAM.}
\item{fit_function}{Custom function(x, metadata) returning each frozen feature exactly once in a data.frame with feature_id, numeric score and effect. It must refit on the supplied retained samples; use for explicitly declared repeated-measure or multi-group models. Control its external dependencies and package versions separately.}
\item{unit_column}{Metadata column identifying whole deletion units for custom refits. At least two units remain.}
\item{budgets}{Nonnegative integer donor-deletion budgets for a curve of precomputed paired effects.}
\item{...}{Arguments forwarded to extremarank_effects by robustness_curve.}
\item{path}{New or empty output directory for a result RDS, flat TSV tables, witness records and model scores.}
}
\details{
The paired conditional extrema algorithm evaluates sorted contiguous blocks and prefix/suffix endpoint subsets. Decisions use exact dyadic integer moments of the prepared finite IEEE-754 binary64 values. Welch bounds are conservative rational enclosures; joint rankings can remain unresolved even if individual scalar endpoints are exact. Zero mean scores are descriptively zero; nonzero zero-variance scores are signed infinity. Displayed t statistics use floating point; ranking and proof comparisons do not.

Native CERTIFIED means every allowed deletion up to the declared budget preserves the queried property. REFUTED requires a fully recomputed actual witness. UNRESOLVED preserves conservative minimum-change intervals. All returned witness indices are one-based. For a stable property the lower bound budget + 1 means no change was found or possible within the budget, not that a change at budget + 1 exists.

Overall status describes whole Top-K membership; per-candidate direction statuses are separate. Observed rank extrema include the baseline; selection fractions concern valid executed deletion refits and are not probabilities.

Model refits return NOT_EVALUATED if no deletion was executed. Capped or failed stable-looking runs return PARTIALLY_EVALUATED. Complete evaluable unchanged runs return OBSERVED_STABLE; actual changes return OBSERVED_CHANGED. A failed baseline returns NOT_EVALUABLE. They never claim an exact certificate for unseen deletions. Failed scenarios and their invalid feature IDs remain available. Raw scores are retained when memory or disk storage is selected; none explicitly discards them. The ranking universe is not repaired after observing failures. This package audits ranking sensitivity; it does not establish differential-expression validity, FDR, causality or clinical utility.

Pseudobulk aggregation accepts sparse matrices without making a dense cell-by-gene copy. Nonnegative integer inputs and every aggregated entry must lie within the exact binary64 integer range. The smaller feature-by-unit output is materialized for subsequent analysis.
}
\value{
Audit and refit functions return an extremarank_result list. Native results contain features (a data.frame with per-property statuses, bounds and list-column witnesses), audit, diagnostics, preparation and elapsed_seconds. Refit results contain fit_status, scores, scenarios, versions, enumeration_capped and minimum_observed_change. Printing gives a compact summary. Model results include per-candidate observed changes/witnesses, exact planned coverage, disk score paths and checkpoint contracts. preflight_extremarank returns deletion-count and approximate storage estimates without fitting. robustness_curve returns budget-by-status counts and full results.

extremarank_extrema returns exact minimum and maximum descriptions, optional witness indices, exact moment strings, candidate counts and total subset counts. pseudobulk_extremarank returns an extremarank_pseudobulk with named data entries (counts and metadata), excluded_units and provenance. pseudobulk_experiment returns a SummarizedExperiment. write_extremarank returns the normalized directory invisibly. Disk scores are copied into a self-contained export with checksums. read_extremarank reloads and validates the disk paths relative to the moved export directory.
}
\examples{
effects <- matrix(c(5, 4, 6, 5, 1, 2, 1, 2, -2, -1, -3, -2),
                  nrow = 4, dimnames = list(paste0("donor", 1:4), c("A", "B", "C")))
result <- extremarank_effects(effects, k = 1, budget = 1)
result
result$features[, c("feature_id", "membership_status", "direction_status")]
extremarank_extrema(c(-3, -1, 2, 4, 6), retain = 3, included = 1)

counts <- matrix(1:24, nrow = 3,
                 dimnames = list(c("A", "B", "C"), paste0("cell", 1:8)))
cells <- data.frame(cell_id = colnames(counts),
  sample_id = rep(paste0("sample", 1:4), each = 2),
  donor_id = rep(c("01", "01", "02", "02"), each = 2),
  group = rep(c("control", "treated", "control", "treated"), each = 2),
  cell_type = "B cells")
pb <- pseudobulk_extremarank(counts, cells, min_cells = 2)
pb$data[["B cells"]]$counts
}
\seealso{\url{https://github.com/heise3/ExtremaRank}}
\keyword{models}
