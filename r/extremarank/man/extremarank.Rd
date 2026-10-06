\name{extremarank}
\alias{extremarank}
\alias{extremarank_effects}
\alias{extremarank_extrema}
\alias{refit_extremarank}
\alias{pseudobulk_extremarank}
\alias{pseudobulk_experiment}
\alias{write_extremarank}
\title{Native R Sample Deletion Robustness for Omics Rankings}
\description{Audit candidate Top-K membership and effect signs under shared biological sample deletion. Native R calls an arbitrary precision C++ kernel directly. No Python interpreter, subprocess bridge or reticulate installation is used. Model refits have a separate observed sensitivity scope.}
\usage{
extremarank(x, metadata = NULL, target, reference,
  design = c("paired", "welch"), k = NULL, budget = 2L,
  direction = c("up", "down", "absolute"),
  transform = c("none", "log2", "cpm-log2"), features = NULL,
  all_features = FALSE, max_nodes = 10000L, max_scenarios = 10000L,
  diagnostics = TRUE, assay = NULL, missing = c("error", "drop_features"),
  incomplete_pairs = c("error", "drop"), min_cpm = 0, pseudocount = 1)

extremarank_effects(effects, k = NULL, budget = 2L,
  direction = c("up", "down", "absolute"), features = NULL,
  all_features = FALSE, max_nodes = 10000L, max_scenarios = 10000L,
  diagnostics = TRUE)

extremarank_extrema(x, retain, included = numeric())

refit_extremarank(x, metadata = NULL, target, reference,
  model = c("limma", "limma-voom", "edgeR", "DESeq2"),
  design = c("paired", "welch"), k = NULL, budget = 1L,
  direction = c("up", "down", "absolute"), categorical = character(),
  numeric = character(), min_total = 10, min_count = 1,
  min_count_samples = 0L, max_refits = 1000L, assay = NULL,
  incomplete_pairs = c("error", "drop"))

pseudobulk_extremarank(x, metadata = NULL, assay = NULL,
  min_cells = 10L, covariates = character(), merge_technical = FALSE,
  missing_cell_type = c("error", "drop"))

pseudobulk_experiment(x, cell_type)

write_extremarank(x, path)
}
\arguments{
\item{x}{A numeric feature-by-sample matrix, numeric Matrix, or SummarizedExperiment. Declare its assay explicitly. Direct SingleCellExperiment audit input is rejected: aggregate raw cells first. For pseudobulk, columns are cells and raw nonnegative integer counts are required. For extrema, a finite numeric vector of optional values. For export or pseudobulk conversion, the corresponding result object.}
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
\item{diagnostics}{Also compute all feasible one-unit deletion rankings. Their empirical selection fractions are descriptive, not probabilities or multi-deletion certificates.}
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
\item{path}{New or empty output directory for a result RDS, flat TSV tables, witness records and model scores.}
}
\details{
The paired conditional extrema algorithm evaluates sorted contiguous blocks and prefix/suffix endpoint subsets. Decisions use exact dyadic integer moments of the prepared finite IEEE-754 binary64 values. Welch bounds are conservative rational enclosures; joint rankings can remain unresolved even if individual scalar endpoints are exact. Zero mean scores are descriptively zero; nonzero zero-variance scores are signed infinity. Displayed t statistics use floating point; ranking and proof comparisons do not.

Native CERTIFIED means every allowed deletion up to the declared budget preserves the queried property. REFUTED requires a fully recomputed actual witness. UNRESOLVED preserves conservative minimum-change intervals. All returned witness indices are one-based. For a stable property the lower bound budget + 1 means no change was found or possible within the budget, not that a change at budget + 1 exists.

Model refits return OBSERVED_STABLE, OBSERVED_CHANGED, PARTIALLY_EVALUATED or NOT_EVALUABLE. They never claim an exact certificate for unseen deletions. Every failed scenario and any nonfinite raw scores remain available. The ranking universe is not repaired after observing failures. This package audits ranking sensitivity; it does not establish differential-expression validity, FDR, causality or clinical utility.

Pseudobulk aggregation accepts sparse matrices without making a dense cell-by-gene copy. Nonnegative integer inputs and every aggregated entry must lie within the exact binary64 integer range. The smaller feature-by-unit output is materialized for subsequent analysis.
}
\value{
Audit and refit functions return an extremarank_result list. Native results contain features (a data.frame with per-property statuses, bounds and list-column witnesses), audit, diagnostics, preparation and elapsed_seconds. Refit results contain fit_status, scores, scenarios, versions, enumeration_capped and minimum_observed_change. Printing gives a compact summary.

extremarank_extrema returns exact minimum and maximum descriptions, optional witness indices, exact moment strings, candidate counts and total subset counts. pseudobulk_extremarank returns an extremarank_pseudobulk with named data entries (counts and metadata), excluded_units and provenance. pseudobulk_experiment returns a SummarizedExperiment. write_extremarank returns the normalized directory invisibly.
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
