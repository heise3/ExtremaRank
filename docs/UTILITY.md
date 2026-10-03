# Practical use and evidence limits

ExtremaRank answers a concrete sensitivity question: **would the entire selected feature set change if up to `b` biological donor pairs were absent, after recomputing the paired mean and sample variance?** A result is either a complete stability certificate, one shared donor deletion witness, or unresolved at the declared search budget. This can help a researcher assess whether a shortlist is sufficiently stable to justify follow-up work and identify which omitted donor pairs would change it.

The useful deliverable is the actual omitted-pair witness and recomputed feature set. A sensitive list is not automatically incorrect, and an influential donor is not automatically an outlier. Deleting the witness to improve a preferred result would be a different analysis requiring scientific justification.

## Input and ranking contract

Supply a rectangular **donor × feature** matrix of finite paired differences and a unique persistent identifier for each feature. For example, a row can be the difference between target and reference log-expression measurements from one subject. All paired measurements belonging to the same subject are represented by one row, so deletion removes the biological pair together. The software does not infer or repair pair identity from sample labels.

Freeze the contrast direction, paired transformation, normalization, filtering, feature universe, `K`, deletion budget and tie rule before inspecting the audit. The provided real-data preparation uses `log2(CPM+1)` target-minus-reference differences, a declared CPM filter and provenance files. This simple transformation makes a reproducible descriptive example; it is not an assertion that CPM plus paired t is the best RNA-seq inference pipeline.

The API supports three feature ranking directions:

| `direction` | Selected features | Biological use |
| --- | --- | --- |
| `up` | Largest signed refitted paired t | Positive target-minus-reference shortlist |
| `down` | Most negative signed refitted paired t | Negative target-minus-reference shortlist |
| `absolute` | Largest absolute refitted paired t | A two-sided shortlist under this statistic |

The underlying exact signed envelope is converted conservatively for absolute ranking. If a signed interval crosses zero, its lower absolute bound is zero even when no discrete subset attains zero. This can reduce pruning; it cannot support a false certificate. An exact absolute-minimum subset theorem is not claimed.

All genes are ranked on the **same retained donor subset**. At any fixed retained count, the signed score `S/sqrt(Q)` has the same ordering as the refitted paired t under the stated conventions. Exact ties use lexical feature-ID order, including numeric IDs stored as strings. The target is the initial **unordered top-K set**: reordering members within that set is permitted. The API reports the ranked initial and witness sets for readability.

Missing values, infinity, nonrectangular matrices and nonunique feature IDs are rejected. No imputation or gene-specific complete-case deletion is performed. At least two pairs must remain, so `0 <= b <= D-2`. There are no observation weights, covariates, batch corrections, mixed models or moderated variances in this target.

Classical paired t is undefined for a retained all-zero vector and has zero estimated standard error for a nonzero constant vector. This release explicitly extends descriptive rankings: all-zero has score zero; nonzero constants have extended t `+/-infinity`. These conventions are included in the theorem and tests. They do not supply valid p values for singular data.

## Reading the result

| Status | What has been established | Appropriate next action |
| --- | --- | --- |
| `CERTIFIED` | Every shared deletion set of size at most `b` preserves the initial feature set, for the frozen input and specified ranking. | Record the scope and budget with the shortlist. This is sensitivity evidence, not biological validation. |
| `REFUTED` | One listed deletion set changes the set; the witness ranking is recomputed under that common complement. | Inspect the paired metadata and changed features. Retain the witness as an auditable sensitivity result. |
| `UNRESOLVED` | The declared node budget was reached before all remaining possibilities were closed. | Increase the budget if worthwhile, or report the limit. Do not interpret unresolved as stable or unstable. |

`max_nodes` is a node limit, not a wall-clock deadline. `retained_counts_completed` records which exact retained counts have been fully audited before the result. `bound_candidates` counts candidate evaluations; it is not the number of biological samples or independent experiments. The default witness proposals can find counterexamples earlier, but every proposal is verified by recomputing all feature rankings on a common retained subset.

The certificate is exact for the **binary64 numbers supplied to the API**. The implementation forms integer-scaled moments and uses exact sign-aware cross products for decisions. Floating display values can be approximate. Upstream measurement uncertainty, rounding and preprocessing decisions are outside this guarantee.

## What the present checks show

[test_audit.py](../tests/test_audit.py) uses an independent `Fraction` implementation to enumerate shared deletion subsets. It checks complete audit statuses and every refutation witness; it also covers lexical ties, zero data, opposite signs, offsets, cancellation, subnormal values, disparate scales and budget-limited unresolved outcomes. The independent oracle checker in [independent_theorem_check.json](../results/independent_theorem_check.json) completed 200,000 exact small problems with 1,809,171 exhaustive subset evaluations and no failures. These are finite verification results, alongside the proof in [THEOREM.md](THEOREM.md).

[search_ablation.py](../benchmarks/search_ablation.py) compares only the bound in the same search with identical branching, `witness_trials=0` and a 100,000-node limit. All cases are declared deterministic **synthetic** constructions. Ground truth enumerates every common deletion set independently. The six signed-up cases and two explicitly constructed down/absolute cases all agreed with that ground truth for every bound method.

| Fixed synthetic case | Outcome | Two-family nodes | Moment-box nodes | No-bound nodes |
| --- | --- | ---: | ---: | ---: |
| D=12, G=8, K=2, b=2, finite-variance stable winners | Certified | 2 | 30 | 154 |
| D=16, G=12, K=2, b=3, finite-variance stable winners | Certified | 3 | 261 | 1,389 |
| D=12, G=8, K=2, b=2, fixed mixed-sign seed | Refuted | 13 | 13 | 13 |
| D=16, G=12, K=2, b=2, fixed mixed-sign seed | Refuted | 17 | 17 | 17 |
| D=4, G=2, K=1, b=1, variance reversal | Refuted | 4 | 4 | 4 |
| D=4, G=2, K=1, b=2, proportional exact ties | Certified | 18 | 18 | 18 |
| D=12, G=8, K=2, b=2, negative stable winners, down | Certified | 2 | 30 | 154 |
| D=12, G=8, K=2, b=2, negative stable winners, absolute | Certified | 2 | 30 | 154 |

These results isolate a real pruning gain on the constructed stable cases: exact coupled moments shrink a bound that independent numerator/denominator intervals leave wide. The gain is conditional on the data. In easy refutations and correlated ties, envelope computation brings overhead with no node reduction. On these tiny inputs, the no-bound search can be faster even when it visits more nodes. The machine, timings, source constructions, hashes, status equality and witness checks are recorded in [search_ablation.json](../results/search_ablation.json). One local timing run is illustrative and does not establish a general performance claim.

The four-donor example makes the practical target difference explicit. For `a=(2,2.1,1.9,10)` and `h=(1,1.1,0.9,2)`, every donor supports `a>h`. Nevertheless, baseline signed t ranks `h` above `a` (approximately 4.935 versus 2.000). Removing the fourth donor reverses the ranking (approximately 17.321 versus 34.641). Sample variances remain positive in every delete-one subset. A certificate for an additive mean margin would miss this question because the rank change comes from refitting variance.

The proportional example provides the opposite caution. If `a=2h`, every shared retained subset gives exactly tied scores and lexical ID `a` always wins. The independently worst score for `a` can still fall below the independently best score for `h`, attained on different donor subsets. Such overlapping ranges are a reason to continue the common-subset search, not a refutation.

## Practical scope and next actions

The most natural first users have genuine paired small-to-moderate studies and a frozen feature shortlist. The one-feature theorem replaces combinatorial subset enumeration with a linear candidate scan after sorting. The full set audit can still have exponential search cost when gene-wise envelopes overlap. This Python prototype also sorts and encodes conditional values at each node and pays exact-integer bit costs; a node reduction is therefore not automatically an equal runtime reduction. No claim of routine completion for arbitrary genome-wide data or large deletion fractions is justified.

The package contains provenance-backed preparations of GSE87290 (14 complete PBMC pairs, one incomplete pair explicitly excluded) and GSE50760 (18 normal-colon/primary-tumor pairs). Their current real-data results belong to [real_benchmark.json](../results/real/real_benchmark.json); prepared inputs alone are not executed audit results. The benchmark contract, transform, pair map and hashes are retained in each preparation file. A first-100-eligible-gene comparison tests a declared restricted feature universe; it must not be described as a genome-wide top-K certificate.

Useful reproduction steps are to run the complete tests, repeat the synthetic ablation, reproduce paired input preparation, and independently recompute each real witness. A practical extension would make the changed features and donor metadata easier to inspect and benchmark the node oracle on declared larger problems. Covariate-adjusted or moderated targets would require a new derivation; copying this certificate label onto such methods would be unsupported. Publication-level priority also requires closing the full-text gaps recorded in [PRIOR_ART.md](PRIOR_ART.md).

The current evidence supports a concrete exact optimization contribution with an auditable bioinformatics use case. It does not show biological truth, better differential-expression discovery, out-of-sample validity, causal effects, clinical utility or guaranteed GitHub adoption.
