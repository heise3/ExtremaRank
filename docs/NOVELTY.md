# ExtremaRank: technical contribution and claim ceiling

ExtremaRank's intended contribution is an exact conditional extremum algorithm for the signed self-normalized donor score that ranks features identically to refitted ordinary paired t. It replaces enumeration of all fixed-cardinality donor subsets for a single feature with two explicit ordered subset families, and supplies that oracle to a search over shared donor deletions. The contribution is an algorithmic statement that can be proved and falsified, rather than a new name for robustness auditing.

## The specific algorithmic advance

At a search node, some donors are already retained, some are already deleted, and the remainder are undecided. For each feature, fixed retained donors contribute first and second moments `(a,b)`. Choosing `m` of `q` undecided donors changes both the numerator and the normalization:

\[
f(A) = \frac{a+\sum_{i\in A}x_i}{\sqrt{b+\sum_{i\in A}x_i^2}},\qquad |A|=m.
\]

The two-family theorem states that a minimizing and a maximizing subset can be found in the union of:

1. Every contiguous length-`m` block in the value-sorted undecided donors.
2. Every `j` smallest plus `m-j` largest subset, for `j=0,...,m`.

For `1 <= m <= q`, the union contains at most `(q-m+1)+(m+1)=q+2` candidates before removing duplicates. Prefix sums of values and squared values evaluate each candidate's moments in constant arithmetic operations. A fixed-`m` oracle therefore scans `O(q)` candidates after sorting, instead of `binomial(q,m)` subsets. The edge cases `m=0` and `m=q` have only one distinct subset.

At the root, scanning deletion counts `0,...,budget` gives `O(D*(budget+1))` candidate evaluations per feature after `O(D log D)` sorting. For `G` features the root cost is `O(G D log D + G D (budget+1))` arithmetic operations. This is not a polynomial-time claim for the entire exact topK search: coupled rankings may still require exponentially many search nodes. Exact rational arithmetic also has bit cost that an arithmetic-operation count does not capture.

The conditional `(a,b)` version is essential. A root-only trimmed scan cannot automatically give valid bounds after a search fixes donors. An exact node oracle permits safe pruning at every node and preserves completeness when the shared search is allowed to finish.

## Mathematical reason for the two families

For positive extrema, the ratio identity

\[
S^2/Q=\max_{\lambda>0}\{2\lambda S-\lambda^2 Q\},\qquad S>0,\ Q>0,
\]

exposes a positive maximum through donor scores `2 lambda x - lambda^2 x^2`. These form a concave quadratic in the donor value, so a maximizing selection can be contiguous around its center. The auxiliary-variable identity is established fractional-programming algebra; see [Shen and Yu's quadratic transform](https://www.comm.utoronto.ca/~weiyu/FP_part1.pdf). ExtremaRank does not claim to invent it.

The opposite family follows from an exchange argument for a positive minimum when every feasible numerator is nonnegative. With every other retained value fixed, the score as a function of one selected value is `(A+x)/sqrt(B+x^2)`. Its derivative has numerator `B-A*x`; an interior stationary point with positive score is a maximum. An interior selected value can therefore be replaced by an outer unselected endpoint without increasing the score. Repeating this exchange grows the retained prefix or suffix until an endpoint-family selection is reached. Negation handles the negative extrema; zero extrema require the explicit sign/denominator edge cases in the proof.

The complete proof in [THEOREM.md](THEOREM.md) handles zero denominators, zero sums, duplicates, fixed included donors and the statistic's degeneracy policy. Randomized exhaustive checks support an implementation; they do not replace the proof.

Contiguous scans already solve the univariate minimum-variance subset problem, as documented in the [MCD authors' review](https://wis.kuleuven.be/statdatascience/robust/papers/2010/wire-mcd.pdf). ExtremaRank claims neither that ordered windows are new nor that the convexity ingredient is new. The screened technical content is the complete signed-ratio characterization with both families and arbitrary forced retained values.

## Why this is beyond the additive DonorCert prototype

An additive signed mean ranking has a sorted-prefix worst-deletion margin. A paired t ranking additionally refits sample variance, so removing a donor can change a feature's normalization and reverse its rank even if all donors support the same mean ordering. The previous additive algorithm does not certify this target.

For example, take paired differences `a=(2,2.1,1.9,10)` and `h=(1,1.1,0.9,2)`. Every row has `a > h`, so mean ordering is preserved under any common nonempty deletion. Baseline signed t is approximately `2.000` for `a` and `4.935` for `h`. After deleting donor 4, it is approximately `34.641` for `a` and `17.321` for `h`. All delete-one subsets in this example have positive sample variance. This is a change in a variance-refitted target, rather than a mean audit presented under another label.

The oracle's exact target is the signed score computed from frozen donor-feature differences. It supports audits of up, down, and absolute paired-t rankings. For the absolute-ranking audit, converting a signed envelope that crosses zero to a lower magnitude bound of zero is conservative: that zero need not be attained by a feasible discrete subset. No exact absolute-minimum oracle is claimed. This is a ranking audit; it does not establish a new differential-expression significance test, moderated t method, negative-binomial fit, FDR control, or a certificate for preprocessing refitted after deletion.

## Established components and the publication ceiling

The [prior-art screen](PRIOR_ART.md) documents the evidence and retrieval limits. In particular, earlier work already provides ordered subset statistics for gene expression, optimization over trimming proportions, data-removal audits in differential expression, exact deletion optimization for regression, and the topK boundary-pair reduction. Those components must be credited.

The project's precise research claim is:

> We develop and implement a conditional two-family extremum oracle for signed self-normalized donor scores, and use its exact bounds to accelerate a complete shared-donor audit of paired-t topK stability.

This formulation says what the method contributes. It makes no global mathematical-priority assertion. The targeted screen found no primary source proving the same conditional two-family characterization, but absence from a search is not proof of novelty. The closest adaptive-trimming paper's full text remains unverified; author software documentation verifies a different trim-grid procedure. See [PRIOR_ART.md](PRIOR_ART.md) for the exact distinction.

## Certification and numerical claims

One-feature extrema may be attained by different donor subsets for different features. Separating their ranges is a valid sufficient certificate. Overlapping their ranges is not a counterexample and does not show an attainable worst rank. The shared-subset search must close the remaining branches or return an unresolved result. A refutation must name one shared deletion set and independently recompute the complete topK ranking on its complement.

The mathematical theorem concerns exact arithmetic. The implementation lifts each finite binary64 input to an exact integer on a shared power-of-two scale, forms exact first and second moments, and decides extrema and score order by integer cross products. Floating square roots are used for display only; they do not determine an audit result. Input conversion freezes the supplied binary64 values, so the certificate does not cover uncertainty or rounding in upstream measurement or preprocessing. There is no faster floating decision path that would require a separate enclosure proof. A passing finite test suite remains evidence about tested cases, rather than a proof of program correctness.

Ordinary paired t is undefined for all-zero retained values and has zero estimated standard error for a nonzero constant retained vector. The implementation uses a documented descriptive extension: the all-zero vector has score zero; a nonzero constant vector has signed self-normalized score `+/-sqrt(n)` and extended t `+/-infinity`. The theorem includes those cases. These values are ranking conventions and do not create a valid classical t test for a singular subset. No stabilizing constant or silent subset dropping is applied.

## Evidence that makes the contribution worth reproducing

The algorithmic contribution becomes reviewable through four discriminating checks:

1. **Exact oracle correctness.** Exhaustively enumerate small rational problems with fixed included donors and every feasible `m`. Verify both extrema and attaining subsets. Cover mixed signs, zero sums, duplicates, unequal scales and degeneracies under the chosen target rule.
2. **Candidate-family necessity.** Include examples where an interior block is needed and where a disconnected prefix-plus-suffix selection is needed. Compare the full oracle with block-only, endpoint-only and independent moment-box bounds.
3. **Search gain.** Use the same shared-donor search, branch ordering and stopping budget for exact envelopes and moment boxes. Report certified/refuted/unresolved outcomes, visited nodes, runtime, memory and verified witnesses. This isolates the oracle's contribution from implementation speed or a lucky early counterexample.
4. **Biological utility.** Audit real paired transcriptomic differences with provenance and frozen preprocessing. Report which initial gene sets are stable and which have actual deletion witnesses. Separate that audit from any claim of biological truth, population generalization or improved DE significance.

Small genuine paired studies provide exact validation and an interpretable demonstration. Larger declared simulations can demonstrate computation beyond exhaustive enumeration. Copies of the same biological donors are not evidence from additional independent donors. A successful biological example alone does not establish the computational advance; a synthetic theorem test alone does not establish biological usefulness.

The intended result is therefore a concrete algorithm with a polynomial-size one-feature candidate family and a transparent exact/anytime shared-subset audit. Publication-level claims require the complete theorem, numerical correctness appropriate to the implementation, primary-method collision checks, and evidence of measurable gain. GitHub visibility and star counts cannot be guaranteed by these results.

The current synthetic search ablation in [search_ablation.json](../results/search_ablation.json) holds branch order, witness trials and stopping limits fixed. Its two stable finite-variance constructions show fewer nodes with the exact envelope than with independent moment boxes; its easy counterexample and proportional-tie cases show no node benefit and can favor exhaustive search in wall time. These deliberately synthetic results isolate the oracle's effect and do not establish a general runtime advantage or biological relevance. [UTILITY.md](UTILITY.md) gives the practical use contract and the evidence limits.
# v0.2 extension note

The Welch mode uses routine moment algebra and bounded exhaustive enumeration.
Matrix adapters, delete-one diagnostics and external ranking-table comparison
are usability extensions, not additional mathematical novelty claims. The
conditional two-family theorem and the priority limits below concern paired
signed self-normalized scores only.

## v0.3 extension boundary

The implementation adds shared per-feature membership/sign searches, certified
lower/upper minimum-change intervals, and conservative exact-rational Welch
conditional enclosures. The minimum-variance window fact, endpoint membership
proofs, case-deletion sensitivity, and generic branch-and-bound have prior art.
Automatic R refits, strict donor-level pseudobulk, H5AD/10X ingestion, and
provenance checks are practical engineering extensions, not new DE estimators
or first single-cell robustness claims. The exact paired conditional two-family
characterization remains the specific proposed mathematical contribution.
No publication-level or global-priority novelty claim is implied by a release.
See [bound derivation](ROBUSTNESS.md) and the preserved [prior-art screen](PRIOR_ART.md).
