# Prior-art screen for ExtremaRank

Screen date: 2026-10-03. This is a targeted screen, not a proof that no earlier result exists. The paper and software sources below are primary sources. Verification status is recorded separately from the inference drawn from each source.

## The proposition being screened

For one feature, let `I` be donors already fixed as retained and `U` the undecided donors. Set

\[
a=\sum_{i\in I}x_i,\qquad b=\sum_{i\in I}x_i^2,
\qquad f(A)=\frac{a+\sum_{i\in A}x_i}{\sqrt{b+\sum_{i\in A}x_i^2}},
\quad A\subseteq U,\ |A|=m.
\]

The proposed exact oracle obtains both signed extrema of `f` from the union of two families after sorting `U` by value: every contiguous block of length `m`, and every subset consisting of `j` smallest plus `m-j` largest values. It is conditional on the fixed retained donors, rather than being only a root-level trimming procedure.

For a retained set with common size `n`, ordinary paired t is

\[
t=\frac{\sqrt{n-1}\,S}{\sqrt{nQ-S^2}},\qquad S=\sum x_i,\quad Q=\sum x_i^2.
\]

When sample variance is positive, this is a strictly increasing transformation of `S/sqrt(Q)`. Thus the oracle addresses a score that gives the same signed within-subset feature ranking as refitted paired t. Changing the retained count changes that transformation; a displayed range of actual t values must account for the count. Undefined all-zero subsets and nonzero constant subsets require an explicit target convention before invoking the theorem.

The downstream question is whether the *entire initial topK set* survives every shared deletion of at most `budget` whole paired donors. Independent feature extrema can give a valid sufficient certificate but may be attained by different donor subsets. Their overlap does not give a refutation. Shared-subset search and independently recomputed witnesses resolve the remaining cases.

## Closest biological-statistics methods

### Adaptive trimmed t statistics: verified author software documentation; paper full text unavailable

Gleiss, Sanchez-Cabo, Perco, Tong and Heinze, *Adaptive trimmed t-statistics for identifying predominantly high expression in a microarray experiment*, Statistics in Medicine 30(1), 52–61 (2011), [DOI 10.1002/sim.4093](https://doi.org/10.1002/sim.4093), [PubMed record](https://pubmed.ncbi.nlm.nih.gov/20963766/).

**FULL TEXT NOT VERIFIED.** Wiley DOI, PDF and ePDF routes were attempted and did not return the paper. A publication list is not a substitute for the methods. However, the authors' institution provides substantive [t.opt software documentation](https://data-science.meduniwien.ac.at/en/services/klinische-biometrie/software/statistische-software/topt/), which was read directly.

That documentation specifies a two-group t statistic using trimmed empirical means and trimmed empirical variances. It chooses a trimming proportion from a predetermined grid by the smallest raw t-distribution p value, repeats the optimization within each permutation, and can stabilize the pooled standard error with a SAM-like constant. The paper abstract states the biological goal as detecting predominantly higher tumor expression.

This overlaps with optimizing a statistic after removal and scanning ordered trims. Neither the verified documentation nor the abstract supplies an all-cardinality-subset exactness theorem, a conditional fixed-donor envelope, or the two-family characterization. This is a limited finding about the material read; it is **not** a claim that the inaccessible paper cannot contain related mathematics. The primary documentation is sufficient to reject any claim that ExtremaRank first recomputes t statistics after trimming.

### MOST: full formula verified

Heng Lian, *MOST: detecting cancer differential gene expression*, Biostatistics 9(3), 411–418 (2008), [DOI 10.1093/biostatistics/kxm042](https://doi.org/10.1093/biostatistics/kxm042), [author preprint PDF](https://arxiv.org/pdf/0709.1307).

The full preprint was read, particularly equations (6)–(9), PDF pages 3–4. MOST orders disease-group expression values from high to low, sums the largest `k` deviations from the control median, divides by a robust scale derived from both groups, then centers and scales using moments of standard-normal order statistics. It maximizes this adjusted statistic over `k`.

Thus it is not ordinary paired t recomputed from first and second moments of each retained donor subset. Its scale is not the retained-subset sample variance, its family is an ordered disease-group tail scan, and its purpose is heterogeneous activation detection. It does not establish the screened conditional two-family oracle. Nevertheless, ordered subset statistics and their use for cancer marker discovery clearly predate this project.

### Continuous fragility: abstract verified

Caldwell, Youssefzadeh and Limpisvasti, *A method for calculating the fragility index of continuous outcomes*, Journal of Clinical Epidemiology 136, 20–25 (2021), [DOI 10.1016/j.jclinepi.2021.02.023](https://doi.org/10.1016/j.jclinepi.2021.02.023), [PubMed abstract](https://pubmed.ncbi.nlm.nih.gov/33684509/).

The method assesses fragility through iterative substitution/reassignment and can use simulations from summary statistics. It concerns a significance threshold rather than exact shared-donor deletion and feature ranking. The existence of continuous-outcome fragility methods prevents framing this project as the origin of continuous-test sensitivity analysis. The paper's full methods were not verified in this screen.

## Exact optimization and robustness auditing

| Source and material verified | Relevant prior result | Scope difference from the screened oracle |
| --- | --- | --- |
| [Shen and Yu, Fractional Programming Part I, 2018, full author PDF](https://www.comm.utoronto.ca/~weiyu/FP_part1.pdf), Theorem 1, equation (7), PDF page 3 | General quadratic ratio transform; maximizing `A/B` can use an auxiliary quadratic expression. | The identity behind the positive-maximization proof is established optimization machinery. The reviewed theorem does not give the signed conditional block/endpoint subset families for this donor statistic. |
| [Freund and Hopkins, Towards Practical Robustness Auditing for Linear Regression, 2023, full HTML](https://arxiv.org/html/2307.16315v1), section 2 | Exact case-removal optimization for OLS coefficient sign, including spatial branch-and-bound with McCormick bounds. | Exact deletion B&B is established. The reviewed target is a coefficient sign, without the screened one-gene studentized oracle. |
| [Rubinstein and Hopkins, ICLR 2025, full proceedings PDF](https://proceedings.iclr.cc/paper_files/paper/2025/file/a31446105fb4789ca1513f55ea14aa44-Paper-Conference.pdf), Problem 1 and Theorem 1.1 | ACRE and OHARE produce valid upper/lower bounds for refitted OLS coefficient changes; favorable-data assumptions govern tightness. | Their formal target is coefficient change, rather than a recomputed paired standard error or signed t-rank envelope. |
| [Konrad and Kuschnig, Finding Most Influential Sets, ICML 2026, proceedings abstract](https://proceedings.mlr.press/v306/konrad26a.html) | Exact global optimization is available for specified linear-fractional influential-subset objectives. | The screened ratio has a square-root second moment. Proceedings abstract verified; full paper not retrieved, so this comparison is limited. |

No claim should be made that a quadratic transform, ratio optimization, branch-and-bound, or exact case-removal optimization is new. The candidate contribution is the complete signed subset-family theorem, its conditional implementation, and its measurable improvement over a general moment-bound search.

## Ranking and bioinformatics removal audits

| Primary source | Verified overlap and claim boundary |
| --- | --- |
| [Shiffman, Giordano and Broderick, Could dropping a few cells change the takeaways from differential expression?, full preprint HTML v1](https://arxiv.org/html/2312.06159v1) | Influence-based deletion sensitivity for differential expression, including subject/sample observations and downstream gene-set rankings, with actual refit checks. This prevents a first donor/cell deletion robustness claim. No screened exact paired-t oracle was identified in the material read. |
| [Huang, Shen, Wei and Broderick, Dropping Just a Handful of Preferences Can Change Top Large Language Model Rankings, full preprint HTML v2](https://arxiv.org/html/2508.11847v2) | TopK membership stability uses comparisons across the original topK boundary; influence/one-step approximations propose deletions and refitting checks them. The universal boundary-pair reduction is prior art. |
| [Alzhrani, Decision Reliability Profiles for Auditing Targeted Case-Removal Sensitivity in AI Model Comparisons Under Additive Metrics, 2026](https://doi.org/10.3390/a19080676) | Exact additive pairwise reversal budgets sort winner-supporting case contributions. Paired statistical tests are separately reported. Additive winner margins differ from the screened variance-refitted statistic. Publisher-indexed methods were read; full downloadable PDF not verified. |

## Conditional self-normalized subset theorem screen

### Univariate MCD: contiguous variance-minimizing subsets are established

Hubert and Debruyne, *Minimum covariance determinant*, WIREs Computational Statistics 2, 36–43 (2010), [full author-hosted review](https://wis.kuleuven.be/statdatascience/robust/papers/2010/wire-mcd.pdf), DOI [10.1002/wics.61](https://doi.org/10.1002/wics.61). The full PDF was read, including the univariate subsection on PDF page 5.

For univariate data, the review describes finding the size-`h` subset with smallest variance by scanning contiguous sorted subsets, computing means and variances recursively, in `O(n log n)` total time. It connects this objective to least trimmed squares. That prior result directly precludes a first contiguous-subset reduction or first sorted-window algorithm claim.

The target minimizes `Q-S^2/h`, rather than extremizing the signed ratio `S/sqrt(Q)`. The reviewed univariate result does not state the complementary prefix-plus-suffix family, all signed extrema, or the forced-inclusion first/second-moment offsets used here. This distinction supports a specific technical comparison, not proof of priority for the present theorem.

Search families included ordinary/paired/studentized t subset extrema; self-normalized sum subset selection and deletion; minimum/maximum coefficient of variation subsets; contiguous blocks and prefix/suffix extrema; adaptive trimmed t; MOST; and exact fractional subset optimization. These searches did not locate a primary source proving the same union-of-two-families characterization with fixed included first and second moments.

Several apparent matches concern a different object: selecting covariates rather than donors, normal-population ranking, choosing survey strata, probability/concentration bounds for self-normalized sums, or scanning intervals already imposed by the problem. They do not establish exact extrema over arbitrary fixed-cardinality subsets. They were not used as evidence for novelty.

The positive-maximization argument is closely related to the established quadratic transform. The negative/positive sign split, opposite endpoint family, and extension with fixed included moments must therefore be the theorem's explicit content. A proof of only the positive root maximum would provide a substantially weaker contribution.

## Research verdict and remaining checks

The specific candidate is not occupied by any result identified in this screen. That statement is narrower than mathematical priority. There remains a real verification gap for the Gleiss paper and for the full Konrad–Kuschnig paper. A manuscript submission should obtain and check those methods and references before making a priority assertion.

The software can responsibly explain its concrete technical contribution without a global first claim: it replaces a combinatorial one-feature envelope with a proved finite candidate family and uses that exact conditional envelope inside a shared-donor ranking audit. Scientific value still requires actual benchmark evidence showing that the envelope is correct, tighter than moment boxes, and useful in a nontrivial shared-subset search. A paper citation screen alone does not establish those results.
