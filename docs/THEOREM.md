# Exact extrema of a signed self-normalized sum

This document proves a finite optimization result for a descriptive ranking
score. It does not establish a sampling distribution, differential-expression
validity, or the novelty of the result in the mathematical literature.

## Problem and conventions

Let `I` be a fixed multiset of real values that must be retained. Let
`u_1 <= ... <= u_N` be the undecided values in sorted order, breaking equal-value
ties by their original indices. Select exactly `m` undecided positions, where
`0 <= m <= N`, and let `n = |I| + m >= 2`. For a selection `A`, define

\[
 S(A)=\sum_{v\in I}v+\sum_{j\in A}u_j,\qquad
 Q(A)=\sum_{v\in I}v^2+\sum_{j\in A}u_j^2.
\]

The signed self-normalized score is

\[
 R(A)=
 \begin{cases}
 S(A)/\sqrt{Q(A)},&Q(A)>0,\\
 0,&Q(A)=0.
 \end{cases}
\]

Because `Q` is a sum of squares of the actual retained values, `Q=0` implies
that every retained value and `S` are zero. Assigning score zero to this case is
an explicit descriptive extension. A classical t-statistic is undefined for
that all-zero vector.

For `Q>0` and nonzero sample variance, the ordinary one-sample studentized
score of these retained values is

\[
 T(A)=\frac{S(A)\sqrt{n-1}}{\sqrt{nQ(A)-S(A)^2}}
     =\frac{R(A)\sqrt{n-1}}{\sqrt{n-R(A)^2}}.
\]

For fixed `n`, this is strictly increasing in `R` in the interior
`(-sqrt(n), sqrt(n))`. Nonzero constant vectors attain the endpoints
`R=+sqrt(n)` or `R=-sqrt(n)`; their extended studentized scores are respectively
positive or negative infinity. Thus extrema and rankings of signed `R` give
extrema and rankings of signed studentized scores at a fixed retained count,
with the stated descriptive conventions. Absolute-score ranking is a
different target and is not assumed here.

## The two candidate families

The contiguous family consists of every length-`m` block in sorted order:

\[
 C_j=\{j+1,\ldots,j+m\},\quad j=0,\ldots,N-m.
\]

The extreme family consists of a prefix and a suffix of combined size `m`:

\[
 E_j=\{1,\ldots,j\}\cup\{N-(m-j)+1,\ldots,N\},
 \quad j=0,\ldots,m.
\]

The two components of an extreme selection are disjoint since `m <= N`.
Empty components are allowed. When `m=0` or `m=N`, repeated representations
refer to the same selection and can be deduplicated.

**Theorem.** For arbitrary fixed retained values `I`, both the minimum and the
maximum of `R(A)` over all size-`m` undecided selections are attained in
`{C_j} union {E_j}`. There are at most `N+2` candidate representations.

## Lemma 1: a positive maximum has a contiguous witness

Suppose some feasible selection has positive sum. Choose a selection with
maximum positive `R`, and denote its moments by `(S_*,Q_*)`. Then `S_*>0` and
`Q_*>0`. On positive-sum selections, maximizing `R` is equivalent to maximizing

\[
 F(S,Q)=S^2/Q.
\]

This function is convex for `Q>0`. Its supporting linear function at the chosen
optimum, with the constant term omitted, is

\[
 L(S,Q)=2\frac{S_*}{Q_*}S-\frac{S_*^2}{Q_*^2}Q.
\]

Every other positive-sum selection satisfies `F(S,Q) <= F(S_*,Q_*)`.
The convex supporting inequality

\[
 F(S,Q)\ge F(S_*,Q_*)+
 \nabla F(S_*,Q_*)\cdot((S,Q)-(S_*,Q_*))
\]

therefore implies `L(S,Q) <= L(S_*,Q_*)`. Every nonpositive-sum selection with
`Q>0` has `L(S,Q) <= 0`, whereas `L(S_*,Q_*)=S_*^2/Q_*>0`. The all-zero
selection also has `L=0`. Consequently the chosen positive optimum maximizes
`L` over **all** feasible selections, not merely the positive ones.

Fixed values `I` contribute a constant to `L`. Selecting undecided values to
maximize `L` is consequently equivalent to selecting the `m` smallest costs

\[
 (u_j-c)^2,\qquad c=Q_*/S_*.
\]

All values with cost strictly smaller than the cutoff are mandatory and form
a contiguous middle group in sorted order. At the cutoff, equal costs involve
at most two distinct values, one on each side of `c`; if the cutoff is zero,
there is only one value. Equal copies within a value group are interchangeable.

It remains to ensure that a contiguous witness attains the original nonlinear
maximum, rather than only maximizing `L`. After including the mandatory group,
vary the number selected from the left cutoff group. The feasible left counts
form an integer interval. Their `(S,Q)` moments are affine functions of that
count. Every such selection lies on the maximizing support face, so it has
`L=L(S_*,Q_*)>0`; hence its sum is positive and its squared sum is nonzero.
Convexity of `F` along this affine segment implies that its maximum occurs at
an endpoint left count. The original optimum is on this segment, so an
endpoint attains at least its objective; global optimality makes the objective
equal. At an endpoint, one cutoff group is exhausted or the other is unused.
Choosing the appropriate end of each equal-value group makes the selected
positions a contiguous length-`m` block. This proves the lemma, including
cutoff and duplicate-value ties.

Applying Lemma 1 to the negated retained and undecided values proves that a
negative minimum of the original score has a contiguous witness as well.
Contiguity is preserved when sorted order is reversed.

## Lemma 2: an all-positive minimum has an extreme witness

Suppose every feasible selection has strictly positive sum. Start with a
selection minimizing `R`. If its undecided positions already form a prefix
plus suffix, the conclusion follows.

Otherwise, let `l<h` be the smallest and largest unselected sorted positions.
Some selected position `j` lies strictly between them. Hold `I` and the other
`m-1` selected values fixed; their sum and squared sum are `a` and `b`. The score
as a function of the replacement value `x` is

\[
 f(x)=(a+x)/\sqrt{b+x^2}.
\]

Replacing the interior selected value by either `u_l` or `u_h` is feasible.
Both replacement sums are positive by hypothesis, so the numerator is
positive throughout the whole interval `[u_l,u_h]`. The denominator is nonzero
there. Where differentiable,

\[
 f'(x)=\frac{b-a x}{(b+x^2)^{3/2}}.
\]

If `a>0`, the derivative is nonnegative before `b/a` and nonpositive after it:
the only possible interior turning point is a maximum. If `a=0`, the function
is nondecreasing. If `a<0`, positivity of `a+x` implies `x>-a>0`, so
`b-a x>0` and the function is increasing. The exceptional case `a=b=0` has
`f(x)=1` on this positive-numerator interval. Thus there is no interior
minimum, and at least one endpoint replacement has score no larger than the
current score. Since the current selection is globally minimizing, the
replacement preserves the global minimum.

**Termination is finite.** Replacing the selected interior index `j` by one of
the unselected endpoint indices removes that endpoint from the unselected set
and adds `j` inside its previous span. Therefore the integer span
`largest_unselected_index - smallest_unselected_index` strictly decreases.
If the endpoint and interior values happen to be equal, the moments and score
are unchanged, but the index span still strictly decreases. There can be only
finitely many exchanges. When no selected position remains inside the span,
the unselected positions form one contiguous interval, and its complement is
exactly a prefix-plus-suffix selection. With zero or one unselected position,
this form already holds. This proves the lemma for arbitrary fixed `I` and
duplicate values.

Negating all values proves that when all feasible sums are strictly negative,
a maximum has a prefix-plus-suffix witness.

## Completion of the theorem and zero sums

If feasible sums take both strictly positive and strictly negative values, the
maximum is positive and the minimum negative. Lemma 1 and its sign reversal
cover both with contiguous witnesses.

If all sums are strictly positive, Lemma 1 covers the maximum and Lemma 2 the
minimum. If all are strictly negative, their sign reversals cover both.

If all sums are nonnegative and zero is attainable, the minimum score is zero.
The `m` smallest undecided values minimize the sum, so they attain that zero;
they belong to both candidate families. Any positive maximum is covered by
Lemma 1. If every sum is zero, all scores are zero. The analogous argument
with the `m` largest values covers the all-nonpositive case. These cases also
include `Q=0`. If `m=0` or `m=N`, there is a single feasible selection, making
the conclusion immediate. This exhausts the possibilities.

## Computational and certification consequences

Sorting takes `O(N log N)` per feature. Prefix sums of `u` and `u^2` evaluate
each candidate's coupled moments in constant arithmetic work. Scanning at most
`N+2` candidates then computes exact mathematical extrema for a fixed `m` in
`O(N)` arithmetic operations. The bit cost of exact arithmetic is additional.
For several admissible retained counts, apply the theorem separately to each
count and combine the appropriate bounds. At a root deletion budget `b`, this
gives `O(N log N + bN)` arithmetic operations per feature.

Exact **marginal** feature extrema do not make independently attained extrema
simultaneously feasible across genes. A whole-top-K certificate may safely
prune a shared-donor search node if each selected feature's lower score bound
outranks every outsider's upper bound for every admissible retained count.
Overlap requires further shared-donor search. A refutation requires one common
donor deletion set that actually changes the fixed top-K set.

For supplied binary64 values, exact moments mean the exact dyadic values of the
entries, not their floating-point reductions. Comparing two nonzero scores of
the same sign needs no square root: compare `S_a^2 Q_b` and `S_b^2 Q_a`, reversing
the result for negative sums. Signs and zero scores are handled first.

## Evidence and prior-art boundary

The convexity ingredient is standard; see Boyd and Vandenberghe,
[Convex Optimization, quadratic-over-linear functions](https://www.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf).
Contiguous-subset reduction is also classical for univariate minimum covariance
determinant variance minimization; see the original MCD authors'
[review](https://wis.kuleuven.be/statdatascience/robust/papers/2010/wire-mcd.pdf).
Those facts are acknowledged ingredients. This document establishes the
specific two-family signed-ratio result with forced inclusion; it does not
assert that no previous publication contains an equivalent theorem.

The independent checker in `benchmarks/independent_theorem_check.py` directly
enumerates all feasible subsets and the two candidate families without importing
the implementation. It lifts supplied binary64 entries to exact integers with
a common power-of-two denominator, so the verification comparisons are exact.
Its recorded JSON provides the seed, case count, coverage, and failure count.
Finite verification supports the implementation and proof review; it does not
replace the proof or establish research priority.
