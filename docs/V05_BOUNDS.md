# Fixed-value conditional variance floor

For a retained set of n values, the integer moment numerator of the sample
variance is `n * sum(x^2) - sum(x)^2 = sum_{i<j}(x_i-x_j)^2`.
Partition retained values into a fixed set F of size f and m selected optional
values O. The pairwise numerator decomposes into:

1. fixed-fixed: `f*Q_F - S_F^2`;
2. fixed-optional: a sum of `f*x^2 - 2*S_F*x + Q_F` over selected optional x;
3. optional-optional: `m*Q_O - S_O^2`.

The first is known exactly. A lower bound for the second is the sum of the m
smallest fixed-optional costs over all optional values. A lower bound for the
third is the minimum pairwise numerator over sorted contiguous m-value windows
(zero for m <= 1). The sum of these three bounds is a lower bound for every
admissible retained set. Separate minima need not be jointly attainable.

The implementation takes the maximum of this floor and the pre-existing
unconstrained-window and moment floors. It therefore never weakens the old
variance lower enclosure. Conversion to variance of the mean divides by
`n^2*(n-1)` with exact rational arithmetic. Independent exhaustive validation
checks every feasible fixed/optional configuration in the small test suite,
including tiny and extreme binary64 values. This is a valid enclosure, not a
claim of exact joint mean/variance extrema or global novelty priority.
