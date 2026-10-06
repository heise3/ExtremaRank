# Individual certificates and minimum deletion changes (v0.3)

`extremarank robustness` accepts the same matrix/metadata preparation contract as
`study`. It audits initial Top-K candidates by default. `--features all` queries
every feature; `--feature-id ID` may be repeated to query a declared subset.
Every queried candidate is compared against the **entire frozen universe**.

```bash
extremarank robustness data/prepared/Golub/matrix.csv.gz \
  --metadata data/prepared/Golub/metadata.tsv --target AML --reference ALL \
  --design welch --budget 2 --top-k 20 --output feature-audit
```

`study --feature-audit` and the original paired-effects CLI `--feature-audit`
provide the same property reports. A `study` feature selector automatically
enables this audit. Ranking directions are `up`, `down`, and `absolute` with
ascending lexical feature-ID ties. Preparation is frozen; normalization is not
refitted in native audits. Paired deletion removes a complete donor pair.
Independent-group deletion removes one independent sample and must retain at
least two observations per group.

## Interpretation

`feature_robustness.json` and `.tsv` contain two properties per queried feature:

- **Membership:** does its initial inside/outside Top-K status survive every
  feasible shared deletion of at most the requested budget?
- **Direction:** does its initial target-minus-reference mean-effect sign
  (positive, zero or negative) survive those deletions? Reaching zero counts
  as a change from a positive or negative effect.

Each property is `CERTIFIED`, `REFUTED` by an actual common deletion, or
`UNRESOLVED` under the declared search limits. Per-feature certification can
succeed while the whole Top-K set is refuted. Certification of an outside
feature concerns remaining outside, not becoming a biological negative.

`minimum_change_lower_bound = L` excludes changes at every feasible deletion
count below L. An observed witness of size U supplies
`minimum_change_upper_bound = U`. Only L=U receives `exact_minimum_change`.
With no change through budget B, the lower bound is B+1 and the upper bound
is null; this is **not** a claim that a change exists at B+1. An unresolved
search retains its valid interval. `certified_through` is the number of
excluded deletion counts. `topk_minimum_change` concerns any set change and
needs all original members queried to derive its full lower bound.

The search proceeds by increasing deletion count and shares evaluated subsets
and conditional nodes across properties. It tries a few actual subsets,
checks sufficient interval proofs, and enumerates small residual cardinalities
(at most 2,000 feasible subsets). Larger cases use branch-and-bound. These are
cost choices; no approximation changes the target. `--max-nodes` and
`--max-subsets` are global per audit, not per feature, and not timeouts.
All-feature audits can be costly. The Python Welch API retains enumeration
as its default for compatibility; `study` defaults to the new safe search.
`--welch-search enumeration` selects the old reference path.

## Why the new bounds are safe

The paired case uses the existing conditional two-family theorem in
[THEOREM.md](THEOREM.md). The Welch case uses enclosing rational intervals,
not a claimed two-group extension of that extremum theorem.

Within one group, let forced retained values have sums Sf/Qf and select m
optional values from a sorted list. With n retained observations:

1. Smallest/largest m values bound the sum S. Smallest/largest m squared
   values independently bound Q.
2. These rectangles enclose V=nQ−S². Thus
   `max(0,n*Qlo−max(Slo²,Shi²))` is a variance-numerator floor, and
   `n*Qhi−min(S² over [Slo,Shi])` is a ceiling.
3. Relax the forced-inclusion constraint and consider every n-element subset
   of the allowed values. Its minimum V is found among sorted windows of n
   values and gives another safe floor for the constrained problem.

The window fact follows from minimizing centered squared deviations. For a
fixed center, the n closest observations can be chosen as a contiguous sorted
window. A globally minimum subset uses its own mean as center. Replacing it
by a closest window cannot increase its squared error; optimizing that
window's center can only decrease it further. Ties admit a window choice.
This is classical univariate minimum-variance trimming, not a priority claim.

The group mean is S/n and its contribution to Welch standard-error squared is
V/[n²(n−1)]. Combining the group boxes bounds the mean difference Δ and
standard-error squared W. Signed Δ/sqrt(W) increases with Δ, and its
monotonicity in W depends on the sign of Δ. The implementation selects the
appropriate endpoints with exact rational comparisons. Absolute ranking uses
zero as a conservative lower bound if the signed interval crosses zero.
Mean/variance endpoint combinations need not be jointly attainable.

For an initial member j, fewer than K possibly outranking features proves
membership. For an initial outsider, at least K definitely outranking
features proves exclusion. Endpoint ordering includes the lexical tie rule.
A conditional node is discarded only after these sufficient proofs hold.
**Overlapping intervals never refute a property**; only an evaluated shared
deletion does. Fully eliminating every smaller cardinality justifies an exact
minimum. Resource caps leave unexamined nodes unresolved.

Decisions use exact integer dyadic moments for supplied finite binary64 inputs.
The descriptive conventions for zero means and zero variance remain unchanged.
A canonical hash binds shape, ordered IDs, native design, normalized group
labels and row-major little-endian binary64 values. Unit IDs are exported.
Displayed decimal scores are not used for decisions.

## Validation and limits

The definition-based independent checker uses Fraction means and centered
sample variances, complete small deletion universes, forced-retained/deleted
nodes, lexical ties and deliberately capped searches. It tests both native
designs and all three ranking directions. See
[the checker](../benchmarks/independent_robustness_check.py),
[results](../results/independent_robustness_v03.json), and
[real validation](../results/validation_v03_real.json).

Favorable separated Welch examples demonstrate pruning. Difficult overlapping
or tied intervals can still require exponential work. Native certificates
concern frozen ranking and mean-effect signs; they certify neither p values,
FDR, discovery truth, prediction accuracy nor clinical usefulness. These
search and membership principles have prior art; the paired theorem remains
the specific research contribution described in [NOVELTY.md](NOVELTY.md).
