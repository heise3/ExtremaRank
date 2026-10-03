# ExtremaRank

[Repository](https://github.com/heise3/ExtremaRank) ·
[Releases](https://github.com/heise3/ExtremaRank/releases) ·
[Report an issue](https://github.com/heise3/ExtremaRank/issues)

**Which of your paired omics candidates survive removing a few donors?**

ExtremaRank computes exact worst-case studentized scores and audits whether a
Top-K feature set can change under shared donor deletion. It recomputes both
mean and variance. A counterexample includes the actual donor IDs and the
replacement feature list.

The core has **no third-party dependencies**. Python 3.10+ is sufficient.
CSV/TSV inputs make it usable with paired transcriptomics, proteomics or
metabolomics whenever an ordinary paired-t ranking is the desired target.

[中文结果与使用说明](docs/REPORT_CN.md) · [Usage](docs/USAGE.md) ·
[Proof](docs/THEOREM.md) · [Prior art](docs/PRIOR_ART.md) ·
[Utility and limits](docs/UTILITY.md) · [Interactive report](results/report.html)

## Install and run

```bash
python -m pip install .
extremarank examples/variance_reversal.tsv --budget 1 --top-k 1 --output my-audit
```

For your data, put donor IDs in the first column and one feature per remaining
column. Values are **already prepared within-donor target-minus-reference
effects**. Freeze normalization, pairing and the feature universe before
auditing. The software checks unique IDs and finite rectangular data.

```text
donor_id    gene_A    gene_B
donor_1     2.0       1.0
donor_2     2.1       1.1
donor_3     1.9       0.9
donor_4    10.0       2.0
```

```bash
# Both up- and down-regulated candidates, ranked by absolute paired t.
extremarank my-effects.csv --direction absolute --budget 2 --top-k 20 --output my-audit

# Exact gene envelopes without running the potentially larger shared search.
extremarank my-effects.csv --budget 4 --skip-topk --output gene-envelopes
```

The output directory contains:

- `summary.tsv`: baseline rank, worst/best signed t, robust positive/negative
  direction, and donor witnesses with recomputed mean and variance.
- `envelopes.tsv`: exact signed score extrema for each deletion count.
- `audit.json`: `CERTIFIED`, `REFUTED` with one shared deletion witness, or
  `UNRESOLVED` at the declared node limit; input hash and parameters included.

`CERTIFIED` concerns the unordered Top-K set over **every deletion of at most
the budget**. It does not say the internal order is unchanged. All-zero values
use descriptive score zero; nonzero constant vectors use extended t of ±∞.
These conventions are recorded explicitly and do not define inferential p-values.

## The algorithmic contribution

For retained effects, let `S = sum(x)` and `Q = sum(x²)`. At a fixed retained
count, signed paired-t ranking is equivalent to `R = S/sqrt(Q)`.

Our conditional two-family extremum theorem shows that, even when some donors
are forced to stay, both extrema over all fixed-size subsets occur among:

1. Contiguous blocks of sorted optional donor values.
2. A prefix plus a suffix of the same combined size.

The union has at most `D+2` representations, compared with `choose(D,m)` subsets.
Prefix moments evaluate it in linear arithmetic work after sorting. The
implementation removes duplicate end blocks. Exact dyadic integer sums and
cross-products make every selection/ranking decision exact for the supplied
binary64 input; decimals shown in reports are only displays.

This oracle also bounds each node of a shared-donor search. Independent feature
interval overlap is never treated as a ranking counterexample. The full search
can still be exponential; a finite node limit reports `UNRESOLVED` honestly.
Absolute ranking uses exact absolute maxima but a conservative zero lower bound
when the signed interval crosses zero.

The contribution is the specific conditional signed two-family characterization,
its exact oracle, and its use in this audit. Ordered trimming, robustness audits,
quadratic transforms, and generic branch-and-bound all have prior art. The
[comparison](docs/NOVELTY.md) states those overlaps and retrieval gaps. This is
an independently derived, proved and tested research implementation, without a
claim of global mathematical priority or peer-reviewed novelty.

## Executed public-data benchmark

![Exact candidate counts after deleting four donors](results/candidate_scaling.svg)

Both studies use verified biological pairing and public NCBI count matrices.
The prepared paired log2(CPM+1) effects and the complete provenance are bundled.
Deletion budgets `0..4`, the filter and the matched first 100 eligible features
were declared before evaluating outcomes.

| Dataset | Genuine pairs | Eligible genes | Whole-matrix scalar extrema, 5 counts | Matched speed vs exact enumeration |
|---|---:|---:|---:|---:|
| GSE87290, PBMC baseline → LPS | 14 | 13,328 | 1.385 s | 22.24× |
| GSE50760, normal colon → primary CRC | 18 | 16,597 | 2.077 s | 55.29× |

These are one-machine timings on arm64 macOS/Python 3.12. Matched timings use
the same exact arithmetic, include sorting, and compare 100 fixed features at
each of five counts: **0 extrema mismatches in 1,000 calls**. Whole-matrix times
exclude preparation, serialization, ablations and joint Top-K search. They are
not speed comparisons against approximate DE tools.

The full public CLI, including input parsing, scalar envelopes, shared Top-20
audit and output writing, took 2.85 s / 3.94 s on the two matrices. Executed TSV
tables in this distribution are losslessly compressed as `.tsv.gz`; the
archive manifest records both byte hashes. Rerunning the CLI/benchmarks writes
ordinary TSVs. The offline HTML embeds its displayed results and needs neither.

At four deletions, block-only candidates miss true lower extrema for 15.05% /
22.27% of features; endpoint-only candidates miss upper extrema for 4.88% /
6.46%. The two-family result therefore changes actual answers in both datasets.

The full eligible Top-20 set has verified one-donor counterexamples in **all
three ranking directions in both studies**. For upward ranking, deleting
`PBMC_replicate_5` changes one list member, and deleting `AMC_6` changes four.
Each baseline and replacement ranking was replayed independently with Fraction
arithmetic. This demonstrates a usable sensitivity report, not biological
ground truth or a reason to remove those donors from the experiment.

Independent checks include 200,000 exact scalar scenarios with forced inclusion
and exhaustive shared-subset ranking tests. See the JSON files in `results/`.

## Reproduce offline

```bash
python -m unittest discover -s tests -v
python benchmarks/independent_theorem_check.py
python benchmarks/search_ablation.py
python benchmarks/scaling.py

# Only raw-count preparation / full real-data reproduction needs NumPy.
python -m pip install '.[benchmark]'
python data/fetch_sources.py --offline
python benchmarks/real_data.py
python benchmarks/real_topk.py
```

No GPU, FASTQ processing, private database login, or network access is needed
for these bundled reproductions. Missing source files can be restored with
`data/fetch_sources.py`, which checks pinned SHA256 values and refuses drift.

## Scope

Use the tool for fixed paired effects and ordinary paired-t feature ranking,
especially small paired studies where a candidate shortlist should be audited
before follow-up work. It provides deterministic worst-case input sensitivity.
It does not replace limma/DESeq2, estimate FDR, model cell types, certify
moderated t or covariate-adjusted models, or re-estimate upstream normalization.
It does not measure prediction accuracy or guarantee biomarker validation.

Software: MIT. Third-party data attribution: [data/DATA_LICENSE.md](data/DATA_LICENSE.md).
Questions, counterexamples to the theorem, independently verified datasets and
new target models are useful contributions; see [CONTRIBUTING.md](CONTRIBUTING.md).
