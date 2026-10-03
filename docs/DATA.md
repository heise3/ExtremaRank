# Real paired RNA input and computational benchmark

The repository includes two small public processed RNA datasets and checksummed
prepared paired effects. They serve as real inputs for an exact deletion
robustness audit. They do not provide biological ground truth for gene discovery.

## Sources and biological units

| Dataset | Source count matrix | Complete paired units used | Eligible genes | Fixed contrast |
| --- | --- | --- | ---: | --- |
| [GSE87290](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87290) | 39,376 genes × 29 samples | 14 independent human donors | 13,328 | Intravenous LPS minus baseline |
| [GSE50760](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760) | 39,376 genes × 54 samples | 18 independent patients | 16,597 | Primary colorectal cancer minus normal colon |

These are **NCBI-generated raw gene counts**, aligned to GRCh38.p13. They are
distinct from the original submitter's Cuffdiff summaries or FPKM files. See
[NCBI's count pipeline documentation](https://www.ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html).
The source bytes and official GEO SOFT metadata were retrieved on 2026-10-03 UTC;
`data/sources.json` records their URLs, exact sizes and SHA256 hashes.

GSE87290's official record describes 15 healthy subjects sampled before and
after a 1 ng/kg LPS challenge. Repeated `PBMC_replicate_n` identifiers in original
sample titles define the pairs. The source paper describes the in vivo
inflammatory comparison ([Lin et al., 2016](https://pubmed.ncbi.nlm.nih.gov/27230130/));
the original GEO tissue annotation is PBMC. The NCBI count matrix lacks
`GSM2326903`, the LPS sample from donor 2. That donor's available baseline
`GSM2326888` is also excluded before preparing any effects. The benchmark thus
uses **14 complete pairs**, with matching within-pair sex, race and responder
group metadata. These are selected high/low responders in the original study;
this is not a representative population sample.

GSE50760 contains three tissues from each of 18 patients with colorectal cancer:
normal-looking surrounding colonic epithelium, primary cancer and liver
metastasis ([Kim et al., 2014](https://pubmed.ncbi.nlm.nih.gov/25049118/)). Donor
identifiers are `AMC_2`, `AMC_3`, `AMC_5`, `AMC_6`, `AMC_7`, `AMC_8`, `AMC_9`,
`AMC_10`, `AMC_12`, `AMC_13`, `AMC_17`, `AMC_18`, `AMC_19`, `AMC_20`, `AMC_21`,
`AMC_22`, `AMC_23`, and `AMC_24`. Titles ending `-1`, `-2`, and `-3` indicate
primary cancer, normal colon, and metastasis, respectively. Only the normal and
primary pair is used. All 36 required sample columns occur in the count matrix.
The original GEO titles and source fields are authoritative for this pairing;
the EBI mirror labels the `AMC_2-3` metastasis incorrectly as normal.

The original source files retain their source status. The software license does
not relicense third-party data. Cite the dataset generators and GEO/NCBI when
using these inputs; no per-file MIT or CC license is inferred from public
download access.

| Count file | Bytes | SHA256 |
| --- | ---: | --- |
| GSE87290_raw_counts_GRCh38.p13_NCBI.tsv.gz | 1,333,166 | `030e3a59ae8d3dcc8d3861d468bb11cb2acb19af70c635dfacdb760581c14bd2` |
| GSE50760_raw_counts_GRCh38.p13_NCBI.tsv.gz | 2,358,563 | `ffa68762e84954bd038fda10c408b0eecc0ba725f181dd15c8dbcbbfc44ebeb1` |

Both files have the same ordered universe of 39,376 NCBI GeneIDs. The SHA256 of
one GeneID per line, with a terminal newline, is
`454511bd3fda79fc846d197e654ff94946ab45d2bdce6bf6ea07d550558e2299`.

## Frozen preprocessing contract

Preprocessing is applied once, before any deletion or extrema evaluation:

1. Select complete biological donor pairs for the specified contrast. Exclude
   incomplete donors and other tissues before preparing effects.
2. For each selected sample, compute its library denominator from **all 39,376
   source genes**. Define `CPM = count × 1,000,000 / library_denominator`.
3. Keep a gene if `CPM >= 1` in at least half of the selected samples: at least
   14 of 28 samples for GSE87290, or 18 of 36 for GSE50760.
4. Compute `log2(CPM + 1)` with a fixed pseudocount of one. Form target minus
   reference within each biological donor.
5. Freeze this gene universe, sample pairing, transform and donor-by-gene
   effect matrix. Preserve source GeneID row order.

Deletion removes entire paired donor rows. The signed paired t score recomputes
both the retained mean and retained sample variance; the normalized input
effects remain fixed. This simple transform demonstrates computational
robustness. It does not refit a negative-binomial count model, estimate
normalization uncertainty or borrow variance across genes.

Prepared files appear under `data/prepared/<accession>/`:

- `effects.csv.gz`: donor rows, GeneID columns, binary64 values written with
  17 significant digits for round-trip recovery; deterministic gzip timestamp.
- `gene_ids.txt`: eligible GeneIDs in original source order.
- `donor_pairs.tsv`: exact source sample IDs and donor pairing.
- `preparation.json`: filter, pairing, exclusions, library denominators,
  transformation, source hashes and prepared-input hashes.

## Reproducing offline

From the repository root, with Python 3.10+ and NumPy available:

```sh
python data/fetch_sources.py --offline
PYTHONPATH=src python benchmarks/real_data.py
```

The default run makes no network requests. It verifies bundled source hashes,
reconstructs the frozen effects and runs the real benchmark. Missing source
files can be restored using `python data/fetch_sources.py`; downloads must match
the pinned hashes. Updated remote metadata that differs from the snapshot is
rejected. A prepared-input-only rerun is available:

```sh
PYTHONPATH=src python benchmarks/real_data.py --prepared-only
```

`--prepare-only` rebuilds the paired effects without running the algorithm.
`--dataset GSE87290` or `--dataset GSE50760` restricts the run.
`--output-dir <directory>` selects a destination for results.

## Fixed comparison and practical audit

The benchmark contract fixes maximum deletion at four donors, using retained
counts `n, n-1, ..., n-4`. The two-family algorithm is run on **every eligible
gene**. The first 100 eligible source rows are selected deterministically for
the exhaustive comparison, before evaluating scores. Both methods use the
same core exact integer score primitives. Timed calls include sorting and
exact arithmetic; preprocessing, serialization and ablations are outside the
matched runtime. Execution order alternates deterministically. Timing is one
matched run and depends on Python and hardware.

Each minimum and maximum is compared with exhaustive enumeration without a
floating tolerance. Equal extremal values may have different valid witnesses.
The benchmark validates exactness on binary64 paired-effect inputs, rather
than asserting exactness of the preceding floating normalization as real
arithmetic.

Both single-family ablations are evaluated against exhaustive enumeration on
the same first 100 genes and against the full two-family extrema on every
eligible gene. The output records wrong lower and upper bounds by gene and
deletion budget. Ablations are diagnostic demonstrations of the need for both
families; neither is a valid certificate method by itself.

For practical inspection, `<accession>_gene_audit.tsv` contains every eligible
gene's baseline signed t, mean and sample variance, its worst signed t under
up to four donor deletions, the mean and variance at that witness, the actual
removed donor IDs, and robust-positive/negative flags. Cross-budget t comparison
uses exact signed rational squared-t comparisons, including retained sample
count. `<accession>_baseline_top20.tsv` gives the 20 largest baseline signed-t
genes and the same audit fields. Ties use lexicographic GeneID order.

A robust-positive flag means the paired mean remains positive under every
allowed whole-donor deletion. It is a descriptive stability property, with no
significance, disease-marker, causal or clinical interpretation. This audit
helps prioritize follow-up by showing exactly which donor removals weaken an
observed ranking or directional effect. All-zero effects receive a descriptive
zero score; constant nonzero effects can display infinite t. Neither convention
is an inferential test.

`results/real/real_benchmark.json` and per-dataset summaries record checksums,
environment, call counts, measured timings, exact-oracle equality and practical
audit summaries. Detailed extrema, matched comparison, ablation failures and
gene audit tables are saved as TSV.

The executed reference run had zero minimum/maximum mismatches in 500 matched
gene-budget calls per dataset. On the recorded host, aggregate matched-call
runtime ratios were 22.24× for GSE87290 and 55.29× for GSE50760; these are measured
single-run values, not portable performance guarantees. Full two-family calls
took 1.38 s and 2.08 s, respectively, excluding preprocessing, ablations and
serialization.

The full-spectrum ablation establishes practical need for the combined family:

| Dataset, exactly four deletions | Contiguous-only wrong minimum | Ends-only wrong maximum |
| --- | ---: | ---: |
| GSE87290 | 2,006 / 13,328 (15.05%) | 651 / 13,328 (4.88%) |
| GSE50760 | 3,696 / 16,597 (22.27%) | 1,072 / 16,597 (6.46%) |

The practical audit found 2,616 robust-positive genes in GSE87290 and 6,538 in
GSE50760. Every baseline top20 gene in these two inputs remains directionally
positive through four donor deletions. Directional robustness does not imply
that its top20 membership is stable; a joint ranking audit is needed for that
separate question.
