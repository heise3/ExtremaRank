# Native R validation, v0.4.0

The standalone package is in [`r/extremarank`](../../r/extremarank). Its R/C++
algorithm has no Python runtime or subprocess bridge. Development validation
uses an independent Python Fraction oracle; installed-package tests run solely
in R and replay frozen rational golden results.

| Verification | Actual result |
|---|---|
| Local standard `R CMD check` (R 4.6.1, Apple clang 17, C++17) | **Status: OK**, zero errors, warnings and notes; [log](R_CMD_check.log) |
| Small exact and extreme-value cases | 1,004 cases, 6,024 design/direction problems, **439,536 independent checks**, zero failures; [record](native_fraction_validation.json) |
| R-only independent enumeration | 5,040 property/cap/witness checks, 100 exhaustive conditional extrema cases, 1,632 frozen Fraction golden property results; sparse inputs, IDs, experiment objects, DFS/indexed proofs and stable 100-sample Welch pruning |
| Real full-feature ranks | **183** baseline and distinct witness ranks matched the independent centered-variance Fraction oracle; [record](real_fraction_validation.json) |
| Native R pseudobulk | **4,347,470** integer sums matched independent R Matrix row sums, across eight cell types; 24,679 provider singlets, six missing cell-type annotations explicitly excluded |
| Established models | Direct R limma, limma-voom, edgeR and DESeq2 synthetic/covariate/pairing/cap tests passed; failed designs were preserved |
| Existing Python regression | 59 tests, no failures; two optional H5AD tests skipped locally in this run and exercised by the dedicated CI job |

Upward candidate results at deletion budget 2:

| Example | Units / features | Individually certified Top-K members | Certified effect signs | Exact minimum whole-Top-K change |
|---|---|---|---|---|
| Golub microarray | 38 samples / 3,051 probes | 6 / 20 | 20 / 20 | 1 |
| Kang B-cell donor effects | 8 donors / 35,635 genes | 6 / 20 | 20 / 20 | 1 |
| CPTAC spike-in measurements | 6 technical runs / 1,097 proteins | 4 / 20 | 20 / 20 | 1 |
| Nutrimouse lipids | 40 mice / 21 features | 4 / 5 | 5 / 5 | **2** |

CPTAC deletions describe technical measurement sensitivity, not biological
replication. Actual R audit properties match the published v0.3 results in all
three ranking directions. Certificates apply to each frozen R binary64 input;
cross-language preprocessing is not claimed to be bit identical. Hexadecimal
input tables permit exact independent replay.

Model outcome details are in [real_validation.json](real_validation.json):

- Golub limma: 38 valid refits; Nutrimouse diet-adjusted limma: 40 valid refits.
- Kang limma-voom: eight valid refits for both the default count filter and the
  previously published count filter.
- Kang edgeR: the baseline ranking is `NOT_EVALUABLE`, with two valid and six
  failed refits, for both filters. Negative raw QLF statistics are retained.
- Kang DESeq2: the default filter produces seven valid and one failed refit.
  An actual ranking change is observed in valid refits, so overall status is
  `OBSERVED_CHANGED`; this does **not** make the failed scenario evaluable.
  The previously published filter (`min_count_samples = 3`, from
  `benchmarks/validate_v03.py`) produces eight valid refits. This is a separately
  labeled replay of an existing contract, not an outcome-based repair of the
  default analysis. Both full outputs are preserved.

Each result directory contains `result.rds`, usable feature or fit-status TSVs,
scores and witnesses. Large TSVs are gzip-compressed without loss. Each dataset
stores one shared `frozen_hex.tsv` or `frozen_hex.tsv.gz` under its upward run.
RDS records retain all nested metadata, exact strings and one-based indices.
Observed model refits never certify unseen deletions. Ranking sensitivity does
not establish FDR, differential-expression validity or clinical usefulness.

Reproduce from this repository with installed R dependencies and optional models:

```sh
R CMD INSTALL r/extremarank
R CMD build r/extremarank
R CMD check --no-manual --install-args=--strip extremarank_0.4.0.tar.gz
Rscript benchmarks/native_r_real.R .
# Independent development oracle; Python is not part of the R runtime.
python benchmarks/validate_native_r.py --cases 1000 --work /tmp/extremarank-r-check
python benchmarks/validate_native_r_real.py --work /tmp/extremarank-real-check
python benchmarks/compact_native_r.py
```

The [native R CI](../../.github/workflows/native-r.yml) checks source compilation
and R-only tests on Linux, macOS and Windows, plus direct Bioconductor models on
Linux. Optional packages can be absent in core-only checks; model tests then run
in the dedicated job. Public source/example attribution remains in
[`DATA_LICENSE.md`](../../DATA_LICENSE.md).
