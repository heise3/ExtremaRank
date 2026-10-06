# Native R validation, v0.5.0

The native R/C++ package runs without Python. Independent Python Fraction
validation is a development oracle; installed-package property tests run in R.

| Verification | Actual result |
|---|---|
| Standard local R CMD check, R 4.6.1 / Apple clang 17 / C++17 | **Status: OK**, zero errors, warnings and notes; [log](R_CMD_check.log) |
| Independent complete deletion oracle after the Welch bound change | **439,536 checks**, 1,004 cases, 6,024 paired/Welch direction problems, zero failures; [record](native_fraction_validation.json) |
| Native R-only tests | 5,040 property/cap/witness checks; 100 conditional extrema cases; 1,632 frozen Fraction golden checks; DFS/indexed proofs and the stable 100-sample Welch case passed |
| New workflows | Exact grouped counts matched independent enumeration; coverage/caps, candidate ranks/witnesses, disk/none modes, moved exports and checksum corruption detection, changed checkpoint rejection, contrasts, custom repeated observations and delayed/HDF5 aggregation passed |
| Parallel reproducibility | Serial versus two-worker PSOCK stochastic refits produced identical scores and candidate results; caller RNG state preserved |
| Four model adapters | limma, limma-voom, edgeR and DESeq2 covariate/paired/capped tests passed; fixed treatment coding also matched with global sum contrasts |
| Frozen real native analyses | Four datasets times three ranking directions: **12 analyses**; both search orders matched the v0.4 candidate properties, baseline Top-K and minimum-change bounds |
| Frozen real model analyses | **Five analyses** matched v0.4 baseline, every valid/failed scenario and every observed Top-K; [record](real_validation.json) |
| Existing Python tests | 59 tests, no failures; four optional tests skipped locally; dedicated CI jobs exercise their supported dependencies |

Model coverage and returned object sizes with explicitly selected
`store_scores="none"`:

| Analysis | Valid / failed / not run deletions | Result | v0.4 memory object | v0.5 summary object |
|---|---|---|---|---|
| golub-limma | 38 / 0 / 0 | `OBSERVED_CHANGED` | 22.38 MB | 0.18 MB |
| nutrimouse-limma | 40 / 0 / 0 | `OBSERVED_CHANGED` | 0.30 MB | 0.14 MB |
| kang-b-limma-voom | 8 / 0 / 0 | `OBSERVED_CHANGED` | 15.68 MB | 1.98 MB |
| kang-b-edgeR | 2 / 6 / 0 | `NOT_EVALUABLE` | 10.41 MB | 1.96 MB |
| kang-b-DESeq2 | 7 / 1 / 0 | `OBSERVED_CHANGED` | 10.42 MB | 1.98 MB |

Sizes are `object.size()` of the returned objects on this host, with decimal MB.
The summary mode deliberately omits complete score tables; this is not a peak
RAM comparison or a claim of a uniformly faster algorithm. Candidate summaries,
actual witnesses, coverage and failed feature IDs are retained. Use `memory`
or `disk` when complete scores must remain available.

The existing Kang edgeR baseline remains `NOT_EVALUABLE`; six failed deletions
are retained. Kang DESeq2 retains one failed deletion despite an observed
Top-K change. No filter was adjusted to rescue a favorable result. Overall
status describes Top-K membership; candidate direction status is separate.

The new Welch floor uses exact pairwise squared-difference lower bounds under
fixed inclusion constraints; [derivation](../../docs/V05_BOUNDS.md). Search-order
heuristics only affect traversal. Some examples take fewer nodes with the old
order, so the new order is optional. No global novelty or universal speed claim
is made. The Nutrimouse upward example still requires exactly two deletions
for a whole-Top-K change, illustrating a change missed by delete-one diagnostics.

Native certificates refer to declared frozen binary64 inputs and deletion
budgets. Fitted-model outputs describe executed refits; failed or capped refits
cannot be treated as certificates. CPTAC here contains technical measurement
runs, not biological replicates. Ranking sensitivity does not establish FDR,
clinical usefulness or a biological finding.

Delayed/HDF5 raw cell inputs are validated and summed in row blocks. The smaller
feature-by-biological-unit output is materialized. Memory guards cover selected
numeric storage, not all model or arbitrary-precision work memory. Model
checkpoints resume cached refits, not native DFS traversal.

Reproduce in a clean checkout after installing the documented R dependencies:

```sh
R CMD INSTALL r/extremarank
R CMD build r/extremarank --no-manual
R CMD check --no-manual --install-args=--strip extremarank_0.5.0.tar.gz
EXTREMARANK_TEST_PARALLEL=true Rscript r/extremarank/tests/workflows.R
# Development-only independent oracle:
python benchmarks/validate_native_r.py --cases 1000 --work /tmp/er-v05 --output /tmp/er-v05-validation.json
```

`Rscript benchmarks/validate_v05.R .` replays the frozen v0.4 real results and
writes five new result directories. Move the existing five model directories
under `results/v05` out of the checkout first; exports refuse to overwrite a
nonempty directory. The [CI workflow](../../.github/workflows/native-r.yml)
checks Linux, both macOS architectures, Windows, R 4.1/oldrel/devel and optional
Bioconductor/delayed workflows. Platform binaries are reinstalled in isolated
libraries before artifact upload. Release assets and SHA-256 values are listed
on the [v0.5.0 release](https://github.com/heise3/ExtremaRank/releases/tag/v0.5.0).
See [Chinese usage](../../docs/V05_CN.md) and [data attribution](../../DATA_LICENSE.md).
