# Executed artifacts

Large executed TSV tables are stored as lossless `.tsv.gz` files to keep the
repository compact. `results/archive_manifest.json` records the original and
compressed SHA256, sizes, and exact decompression verification. The benchmark
JSON records the original executed output names and byte sizes; those tables
can be restored by ordinary gzip decompression. Reproduction writes uncompressed
TSV files, as documented by the CLI. No rows were discarded or rounded during
distribution compression.

- `results/real/`: complete per-gene original-count benchmarks, practical audit,
  exact-cardinality extrema, matched comparison and family ablation failures.
- `results/GSE87290_cli/`, `results/GSE50760_cli/`: public CLI results with exact
  rational moments, donor witnesses, summary and joint audit JSON.
- `results/real_topk.json`: all six real shared-ranking counterexamples and
  independent baseline/witness recomputation status.
- `results/cli_end_to_end.json`: exact inputs and full CLI execution times.
- `results/independent_theorem_check.json`: 200,000 independent exact cases.
- `results/search_ablation.json`: eight declared synthetic joint searches.
- `results/scaling.json`: declared synthetic scalar scaling, with large
  exhaustive counts explicitly distinguished from executed enumeration.
- `results/report.html`: offline interactive display, no external assets.
- `results/candidate_scaling.svg`: standalone editable computational figure.

CSV/gzip input hashes and original source hashes are in `data/`; their
computation target is defined by `data/prepared/*/preparation.json`.
