# v0.6 Mac native R validation

Measured on 2026-10-07: Apple M4 Max, 16 CPU / 40 GPU cores, 128 GiB unified
memory, arm64, R 4.6.1, Matrix 1.7-5, Apple clang 17 / macOS SDK 15.5.
These are measurements of exact integer donor/group/cell-type pseudobulk.
All original features and eligible cells are retained. They do not measure GPU
training, the exact ranking search or arbitrary downstream Bioconductor models.

## Complete pseudobulk call: released v0.5 versus new native CPU

| Input | Cells × features | v0.5 median | v0.6 native median | Speedup | v0.5 peak RSS | v0.6 peak RSS | RSS reduction |
|---|---:|---:|---:|---:|---:|---:|---:|
| Fixed sparse 100k | 100,000 × 4,000 | 0.512 s | 0.169 s | 3.03× | 369.88 MiB | 288.05 MiB | 22.1% |
| Fixed sparse 250k | 250,000 × 8,000 | 3.238 s | 0.456 s | 7.10× | 647.12 MiB | 488.55 MiB | 24.5% |
| Kang real counts | 24,679 × 35,635 | 2.660 s | 0.273 s | 9.74× | 1,088.34 MiB | 730.47 MiB | 32.9% |

The synthetic cases contain 2,000,000 / 4,000,000 stored values and 64 / 128
units. Kang contains 14,188,476 stored values; its existing six missing cell-type
annotations are explicitly dropped and recorded. No new gene/cell filter was
selected after seeing performance. The reference pipeline is the publicly
released native R v0.5.0, not a dense cell-by-gene allocation.

## Mature sparse operator and Metal comparison

| Input | Matrix operator | Native operator | Metal 8 MiB operator | Metal 64 MiB operator | Native complete call | Metal 64 MiB complete call |
|---|---:|---:|---:|---:|---:|---:|
| Sparse 100k | 0.029 s | 0.011 s | 0.015 s | 0.018 s | 0.169 s | 0.178 s |
| Sparse 250k | 0.095 s | 0.029 s | 0.057 s | 0.024 s | 0.456 s | 0.439 s |
| Kang | 0.126 s | 0.087 s | budget too small | 0.075 s | 0.273 s | 0.266 s |

Matrix uses the established full sparse multiplication operator, with native
input/output validation. Membership is prepared outside operator timing. New
native/Metal operators integrate validation. Complete calls include cell-ID and
metadata checks, grouping, aggregation, filtering, splitting and provenance
hashing; operator timings do not. The first operator call, including any Metal shader/pipeline setup, is
reported separately in each JSON and excluded from warm medians.

The small Metal advantage in two cases is insufficient to establish a general
end-to-end speed advantage over the new CPU kernel. `auto` therefore selects
native CPU. A smaller staging budget adds commands and can cost time.

## Bounded Metal buffers

| Input | Whole-input staging requested buffers | Bounded requested buffers | Buffer reduction | Whole-input peak RSS | Bounded peak RSS |
|---|---:|---:|---:|---:|---:|
| Sparse 100k | 32.47 MiB | 8.00 MiB | 75.4% | 357.00 MiB | 318.75 MiB |
| Sparse 250k | 68.85 MiB | 8.00 MiB | 88.4% | 564.75 MiB | 503.98 MiB |
| Kang | 251.30 MiB | 64.00 MiB | 74.5% | 995.80 MiB | 813.27 MiB |

Whole-input staging is an explicit Metal ablation with a 1 GiB budget; it is not
v0.5 and does not represent another tool's GPU memory use. The reported buffer
sum includes persistent two-limb outputs, an overflow flag and input staging.
Kang output alone exceeds 8 MiB, so that budget is rejected as designed.

Requested buffer reductions are **not** reductions in total process memory.
Metal allocator rounding, shaders, framework caches, R objects and copies remain
outside the cap. Equal requested buffers can have different process high-water
marks (100k 64 MiB and whole-input configurations are one example). Apple unified
memory is shared; CPU RAM and GPU memory must not be summed. GPU command seconds
are recorded, but system GPU utilization was not measured.

## Correctness and measurement protocol

- 13,568 aggregated count-value comparisons across native and **actually
  executed** Metal outputs on this Mac. The repeated dense/CSC/CSR/triplet cases
  represent 91,760 input matrix positions; these are two distinct counters.
- Explicit tests for high integer limbs, simultaneous carries, multiple chunks,
  `2^53-1` acceptance, values/sums above the exact range and accumulator overflow
  past 64 bits; invalid counts in excluded annotation cells are also rejected.
- HDF5-backed row blocks agree with the independent dense reference.
- Every timed operator and complete call checks the independently generated
  Matrix count hash, including units later excluded by the min-cell rule.
- Each case/method runs in a new R process. Three warm operator measurements and
  three complete calls are retained in [the raw JSON files](mac/). Legacy v0.5
  measures complete calls only. Inputs are prepared beforehand from the fixed
  contract, and input loading is outside call timing.
- Peak RSS is the **one-process** high-water mark from macOS `/usr/bin/time -l`,
  including R startup, libraries, input, warmup, repeated calls, output and hash
  validation. It is not the incremental memory of one function. `.time.txt`
  files retain the OS measurement. RSS was not repeated across independent
  processes; the percentages describe these runs, not confidence intervals.
- All methods were run sequentially, without concurrent computation benchmarks.
  Runtime medians on a fast machine include millisecond timing quantization and
  normal system noise. Hardware-independent speed guarantees are not claimed.

The primary cases/seed/task were fixed in
[mac_gpu_contract.json](../../benchmarks/mac_gpu_contract.json) before timing.
Kang is an additional real-data check using the existing repository data and
missing-annotation policy. [mac_pseudobulk.R](../../benchmarks/mac_pseudobulk.R)
generates inputs and reports; [the guide](../../docs/MAC_CN.md) supplies commands.

The [Windows script](../../benchmarks/windows_pseudobulk.ps1) reuses the same
inputs and hash checks, with OS peak working-set sampling. Windows hardware
performance and that PowerShell runner remain for user validation. Windows CPU
package compatibility is checked separately in CI. CUDA is not implemented;
NVIDIA GPU memory, speed and utilization have not been measured.
