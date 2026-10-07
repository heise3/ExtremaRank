# Mac 原生 R 优化与 Windows 复现

v0.6 的优化对象是从原始细胞计数生成 donor/group/cell-type pseudobulk。
细胞数、基因数、计数、分组、过滤规则与 v0.5 保持一致。后续精确排序证明
仍使用原有任意精度内核；这次没有实现单细胞深度模型训练的显存优化。

## 在 R 中使用

```r
source("https://raw.githubusercontent.com/heise3/ExtremaRank/v0.6.0/r/install_extremarank.R")
install_extremarank()
library(extremarank)
extremarank_backend_info()

# 新的 CPU 单遍内核是默认值。普通 Windows、Linux、Mac 都可用。
pb <- pseudobulk_extremarank(counts, cells, backend = "native")

# Mac 可选 Metal；整个调用留在 R 会话内。
if (isTRUE(extremarank_backend_info()$metal_available)) {
    pb_gpu <- pseudobulk_extremarank(counts, cells, backend = "metal",
                                   max_gpu_bytes = 64 * 1024^2)
    stopifnot(identical(pb$data, pb_gpu$data))
    pb_gpu$provenance$computation
}
```

`counts` 的行为基因、列为细胞；`cells` 有 `cell_id/sample_id/donor_id/group/cell_type`。
SingleCellExperiment 可直接传入，但需要 `assay="counts"`。稠密矩阵、CSC、CSR、
triplet 稀疏矩阵可直接扫描，HDF5Array/DelayedArray 保持按行分块读取。
新 CPU 内核避免反复切稀疏矩阵和为全体存储值建立 R 校验临时向量，
元数据用稳定的整数三元组分组。

Metal 使用共享内存、整数高低位累加和进位检查。它按预算分块填充输入缓冲区，
最多两个缓冲区交替使用；不丢细胞、不改成浮点近似。输入和每个汇总值必须在
`0 .. 2^53-1` 之内，超限拒绝执行。首次调用会在当前 R 会话编译 shader，
后续调用复用计算管线。预编译 R 包使用系统 Metal/Foundation，无需 Python。
源码安装需要 macOS Command Line Tools 的 C++/Objective-C++ 编译器。

## 如何理解内存和速度

`max_gpu_bytes` 限制本次请求的 **Metal 数据缓冲区字节数**，涵盖整数输出、
溢出标志和输入暂存区。完整矩阵的最低预算约为 `8 * 基因数 * 汇总单元数 + 36`。
预算不足会明确报错。磁盘输入逐块计算时，此最低预算以当前输入块的基因数计算。
该参数不限制 R 对象、输出副本、Metal 框架、分配器按页取整等开销。

Apple Silicon 使用 CPU/GPU 共享的统一内存，不能把 RAM 和“显存”相加。
报告分别给出进程峰值 RSS、请求的 Metal 缓冲区和整条 pseudobulk 调用耗时。
`gpu_command_seconds` 是 GPU 命令时间，不能据此计算系统 GPU 利用率。
Metal 需要输入打包和命令提交，GPU 核心计算快也不保证完整 R 调用更快。
`backend="auto"` 当前选择 `native`；`metal` 为显式选项。

实际数据、哈希、每轮测量和比较见 [v0.6 验证报告](../results/v06/README.md)。
基准契约在计时前固定于 [mac_gpu_contract.json](../benchmarks/mac_gpu_contract.json)。
主基准是 10 万和 25 万细胞；Kang 真实数据为额外验证，沿用已有输入和缺失标注政策。

## Mac 复现

先在两个独立 library 安装 v0.5 和 v0.6，安装函数的 `version`、`lib` 参数可指定。
在仓库根目录运行以下命令；基准脚本需要 Matrix、digest、jsonlite。

```sh
Rscript benchmarks/mac_pseudobulk.R prepare sparse-100k /tmp/sparse-100k.rds
/usr/bin/time -l Rscript benchmarks/mac_pseudobulk.R run /tmp/sparse-100k.rds \
  native /tmp/native.json /path/to/v06-library
/usr/bin/time -l Rscript benchmarks/mac_pseudobulk.R run /tmp/sparse-100k.rds \
  released-0.5 /tmp/v05.json /path/to/v05-library
/usr/bin/time -l Rscript benchmarks/mac_pseudobulk.R run /tmp/sparse-100k.rds \
  metal-8MiB /tmp/metal.json /path/to/v06-library
```

`sparse-250k`、`kang-real` 是其他输入。方法还包括 `matrix-operator`、`metal-64MiB`
和 `metal-whole`（最多 1 GiB 暂存预算的消融）。每个方法使用新进程，输入加载后
先预热，再运行三次算子和三次完整调用，并检查独立 Matrix 汇总哈希。
进程 RSS 包括 R、依赖、输入、预热、输出及哈希校验；不能解释成某一个函数的
独占内存。Kang 的整数输出超过 8 MiB，因此该数据使用 64 MiB 或更高预算。

## Windows 后续验证

Windows 原生 CPU 包和相同精确性测试继续支持。v0.6 的 GPU 路径是 Mac Metal，
`cuda_available` 返回 `FALSE`；NVIDIA CUDA 执行和硬件性能验证仍待后续实现。
以下脚本用于 Windows 上比较新 CPU 与 v0.5，并记录进程峰值工作集，
不会把 CPU 执行报告成 NVIDIA 加速。

```r
source("r/install_extremarank.R")
dir.create("C:/er-libs/v06", recursive=TRUE, showWarnings=FALSE)
dir.create("C:/er-libs/v05", recursive=TRUE, showWarnings=FALSE)
install_extremarank("0.6.0", lib="C:/er-libs/v06")
install_extremarank("0.5.0", lib="C:/er-libs/v05")
install.packages("jsonlite")
```

在仓库根目录用 PowerShell：

```powershell
./benchmarks/windows_pseudobulk.ps1 -CurrentLibrary C:/er-libs/v06 -LegacyLibrary C:/er-libs/v05
```

输出保存到 `results/windows-local/`。脚本复用同一 R 输入生成和哈希比较流程，
并记录 R/Matrix 版本；结果只反映执行机器，不跨硬件直接比较绝对速度。
