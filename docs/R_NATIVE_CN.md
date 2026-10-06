# 原生 R 包：直接在 R 中运行 ExtremaRank

v0.5.0 提供独立的标准 R 包，源码位于 `r/extremarank`。算法在 R 会话中
直接调用编译后的 C++，没有 Python、reticulate、Rscript 子进程或命令行转接。
这是 R 中直接可调用的算法实现。原有 Python 版本继续保留。

v0.5 的预编译安装、覆盖率、逐候选重拟合、磁盘输入、恢复和绘图见
[新增实用功能](V05_CN.md)。

## 安装和第一个完整例子

```r
install.packages(c("Rcpp", "BH", "Matrix", "digest"))
install.packages(
  "https://github.com/heise3/ExtremaRank/releases/download/v0.5.0/extremarank_0.5.0.tar.gz",
  repos = NULL, type = "source")
library(extremarank)

# 预先计算的供者内效应：行是供者，列是特征。
effects <- matrix(c(5, 4, 6, 5, 1, 2, 1, 2, -2, -1, -3, -2), nrow = 4,
  dimnames = list(paste0("donor", 1:4), c("A", "B", "C")))
result <- extremarank_effects(effects, k = 1, budget = 1)
result
result$features[, c("feature_id", "membership_status", "direction_status")]
write_extremarank(result, "r-audit-output")
```

运行不需要 Python。R 包依赖 Rcpp、BH、Matrix 和 digest。
从源码安装需要 C++17 编译器：Windows 安装与 R 版本对应的 Rtools；
macOS 安装 Xcode 命令行工具。GitHub Release 还提供小型标准 R 源码包
`extremarank_0.5.0.tar.gz`，可用 `install.packages(..., repos = NULL, type = "source")`
安装。大型公共示例数据位于完整仓库归档，不放进 R 安装包。
上面的安装方法只下载小型 R 包。也可以安装 `remotes` 后使用
`remotes::install_github("heise3/ExtremaRank", subdir = "r/extremarank", ref = "v0.5.0")`。

目前发布渠道是 GitHub；本版本没有宣称已经进入 CRAN 或 Bioconductor。

## 实际表达矩阵与生信对象

普通矩阵的行是基因、蛋白或其他特征，列是样本。元数据包含
`sample_id`、`group`，配对研究另需 `donor_id`。样本 ID 的集合必须完全匹配，
包会按矩阵列顺序对齐，保留因子、前导零和 Unicode 标识符。

```r
result <- extremarank(expression, samples,
  target = "treated", reference = "control",
  design = "paired", k = 20, budget = 2,
  transform = "cpm-log2")
result$features
result$audit$topk_minimum_change
result$diagnostics$features
```

原始 counts 可显式选择 `cpm-log2`；已准备好的连续表达量、蛋白丰度或
脂质测量通常使用 `transform = "none"`。独立两组改为 `design = "welch"`。
处理 `SummarizedExperiment` 时可省略元数据，从 `colData` 读取，但必须显式
选择 `assay`。原生证书针对 R 完成预处理后冻结的 binary64 数值；不同软件
的对数和文本解析可能存在末位差异，不宣称所有跨语言预处理都逐位相同。

默认检查原始 Top-K 的每个候选。用 `features = c("GENE1", "GENE2")` 指定候选，
或用 `all_features = TRUE` 检查所有特征。未指定的特征仍参与排序竞争。
所有反例索引使用 R 的 **1 起始索引**，同时附带确切的供者或样本 ID。

## 在 R 内重拟合常用模型

```r
install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "DESeq2"))

fits <- refit_extremarank(counts, samples, "treated", "control",
  model = "DESeq2", design = "paired", k = 20, budget = 1,
  min_count_samples = 3)
fits$fit_status
fits$scores$baseline
fits$scenarios
```

支持 limma、limma-voom、edgeR 和 DESeq2。limma 接收连续测量，其余模型接收
原始整数 counts。独立组可声明 `categorical = "batch"`、`numeric = "age"`。
包直接调用相应 R 函数；每次删除后重新估计标准化与方差或离散度。
基线过滤、协变量、分析对比和特征集合保持固定。若删除造成设计不可辨识、
剩余自由度不足或排序分数非有限，保留 `NOT_EVALUABLE` 记录、原始分数、
警告和模型信息，不在看到结果后删基因或截断统计量。

模型重拟合结果属于已运行情景的观察敏感性，不具备原生配对/Welch 证书的
数学保证。`max_refits` 是明确的情景上限，达到后标记 `enumeration_capped`。

## 单细胞供者级分析

```r
pb <- pseudobulk_extremarank(sce, assay = "counts", min_cells = 10)
b <- pseudobulk_experiment(pb, "B cells")
result <- extremarank(b, target = "stim", reference = "ctrl",
  assay = "counts", design = "paired", transform = "cpm-log2",
  k = 20, budget = 2)
```

支持 `SingleCellExperiment`、稀疏 Matrix 和普通原始 counts 矩阵。
细胞元数据需要 `cell_id`、`sample_id`、`donor_id`、`group`、`cell_type`。
汇总按供者、条件和细胞类型进行，删除单位是整个供者。直接把单细胞对象
送入样本删除审计会被拒绝，以免将细胞数误当成生物学重复数。

不同样本 ID 不会被默认合并；确认属于技术重复时才显式启用
`merge_technical = TRUE`。缺少细胞类型注释时默认报错，可预先声明
`missing_cell_type = "drop"`，被排除的细胞会列在 `excluded_cells`。
低细胞数的汇总单位列在 `excluded_units`。汇总过程不将完整稀疏细胞矩阵
转成稠密矩阵，仅生成较小的特征乘供者输出，并检查整数精确表示范围。

## 结果含义与验证

- `CERTIFIED`：预算范围内所有可行删除都保持指定的成员身份或效应符号。
- `REFUTED`：找到并重新计算了实际共同删除反例。
- `UNRESOLVED`：在当前搜索上限下不能证明，但保留有效的最少删除数区间。
- 模型重拟合使用 `OBSERVED_STABLE`、`OBSERVED_CHANGED`、
  `PARTIALLY_EVALUATED`、`NOT_EVALUATED` 或 `NOT_EVALUABLE`，不认证未运行的删除组合。

R 的包测试独立运行，不调用 Python。测试包括完整删除枚举、强制保留值的
极值、搜索上限、Unicode/前导零 ID、稀疏汇总、实验对象与四个模型。
额外的开发验证用 Python Fraction 完整枚举核对 R 原生算法，包括消去误差、
极大/极小值与次正规数；这属于验证工具，不是 R 包运行依赖。

具体结果、复现命令和失败情景见 [v0.4 验证记录](../results/v04/README.md)。
本包证明的是候选排序对删样本的敏感性，不是 FDR 控制、生物标志物有效性、
因果效应或临床效用。
