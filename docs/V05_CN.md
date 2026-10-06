# v0.5：计算范围、候选结果与大数据工作流

这是独立原生 R 包。所有接口直接在 R 会话中运行；不需要 Python。

## 安装

R 4.6 的 Windows x86_64、macOS arm64 和 macOS x86_64 使用 GitHub Release
的预编译包（macOS 为 CRAN R.framework 构建）。其他受支持的 R 版本或
Homebrew R 可安装源码包，需要 C++17 编译器。
在 R 中运行版本固定的安装辅助函数，它会核对下载文件的 SHA-256：

```r
source("https://raw.githubusercontent.com/heise3/ExtremaRank/v0.5.0/r/install_extremarank.R")
install_extremarank()
library(extremarank)
```

或显式安装源码：

```r
install.packages(c("Rcpp", "BH", "Matrix", "digest"))
install.packages(
  "https://github.com/heise3/ExtremaRank/releases/download/v0.5.0/extremarank_0.5.0.tar.gz",
  repos = NULL, type = "source")
```

## 执行前与执行后的范围

```r
# counts 的行是特征、列是样本；samples 包含 sample_id/group/donor_id。
preflight_extremarank(counts, samples, "treated", "control", budget = 2)
result <- extremarank(counts, samples, "treated", "control", budget = 2,
  transform = "cpm-log2", search_order = "influence", max_diagnostics = 50)
summary(result)
as.data.frame(result)
plot(result)
plot(result, type = "influence")
```

预检不拟合模型，返回精确可行删除数和粗略数值存储估计。估计不包含模型
内部对象、ID、任意精度整数和全部工作者内存，不应解释成峰值内存保证。

默认诊断在预算 0 时跳过，其余预算最多执行 100 次 LOO。诊断有独立上限，
记录未运行单位；`diagnostics="all"` 可在预算 0 时显式开启，并仍使用
`max_diagnostics`。`result$timing` 分别记录搜索和诊断耗时。

## 常用模型的逐候选结果

```r
fits <- refit_extremarank(counts, samples, "treated", "control",
  model = "DESeq2", design = "paired", k = 20, budget = 1,
  store_scores = "disk", score_path = "model-scores",
  checkpoint = "fit-cache", workers = 2, progress = TRUE)
summary(fits)
fits$coverage
fits$features
fits$failed_features
```

- `NOT_EVALUATED`：未运行任何删除；只有基线不能称已观察稳定。
- `PARTIALLY_EVALUATED`：尚未找到变化，但存在失败或因上限未运行的情景。
- `OBSERVED_CHANGED`：确实观察到变化；其他失败或未运行情景仍保留。
- `OBSERVED_STABLE`：声明范围内的可行情景全部运行且可评估，未观察到变化。
- `NOT_EVALUABLE`：基线排序不可评估。

整体 status 描述整个 Top-K 的成员变化；方向状态单独报告，汇总中分别列出。
逐候选结果包括成员和方向变化次数、最小已观察删除数、供者或块 ID 反例，
以及包含基线的已观察最好/最差排名。选择比例仅针对有效的已执行删除，不是概率。
`failed_features` 指明每个情景中的非有限分数/效应特征。它们不是未运行情景
的证书；最小已观察删除数也不自动等于真实最小删除数。

`store_scores="memory"` 保留完整分数；`"disk"` 逐情景写 RDS；`"none"`
只保留摘要、候选反例和失败特征 ID，主动丢弃完整分数。导出磁盘结果时：

```r
write_extremarank(fits, "exported-fit")
fits2 <- read_extremarank("exported-fit")
```

磁盘分数复制到导出目录并带校验值，整个目录可移动。重读会验证校验值。

## 断点恢复、并行、对比和整组删除

再次调用相同模型与固定参数，设置 `checkpoint="fit-cache", resume=TRUE`。
可增加 `max_refits`。输入、预算、模型/包版本、过滤、候选集、对比或随机种子
改变，会拒绝复用旧缓存。断点恢复针对模型重拟合，不是原生 DFS 栈恢复。

`workers` 可为 1–4。每个情景有固定 RNG 初始化，调用者的随机状态被恢复；
单线程与并行按同一情景顺序收集。模型依赖外部状态的行为仍须由用户控制。

```r
# budget 计数整个 site 块；一个供者不能跨两个块。
sites <- refit_extremarank(counts, samples, "treated", "control",
  model = "limma-voom", design = "paired", delete_by = "site", budget = 1)
# 预声明单自由度对比；对比系数名是固定 model.matrix 的列名。
reverse <- refit_extremarank(expression, samples, "treated", "control",
  model = "limma", contrast = c(grouptarget = -1))
```

整组计划用精确整数动态规划计数；只枚举删除后满足最少重复数的情景。
使用固定 treatment 编码，不受全局 options(contrasts=...) 影响。
无法识别的协变量、消失的非零对比系数或改变的对比参照水平作为失败保存。
`refit_extremarank_custom()` 接受用户声明的拟合函数与整供者删除单位，可接入
多条件或重复测量模型。函数必须重新拟合收到的保留样本，并返回完整固定
特征集合的 `feature_id/score/effect`。外部依赖和模型设置必须在分析前固定。

## 磁盘单细胞输入与内存

支持 DelayedMatrix/HDF5Matrix 输入和实验对象中的磁盘 assay。原始细胞矩阵
分块验证和汇总，不整体转换为稠密矩阵。较小的特征乘供者输出会被物化。
`block_rows` 控制块大小，输入块自动调整到约 64 MiB；一行本身更大时不能
进一步拆行。`max_dense_bytes` 默认 512 MiB，限制准备后的稠密矩阵、汇总
输出和数值输入的工作者副本；不覆盖全部模型工作内存。

## 算法改进与实际增益

特征排序在一个审计对象中预先缓存，条件搜索复用其排序；影响度启发式只
改变搜索顺序，所有排序决策和剪枝证明仍用精确计算。Welch 方差界加入了
固定值约束下的成对差平方分解：固定—固定、固定—可选、可选—可选。
分别取保守下界，再与已有下界取较大值；不把不能同时达到的端点称为精确
联合极值。参见 [推导](V05_BOUNDS.md)。

`robustness_curve(effects, budgets=0:2)` 对预先计算的供者效应生成预算曲线，
明确区分已证明、找到反例和未解决。随机重抽样和有限 LOO 不能替代穷尽证书。
真实数据与受控验证结果见 [v0.5 验证](../results/v05/README.md)。
