# ExtremaRank v0.3.0：逐特征认证、自动重拟合与单细胞扩展

本轮将此前提出的六项优化全部落实。工具的用途是回答：**删除少数独立样本或完整供者后，候选的排名入选状态及效应方向还能不能保持？** 可以接入常用差异分析流程，作为后续验证之前的敏感性检查。

| 优化 | 已实现行为 | 结果边界 |
|---|---|---|
| 逐特征稳定认证 | 分别认证每个候选的 Top-K 内／外状态和均值差方向；支持初始候选、指定 ID 或全部特征 | 整体名单会变时，仍可找到已认证稳定的候选 |
| 最小删除数 | 输出有效下界、实际反例上界；排除所有更小删除数后才标为精确最小值 | 预算之外的最小值不推断，搜索未完成保留区间 |
| 自动模型重拟合 | limma、limma-voom、edgeR、DESeq2；显式对比、分类／连续协变量、完整供者删除 | 外部模型报告已执行情景的观察结果，不冒充原生认证 |
| 供者级单细胞 | 原始计数聚合；dense、10X MTX 和可选 H5AD；显式计数层、技术重复合并声明 | 细胞不能当独立生物学重复；不自动合并随访，不自动删除不完整配对 |
| Welch 安全加速 | 有理数均值／方差界限、条件剪枝、多候选共享搜索；保留穷举参照 | 困难输入仍可能组合爆炸，到上限输出 UNRESOLVED |
| 真实案例与稳定案例 | 新增供者配对单细胞、蛋白测量、真实小鼠脂质数据及固定稳定案例 | 技术重复、百分比组成、真实生物学重复的含义分别记录 |

## 使用范围与实际价值

已验证的场景包括配对 bulk RNA-seq、独立两组微阵列、供者级单细胞 pseudobulk、技术蛋白组测量和生物学重复的脂质组。连续丰度矩阵可使用原生配对／Welch 或 limma；原始 RNA 计数可使用自动重拟合。已有其他模型或复杂分析流程可以导出完整固定特征集合的重拟合得分，再用 `compare`。

多组或时间序列需要先声明具体两组对比；原生模式不拟合协变量。未支持的复杂相关结构不能被默认为独立样本。缺失丰度默认报错，可在结果分析前明确选择完整特征排除政策；没有隐式插补。实用性来自流程适配、可解释反例与部分稳定认证，而非一个声称普适的新显著性检验。

## 真实结果

固定删除预算为 2，向上排名结果如下：

| 案例 | 删除单位 | 特征数 | 已认证稳定入选 | 已认证稳定方向 |
|---|---|---:|---:|---:|
| Golub ALL/AML 微阵列 | 38 个独立样本 | 3,051 | 6/20 | 20/20 |
| Kang B 细胞 | 8 个配对供者，删除完整 ctrl/stim 配对 | 35,635 | 6/20 | 20/20 |
| CPTAC 蛋白测量 | 6 次技术测量，不能称为 6 个患者 | 1,097 | 4/20 | 20/20 |
| Nutrimouse 脂质组 | 40 只生物学重复小鼠 | 21 | 4/5 | 5/5 |

Nutrimouse 向上 Top-5 的精确最小改变数为 **2**：所有单只小鼠删除都不改变名单，实际删除两只会改变名单。其余三个向上案例存在单单位反例。三个排名方向均执行；完整结果在 [validation_v03_real.json](../results/validation_v03_real.json)。这些结果不能说明哪个分子是真实标志物，也不构成删除影响供者的理由。

固定的 Kang count-model 特征过滤后保留 8,852 个基因。limma-voom 和 DESeq2 各完成 8 次完整供者删除重拟合。edgeR 在默认有符号 sqrt(QL F) 排名下，基线两行出现极小负原始 F，另有 6 次情景无法形成完整有限得分：保留为 **NOT_EVALUABLE**，原始行和失败说明均公开，没有截断 F、换得分或删基因来得到正面结果。另一个独立集成测试中四种模型均完整通过。

Nutrimouse 的 limma 模型明确调整五水平饮食因素，完成全部 40 次单只删除。原生 Welch 案例考察未调整的固定百分比排名；这两个目标不同，不据此比较算法优劣。

## 验证证据

- 59 项测试全部通过，包含 H5AD dense/CSR/CSC 和四种 R 模型的可选集成测试。核心不安装这些依赖也能运行，其可选测试明确跳过。
- 6,000 个原生设计／方向问题，1,434,768 项独立属性、界限与搜索上限检查，零失败。
- 183 次独立 Fraction 全特征排名，覆盖真实基线和所有不同的报告反例，零不一致；没有声称独立穷举了所有真实多删除集合。
- 4,347,470 个单细胞 pseudobulk 整数，与独立 R Matrix 求和逐个一致；122 个单位覆盖 8 种细胞类型。
- 300 个外部模型 Top-K 集合，与 R 导出三方向排序一致。
- 预先固定的分离稳定 Welch 案例：100 样本、预算 3，穷举检查 166,750 个集合，新搜索检查 24 个集合及 9 个界限节点，得到同一认证。困难平局案例在有限上限下保留 UNRESOLVED。该结果不能作为普遍提速倍数。

当前 R 验证环境为 R 4.6.1、limma 3.68.5、edgeR 4.10.5、DESeq2 1.52.0。版本、输入字节哈希和模型警告随每次拟合保留。

## 运行与复现

```bash
python -m pip install .
extremarank robustness data/prepared/Nutrimouse/matrix.tsv \
  --metadata data/prepared/Nutrimouse/metadata.tsv --target ppar --reference wt \
  --design welch --budget 2 --top-k 5 --output my-lipid-audit
```

打开输出 `report.html`，查看逐特征状态、最小删除区间及实际反例。表格和 JSON 保留全部结果；HTML 限制展示前 200 个被查询特征，避免大文件卡顿。

单细胞／重拟合命令见 [REFITS.md](REFITS.md)，证书与界限推导见 [ROBUSTNESS.md](ROBUSTNESS.md)。核心检查：

```bash
python -m unittest discover -s tests -v
python benchmarks/independent_robustness_check.py --cases 1000
python benchmarks/welch_pruning_v03.py
python data/fetch_v03_sources.py --offline
```

完整公开数据复现需 R 及上述模型包、Matrix、SummarizedExperiment、SingleCellExperiment：

```bash
python benchmarks/reproduce_v03.py
# 使用已打包准备输入，省去原始单细胞重新导出：
python benchmarks/reproduce_v03.py --checks-only
```

算法代码仍为 MIT，第三方数据的原始许可及引用独立记录在 [DATA_LICENSE.md](../DATA_LICENSE.md)。发布包包含源数据、冻结合同、执行结果和字节校验；wheel 包只含工具和 R 适配器。GitHub 公共发布方便复现和引用，不保证星标数量，也不替代独立同行审查。
