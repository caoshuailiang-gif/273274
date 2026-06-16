# GSE273274 单细胞 RNA-seq 分析报告

## 数据概览

10x Genomics（CellRanger 5.0.0）人源单细胞数据，两个分组、各 3 个生物学重复：

| 分组 | 样本 | 原始细胞数 |
|------|------|-----------|
| GC | 87_R1_GC | 9,965 |
| GC | 89_R2_GC | 9,614 |
| GC | 91_R3_GC | 15,089 |
| GS | 86_R1_GS | 14,531 |
| GS | 88_R2_GS | 12,578 |
| GS | 90_R3_GS | 8,522 |

合并后：**70,299 细胞 × 36,601 基因**。

## 分析流程（`analyze_scrnaseq.py`）

1. 读取 6 个样本的 10x mtx 矩阵并合并，加入 `sample_id` / `group` / `replicate` 元数据
2. QC：计算每细胞基因数、总 counts、线粒体比例(`MT-`)、核糖体比例(`RPS/RPL`)
3. 过滤：`min_genes=200`、`min_cells=3`、`n_genes<7000`、`pct_mt<20%`
4. 标准化：`normalize_total(1e4)` + `log1p`
5. 高变基因 2000（按样本做 batch_key）→ 回归 total_counts/pct_mt → scale
6. PCA(30) → 近邻图 → Leiden 聚类(res=0.5) → UMAP
7. Marker 基因（Wilcoxon，每簇 vs 其余）
8. 簇 × 分组组成统计
9. GS vs GC 差异表达（以 GC 为参考）

## 主要结果

- **聚类**：识别出 **28 个 Leiden 簇**。
- **组织类型推断**：GS-vs-GC 顶部上调基因为 **PLP1、MBP、PTGDS、CTNNA3、SLC24A2** 等少突胶质/髓鞘标志物；
  部分肿瘤簇高表达 **MDM2、TERT**（胶质瘤常见扩增/激活），cluster 3 高表达 **PLXDC2、DOCK8、ARHGAP15、TBXAS1、APBB1IP**（小胶质/髓系）。
  综合提示样本来源为**脑/胶质瘤组织**。
- **分组特异性组成差异显著**：
  - GC 富集簇：21 (28.3%)、13 (12.6%)、0 (12.5%)、20 (9.4%)、24 (7.1%)
  - GS 富集簇：25 (15.6%)、2 (10.5%)、3 (9.5%)、1 (8.8%)、19 (8.4%)、6 (8.0%)
  - 多数簇几乎为单一分组所独有，说明两组细胞状态/组成差异极大（可能为不同区域、处理或肿瘤亚型）。

## 注意 / QC 提示

- 样本间深度差异较大：`91_R3_GC` 中位 counts 仅 2,295、基因数 1,559，明显低于其它样本（如 `89_R2_GC` 为 8,612 / 3,900）；测序深度差异可能影响下游聚类，必要时可考虑更严格的下限过滤或深度校正。
- 线粒体比例整体较低（中位 <2%），细胞质量良好。
- 簇与分组高度耦合，提示存在较强的批次/分组效应。若关注**跨组保守的细胞类型**，建议加入批次整合（如 Harmony / scVI），当前流程仅用 HVG 的 `batch_key` 做了弱校正。

## 输出文件

图（`figures/`）：
- `umap_leiden.png` / `umap_group.png` / `umap_sample.png` / `umap_qc.png`
- `violin_qc_pre.png`（过滤前 QC 小提琴图）
- `rank_genes_groups_leiden_markers.png`（各簇 marker）
- `cluster_composition_by_group.png`（簇组成柱状图）

表（`results/`）：
- `qc_summary_per_sample.csv`
- `cluster_markers.csv` / `top10_markers_per_cluster.csv`
- `cluster_counts_by_group.csv` / `cluster_fraction_by_group.csv`
- `DE_GS_vs_GC.csv`（GS vs GC 差异表达）

处理后对象：`GSE273274_processed.h5ad`（已被 `.gitignore` 忽略，未纳入版本控制）。

## 复现

```bash
git lfs pull                       # 拉取原始 10x 数据
python -m venv .venv && .venv/bin/pip install scanpy leidenalg igraph
.venv/bin/python analyze_scrnaseq.py
```
