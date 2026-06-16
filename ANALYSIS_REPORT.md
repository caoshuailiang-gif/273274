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

## 细胞类群注释（`annotate_celltypes.py`）

采用**基于 marker 基因的打分法**（`sc.tl.score_genes`）：对脑/胶质瘤组织 11 套经典 marker
基因集逐细胞打分，再按每个 Leiden 簇的平均分取最高者指派细胞类型。

**28 个簇 → 9 种细胞类型**：

| 细胞类型 | Leiden 簇 | 代表 marker |
|----------|-----------|-------------|
| Tumor_glioma（肿瘤细胞） | 0,1,2,5,7,10,13,14,17,24,25,26,27 | EGFR, SOX2, OLIG2, PTPRZ1, VIM |
| Astrocyte（星形胶质） | 8,16,21,23 | GFAP, AQP4, SLC1A3, ALDH1L1 |
| Neuron（神经元） | 9,11,15,19,20,22 | SNAP25, SYT1, RBFOX3, STMN2 |
| Microglia_Myeloid（小胶质/髓系） | 3 | PTPRC, AIF1, CSF1R, P2RY12, C1Q |
| Oligodendrocyte（少突胶质） | 6 | PLP1, MBP, MOG, MAG |
| Endothelial（内皮） | 4 | CLDN5, PECAM1, VWF, FLT1 |
| Pericyte_Mural（周细胞） | 18 | RGS5, PDGFRB, NOTCH3 |
| T_NK_cell（T/NK） | 12 | CD3D, CD3E, CD8A, NKG7 |

**各细胞类型在两组中的比例**：

| 细胞类型 | GC | GS |
|----------|----|----|
| Tumor_glioma | 0.393 | **0.617** |
| Astrocyte | **0.335** | 0.066 |
| Neuron | 0.170 | 0.117 |
| Microglia_Myeloid | 0.057 | 0.095 |
| Oligodendrocyte | 0.026 | 0.080 |
| T_NK_cell | 0.010 | 0.011 |
| Endothelial | 0.006 | 0.010 |
| Pericyte_Mural | 0.003 | 0.005 |

**关键差异**：GS 组**肿瘤细胞占比显著升高**（62% vs 39%），而 GC 组**星形胶质细胞占比明显更高**（34% vs 7%）；
GS 组的小胶质/髓系与少突胶质比例也更高。提示两组在肿瘤负荷与微环境组成上存在系统性差异。

**注释方法的局限**：
- 打分法将整簇指派为单一类型，未做簇内细分；若需更精细分型，可降低/提高聚类分辨率或用参考数据集做标签转移（如 CellTypist、SingleR）。
- 肿瘤细胞与正常星形胶质/OPC 存在 marker 重叠（OLIG2、VIM 等），部分边界簇的归属可结合 CNV 推断（inferCNV）进一步确认恶性状态。

## 输出文件

图（`figures/`）：
- `umap_leiden.png` / `umap_group.png` / `umap_sample.png` / `umap_qc.png`
- `violin_qc_pre.png`（过滤前 QC 小提琴图）
- `rank_genes_groups_leiden_markers.png`（各簇 marker）
- `cluster_composition_by_group.png`（簇组成柱状图）
- `umap_celltype.png` / `umap_celltype_legend.png`（细胞类型注释 UMAP）
- `dotplot__celltype_markers.png`（各类型 marker dotplot）
- `celltype_composition_by_group.png`（细胞类型组成柱状图）

表（`results/`）：
- `qc_summary_per_sample.csv`
- `cluster_markers.csv` / `top10_markers_per_cluster.csv`
- `cluster_counts_by_group.csv` / `cluster_fraction_by_group.csv`
- `DE_GS_vs_GC.csv`（GS vs GC 差异表达）
- `cluster_to_celltype.csv` / `cluster_celltype_scores.csv`（簇→细胞类型映射及打分）
- `celltype_counts_by_group.csv` / `celltype_fraction_by_group.csv`（细胞类型组成）

处理后对象：`GSE273274_processed.h5ad`（已被 `.gitignore` 忽略，未纳入版本控制）。

## 复现

```bash
git lfs pull                       # 拉取原始 10x 数据
python -m venv .venv && .venv/bin/pip install scanpy leidenalg igraph
.venv/bin/python analyze_scrnaseq.py
```
