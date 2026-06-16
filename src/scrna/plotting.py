"""图表输出：QC、降维聚类、标记基因可视化。"""
from __future__ import annotations

import os

import anndata as ad
import matplotlib

matplotlib.use("Agg")  # 无界面环境下保存图片
import scanpy as sc


def _setup(figures_dir: str) -> None:
    os.makedirs(figures_dir, exist_ok=True)
    sc.settings.figdir = figures_dir
    sc.settings.autoshow = False
    sc.set_figure_params(dpi=120, dpi_save=300, frameon=False)


def plot_qc(adata: ad.AnnData, figures_dir: str) -> None:
    """QC 小提琴图与散点图。"""
    _setup(figures_dir)
    sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                 jitter=0.4, multi_panel=True, save="_qc.png", show=False)
    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt",
                  save="_qc_mito.png", show=False)


def plot_embedding(adata: ad.AnnData, figures_dir: str, color=("leiden",)) -> None:
    """UMAP 嵌入图。"""
    _setup(figures_dir)
    sc.pl.umap(adata, color=list(color), legend_loc="on data",
               save="_clusters.png", show=False)


def plot_markers(adata: ad.AnnData, figures_dir: str, n_genes: int = 5) -> None:
    """各簇 top 标记基因点图。"""
    _setup(figures_dir)
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=n_genes,
                                    save="_markers.png", show=False)


def run_plots(adata: ad.AnnData, figures_dir: str) -> None:
    """生成全套图表（缺失字段时跳过对应图）。"""
    if "pct_counts_mt" in adata.obs:
        plot_qc(adata, figures_dir)
    if "X_umap" in adata.obsm:
        plot_embedding(adata, figures_dir)
    if "rank_genes_groups" in adata.uns:
        plot_markers(adata, figures_dir)
    print(f"[plotting] 图表已保存至 {figures_dir}/")
