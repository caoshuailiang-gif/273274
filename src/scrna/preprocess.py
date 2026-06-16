"""标准化与特征选择：归一化、log1p、高变基因、缩放、PCA。"""
from __future__ import annotations

import anndata as ad
import scanpy as sc


def normalize(adata: ad.AnnData, target_sum: float = 1e4, log1p: bool = True) -> ad.AnnData:
    """文库大小归一化 + 可选 log1p。原始计数保留在 layers['counts']。"""
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    if log1p:
        sc.pp.log1p(adata)
    return adata


def select_hvg(adata: ad.AnnData, n_top_genes: int = 2000,
               flavor: str = "seurat_v3") -> ad.AnnData:
    """选择高变基因 (HVG)。

    seurat_v3 在原始计数上计算，因此使用 layers['counts']。
    """
    kwargs = dict(n_top_genes=n_top_genes, flavor=flavor)
    if flavor == "seurat_v3":
        kwargs["layer"] = "counts"
    sc.pp.highly_variable_genes(adata, **kwargs)
    n_hvg = int(adata.var["highly_variable"].sum())
    print(f"[preprocess] 选出 {n_hvg} 个高变基因")
    return adata


def scale(adata: ad.AnnData, regress_out: bool = False,
          max_value: float = 10) -> ad.AnnData:
    """（可选回归混杂变量后）按基因标准化到单位方差。

    完整表达谱保存在 raw，下游缩放仅作用于 HVG 子集。
    """
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()
    if regress_out:
        sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(adata, max_value=max_value)
    return adata


def run_pca(adata: ad.AnnData, n_comps: int = 50, random_state: int = 0) -> ad.AnnData:
    """主成分分析。"""
    n_comps = min(n_comps, adata.n_vars - 1, adata.n_obs - 1)
    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack", random_state=random_state)
    return adata


def run_preprocess(adata: ad.AnnData, cfg: dict, random_state: int = 0) -> ad.AnnData:
    """执行完整预处理流程，由配置驱动。"""
    adata = normalize(adata, cfg["target_sum"], cfg["log1p"])
    adata = select_hvg(adata, cfg["n_top_genes"], cfg["hvg_flavor"])
    adata = scale(adata, cfg.get("regress_out", False), cfg["scale_max_value"])
    adata = run_pca(adata, random_state=random_state)
    return adata
