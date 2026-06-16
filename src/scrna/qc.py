"""质量控制：计算 QC 指标、过滤低质量细胞与基因、可选双细胞检测。"""
from __future__ import annotations

import anndata as ad
import numpy as np
import scanpy as sc


def calculate_qc_metrics(adata: ad.AnnData, mito_prefix: str = "MT-",
                         ribo_prefix=("RPS", "RPL")) -> ad.AnnData:
    """标记线粒体/核糖体基因并计算每个细胞的 QC 指标。"""
    adata.var["mt"] = adata.var_names.str.startswith(mito_prefix)
    adata.var["ribo"] = adata.var_names.str.startswith(tuple(ribo_prefix))

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    print(f"[qc] 中位基因数={np.median(adata.obs['n_genes_by_counts']):.0f}, "
          f"中位线粒体比例={np.median(adata.obs['pct_counts_mt']):.2f}%")
    return adata


def detect_doublets(adata: ad.AnnData, random_state: int = 0) -> ad.AnnData:
    """使用 scrublet 检测双细胞，结果写入 obs['predicted_doublet']。

    scrublet 缺失时跳过并给出提示，不中断流程。
    """
    try:
        sc.pp.scrublet(adata, random_state=random_state)
        n_doublets = int(adata.obs.get("predicted_doublet", []).sum())
        print(f"[qc] scrublet 检测到 {n_doublets} 个潜在双细胞")
    except (ImportError, ModuleNotFoundError):
        print("[qc] 未安装 scrublet，跳过双细胞检测")
    return adata


def filter_cells_and_genes(adata: ad.AnnData, min_genes: int = 200,
                           min_cells: int = 3, max_genes: int = 6000,
                           max_pct_mito: float = 20.0,
                           remove_doublets: bool = True) -> ad.AnnData:
    """按阈值过滤细胞与基因。返回过滤后的 AnnData。"""
    n_start = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    keep = (adata.obs["n_genes_by_counts"] < max_genes) & \
           (adata.obs["pct_counts_mt"] < max_pct_mito)
    if remove_doublets and "predicted_doublet" in adata.obs:
        keep &= ~adata.obs["predicted_doublet"].astype(bool)

    adata = adata[keep].copy()
    print(f"[qc] 过滤: {n_start} -> {adata.n_obs} 个细胞 "
          f"(保留 {adata.n_obs / n_start:.1%})")
    return adata


def run_qc(adata: ad.AnnData, cfg: dict) -> ad.AnnData:
    """执行完整 QC 流程，由配置驱动。"""
    adata = calculate_qc_metrics(adata, cfg["mito_prefix"], cfg["ribo_prefix"])
    if cfg.get("detect_doublets", False):
        adata = detect_doublets(adata)
    adata = filter_cells_and_genes(
        adata,
        min_genes=cfg["min_genes"],
        min_cells=cfg["min_cells"],
        max_genes=cfg["max_genes"],
        max_pct_mito=cfg["max_pct_mito"],
        remove_doublets=cfg.get("detect_doublets", False),
    )
    return adata
