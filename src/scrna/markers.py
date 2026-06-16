"""标记基因 / 差异表达分析。"""
from __future__ import annotations

import os

import anndata as ad
import pandas as pd
import scanpy as sc


def rank_marker_genes(adata: ad.AnnData, groupby: str = "leiden",
                      method: str = "wilcoxon") -> ad.AnnData:
    """对每个簇做一对其余的差异表达，找出标记基因。"""
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, use_raw=True)
    print(f"[markers] 已用 {method} 计算各簇标记基因")
    return adata


def export_markers(adata: ad.AnnData, output_dir: str, sample_name: str,
                   n_genes: int = 25) -> pd.DataFrame:
    """导出每个簇 top 标记基因为长表 CSV，返回该表。"""
    os.makedirs(output_dir, exist_ok=True)
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    rows = []
    for grp in groups:
        for rank in range(min(n_genes, len(result["names"][grp]))):
            rows.append({
                "cluster": grp,
                "rank": rank + 1,
                "gene": result["names"][grp][rank],
                "logfoldchange": result["logfoldchanges"][grp][rank],
                "pval_adj": result["pvals_adj"][grp][rank],
                "score": result["scores"][grp][rank],
            })
    df = pd.DataFrame(rows)
    out_path = os.path.join(output_dir, f"{sample_name}_markers.csv")
    df.to_csv(out_path, index=False)
    print(f"[markers] 已导出标记基因: {out_path}")
    return df


def run_markers(adata: ad.AnnData, cfg: dict, output_dir: str,
                sample_name: str) -> ad.AnnData:
    """执行标记基因流程，由配置驱动。"""
    adata = rank_marker_genes(adata, method=cfg["method"])
    export_markers(adata, output_dir, sample_name, cfg["n_genes"])
    return adata
