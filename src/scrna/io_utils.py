"""数据读写工具。

支持多种常见单细胞输入格式，并按文件名/扩展名自动识别。
"""
from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import scanpy as sc


def load_data(path: str, transpose: bool = False) -> ad.AnnData:
    """根据路径自动识别格式并读取为 AnnData。

    支持：
        - 10x mtx 目录（含 matrix.mtx[.gz] / barcodes / features）
        - 10x h5 文件 (.h5)
        - AnnData 文件 (.h5ad)
        - csv / tsv 表达矩阵
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"输入路径不存在: {path}")

    if p.is_dir():
        # 10x mtx 目录
        adata = sc.read_10x_mtx(p, var_names="gene_symbols", cache=True)
    elif p.suffix == ".h5ad":
        adata = sc.read_h5ad(p)
    elif p.suffix == ".h5":
        adata = sc.read_10x_h5(p)
    elif p.suffix in {".csv", ".tsv", ".txt"}:
        delim = "," if p.suffix == ".csv" else "\t"
        adata = sc.read_csv(p, delimiter=delim)
        if transpose:
            adata = adata.T
    else:
        raise ValueError(f"不支持的输入格式: {p.suffix or p}")

    # 保证变量名唯一，避免后续索引报错
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    print(f"[io] 读取完成: {adata.n_obs} 个细胞 × {adata.n_vars} 个基因")
    return adata


def save_adata(adata: ad.AnnData, output_dir: str, sample_name: str, suffix: str) -> str:
    """保存 AnnData 到输出目录，返回保存路径。"""
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"{sample_name}_{suffix}.h5ad")
    adata.write_h5ad(out_path)
    print(f"[io] 已保存: {out_path}")
    return out_path
