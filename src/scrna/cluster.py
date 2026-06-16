"""降维与聚类：批次整合、邻接图、UMAP、Leiden。"""
from __future__ import annotations

import anndata as ad
import scanpy as sc


def integrate_harmony(adata: ad.AnnData, batch_key: str, n_pcs: int = 50) -> str:
    """用 Harmony 做批次整合，返回供下游使用的表征字段名。

    harmonypy 缺失或 batch_key 不存在时回退到原始 PCA。
    """
    if batch_key not in adata.obs:
        print(f"[cluster] 未找到批次列 '{batch_key}'，跳过整合")
        return "X_pca"
    try:
        sc.external.pp.harmony_integrate(adata, key=batch_key)
        print(f"[cluster] 已用 Harmony 按 '{batch_key}' 整合批次")
        return "X_pca_harmony"
    except (ImportError, ModuleNotFoundError):
        print("[cluster] 未安装 harmonypy，跳过批次整合")
        return "X_pca"


def neighbors_umap(adata: ad.AnnData, use_rep: str = "X_pca", n_pcs: int = 50,
                   n_neighbors: int = 15, random_state: int = 0) -> ad.AnnData:
    """构建邻接图并计算 UMAP 嵌入。"""
    npcs = None if use_rep != "X_pca" else n_pcs
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=npcs,
                    use_rep=use_rep, random_state=random_state)
    sc.tl.umap(adata, random_state=random_state)
    return adata


def leiden(adata: ad.AnnData, resolution: float = 1.0, random_state: int = 0) -> ad.AnnData:
    """Leiden 社区发现聚类，结果写入 obs['leiden']。"""
    sc.tl.leiden(adata, resolution=resolution, random_state=random_state,
                 flavor="igraph", n_iterations=2, directed=False)
    n_clusters = adata.obs["leiden"].nunique()
    print(f"[cluster] Leiden (分辨率={resolution}) 得到 {n_clusters} 个簇")
    return adata


def run_cluster(adata: ad.AnnData, cfg: dict) -> ad.AnnData:
    """执行完整降维聚类流程，由配置驱动。"""
    rs = cfg.get("random_state", 0)
    use_rep = "X_pca"
    if cfg.get("batch_correction", "none") == "harmony":
        use_rep = integrate_harmony(adata, cfg["batch_key"], cfg["n_pcs"])
    adata = neighbors_umap(adata, use_rep=use_rep, n_pcs=cfg["n_pcs"],
                           n_neighbors=cfg["n_neighbors"], random_state=rs)
    adata = leiden(adata, resolution=cfg["leiden_resolution"], random_state=rs)
    return adata
