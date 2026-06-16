"""
(a) inferCNV-based malignant-cell identification for GSE273274.

Uses infercnvpy. Gene genomic positions are taken from the bundled
oligodendroglioma reference (gene order along chromosomes is conserved
hg19/hg38, which is all inferCNV needs). Clearly non-malignant cell types
(immune / endothelial / oligodendrocyte / pericyte) serve as the CNV
reference baseline; cells are then called malignant from their CNV score.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import infercnvpy as cnv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sc.settings.verbosity = 1
sc.settings.figdir = "figures"
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)
sc.settings.set_figure_params(dpi=120, frameon=False)

REFERENCE_TYPES = ["Microglia_Myeloid", "T_NK_cell", "Endothelial",
                   "Oligodendrocyte", "Pericyte_Mural"]

# ----------------------------------------------- load full-gene log-norm data
proc = sc.read_h5ad("GSE273274_processed.h5ad")
adata = proc.raw.to_adata()                 # log1p-normalized, all genes
adata.obs = proc.obs.copy()                 # leiden / cell_type / group / sample_id
print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")

# ----------------------------------------------- attach gene genomic positions
ref = cnv.datasets.oligodendroglioma()
pos = ref.var[["chromosome", "start", "end"]].copy()
pos.index = pos.index.astype(str)
adata.var["chromosome"] = adata.var_names.map(pos["chromosome"])
adata.var["start"] = adata.var_names.map(pos["start"])
adata.var["end"] = adata.var_names.map(pos["end"])
n_pos = adata.var["chromosome"].notna().sum()
print(f"Genes with genomic position: {n_pos} / {adata.n_vars}")

# ------------------------------------------------------------- run inferCNV
cnv.tl.infercnv(
    adata,
    reference_key="cell_type",
    reference_cat=REFERENCE_TYPES,
    window_size=100,
)
cnv.tl.pca(adata, n_comps=30)
cnv.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
cnv.tl.leiden(adata, resolution=0.5)
cnv.tl.cnv_score(adata)
print(f"CNV clusters: {adata.obs['cnv_leiden'].nunique()}")

# ------------------------------------------------- call malignant via CNV score
is_ref = adata.obs["cell_type"].isin(REFERENCE_TYPES)
ref_mean = adata.obs.loc[is_ref, "cnv_score"].mean()
ref_std = adata.obs.loc[is_ref, "cnv_score"].std()
threshold = ref_mean + 2 * ref_std
adata.obs["malignant"] = np.where(adata.obs["cnv_score"] > threshold,
                                  "malignant", "non-malignant")
adata.obs["malignant"] = adata.obs["malignant"].astype("category")
print(f"\nReference cnv_score: mean={ref_mean:.4f} std={ref_std:.4f} -> thr={threshold:.4f}")
print(adata.obs["malignant"].value_counts())

# ------------------------------------------------------------- summary tables
by_type = adata.obs.groupby("cell_type", observed=True).agg(
    n=("cnv_score", "size"),
    mean_cnv_score=("cnv_score", "mean"),
    pct_malignant=("malignant", lambda s: (s == "malignant").mean() * 100),
).round(3)
by_type.to_csv("results/cnv_score_by_celltype.csv")
print("\nCNV score & malignant fraction by annotated cell type:\n", by_type.to_string())

mal_group = pd.crosstab(adata.obs["malignant"], adata.obs["group"])
mal_group_frac = pd.crosstab(adata.obs["malignant"], adata.obs["group"],
                             normalize="columns").round(3)
mal_group.to_csv("results/malignant_counts_by_group.csv")
mal_group_frac.to_csv("results/malignant_fraction_by_group.csv")
print("\nMalignant fraction by group:\n", mal_group_frac.to_string())

# ----------------------------------------------------------------- figures
cnv.pl.chromosome_heatmap(adata, groupby="cell_type", show=False,
                          save="_cnv_celltype.png")
cnv.tl.umap(adata)
cnv.pl.umap(adata, color="cnv_score", show=False, save="_cnv_score.png")
cnv.pl.umap(adata, color="malignant", show=False, save="_cnv_malignant.png")
cnv.pl.umap(adata, color="cell_type", show=False, save="_cnv_celltype.png")

# ------------------------------------- persist calls back to processed object
proc.obs["cnv_score"] = adata.obs["cnv_score"].values
proc.obs["malignant"] = adata.obs["malignant"].values
proc.write_h5ad("GSE273274_processed.h5ad")
adata.obs[["sample_id", "group", "cell_type", "cnv_leiden",
           "cnv_score", "malignant"]].to_csv("results/per_cell_cnv_calls.csv")
print("\nDone. CNV figures in figures/, tables in results/.")
