"""
(b) Harmony integration of non-malignant cells, re-clustering & re-annotation.

inferCNV (step a) found no clear aneuploid/malignant population, so essentially
all cells are non-malignant. We exclude the few CNV-flagged cells, integrate the
rest across the 6 samples with Harmony, then re-cluster and re-annotate with
canonical normal-brain markers (the previous "Tumor_glioma" label is replaced by
neural-progenitor / radial-glia markers, since the CNV signal does not support
malignancy).
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sc.settings.verbosity = 1
sc.settings.figdir = "figures"
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)
sc.settings.set_figure_params(dpi=120, frameon=False)

marker_sets = {
    "Oligodendrocyte":          ["PLP1", "MBP", "MOG", "MOBP", "MAG", "CLDN11", "CNP"],
    "OPC":                      ["PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10", "BCAN"],
    "Astrocyte":                ["GFAP", "AQP4", "SLC1A3", "ALDH1L1", "S100B", "GJA1"],
    "NeuralProgenitor_RadialGlia": ["SOX2", "NES", "VIM", "HES1", "FABP7", "HOPX",
                                    "PTPRZ1", "PAX6"],
    "Neuron_Excitatory":        ["SLC17A7", "SATB2", "RBFOX3", "SNAP25", "NEFL", "STMN2"],
    "Neuron_Inhibitory":        ["GAD1", "GAD2", "DLX1", "DLX2"],
    "Microglia_Myeloid":        ["PTPRC", "AIF1", "CSF1R", "P2RY12", "CX3CR1",
                                 "C1QA", "C1QB", "C1QC", "CD68", "TYROBP"],
    "T_NK_cell":                ["CD3D", "CD3E", "CD2", "CD8A", "IL7R", "NKG7", "GNLY"],
    "Endothelial":              ["CLDN5", "PECAM1", "VWF", "FLT1", "A2M", "ITM2A"],
    "Pericyte_Mural":           ["RGS5", "PDGFRB", "ACTA2", "NOTCH3", "MYH11"],
    "Proliferating":            ["MKI67", "TOP2A", "CENPF", "CCNB1", "UBE2C"],
}

# ----------------------------------------------------- load + select cells
proc = sc.read_h5ad("GSE273274_processed.h5ad")
adata = proc.raw.to_adata()                 # log1p-normalized, all genes
adata.obs = proc.obs.copy()
adata = adata[adata.obs["malignant"] == "non-malignant"].copy()
print(f"Non-malignant cells: {adata.n_obs} x {adata.n_vars}")

# ------------------------------------------ HVG -> scale -> PCA -> Harmony
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample_id",
                            flavor="seurat")
adata.raw = adata
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30)

import harmonypy
ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, ["sample_id"])
Z = np.asarray(ho.Z_corr)
if Z.shape[0] != adata.n_obs:          # harmonypy returns (n_pcs, n_cells)
    Z = Z.T
adata.obsm["X_pca_harmony"] = Z

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep="X_pca_harmony")
sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(adata)
print(f"Leiden clusters after integration: {adata.obs['leiden'].nunique()}")

# ----------------------------------------------- re-annotate (marker scoring)
full = adata.raw.to_adata()
full.obs = adata.obs.copy()
score_cols = []
for ct, genes in marker_sets.items():
    present = [g for g in genes if g in full.var_names]
    sc.tl.score_genes(full, present, score_name=f"score_{ct}", ctrl_size=50)
    score_cols.append(f"score_{ct}")

mean_scores = full.obs.groupby("leiden", observed=True)[score_cols].mean()
mean_scores.columns = [c.replace("score_", "") for c in mean_scores.columns]
cluster_to_ct = mean_scores.idxmax(axis=1)
mean_scores["assigned"] = cluster_to_ct
mean_scores.to_csv("results/harmony_cluster_celltype_scores.csv")
adata.obs["cell_type_int"] = adata.obs["leiden"].map(cluster_to_ct).astype("category")
full.obs["cell_type_int"] = adata.obs["cell_type_int"].values
print("\nCluster -> cell type (integrated):")
print(cluster_to_ct.to_string())

pd.DataFrame({"leiden": cluster_to_ct.index, "cell_type": cluster_to_ct.values,
              "n_cells": adata.obs["leiden"].value_counts().reindex(cluster_to_ct.index).values
              }).to_csv("results/harmony_cluster_to_celltype.csv", index=False)

# --------------------------------------------------- composition by group
ct_frac = pd.crosstab(adata.obs["cell_type_int"], adata.obs["group"], normalize="columns")
pd.crosstab(adata.obs["cell_type_int"], adata.obs["group"]).to_csv(
    "results/harmony_celltype_counts_by_group.csv")
ct_frac.to_csv("results/harmony_celltype_fraction_by_group.csv")
print("\nIntegrated cell-type fraction by group:\n", ct_frac.round(3).to_string())

# ---- quantify sample mixing: per-cluster dominant-sample fraction (vs before)
dom = pd.crosstab(adata.obs["leiden"], adata.obs["sample_id"],
                  normalize="index").max(axis=1)
print(f"\nDominant-sample fraction per cluster (median) AFTER Harmony: {dom.median():.2f}")
dom.round(3).to_csv("results/harmony_cluster_sample_purity.csv")

# ----------------------------------------------------------------- figures
sc.pl.umap(adata, color="cell_type_int", legend_loc="on data", legend_fontsize=7,
           title="Harmony-integrated cell types", show=False, save="_harmony_celltype.png")
sc.pl.umap(adata, color="sample_id", title="After Harmony: samples mixed",
           show=False, save="_harmony_sample.png")
sc.pl.umap(adata, color="group", show=False, save="_harmony_group.png")

ax = ct_frac.plot(kind="bar", figsize=(11, 5))
ax.set_ylabel("Fraction of cells")
ax.set_title("Integrated cell-type composition by group (GC vs GS)")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("figures/harmony_celltype_composition_by_group.png", dpi=120)

adata.write_h5ad("GSE273274_integrated.h5ad")
print("\nDone. Harmony figures in figures/, tables in results/.")
