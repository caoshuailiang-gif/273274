"""
GSE273274 cell-type annotation (marker-based, brain/glioma context).

Loads the processed object, scores canonical marker gene sets per cell with
sc.tl.score_genes, assigns each Leiden cluster the cell type with the highest
mean score, and writes annotated UMAPs, a dotplot, and a cluster->celltype map.
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

# Canonical markers for adult brain / glioma microenvironment
marker_sets = {
    "Oligodendrocyte":      ["PLP1", "MBP", "MOG", "MOBP", "MAG", "CLDN11", "CNP"],
    "OPC":                  ["PDGFRA", "CSPG4", "OLIG1", "OLIG2", "SOX10", "BCAN"],
    "Astrocyte":            ["GFAP", "AQP4", "SLC1A3", "ALDH1L1", "S100B", "GJA1"],
    "Neuron":               ["RBFOX3", "SNAP25", "SYT1", "NEFL", "GAP43", "STMN2", "MAP2"],
    "Microglia_Myeloid":    ["PTPRC", "AIF1", "CSF1R", "P2RY12", "CX3CR1", "C1QA",
                             "C1QB", "C1QC", "CD68", "TYROBP"],
    "T_NK_cell":            ["CD3D", "CD3E", "CD2", "CD8A", "IL7R", "NKG7", "GNLY"],
    "B_Plasma":             ["CD79A", "CD79B", "MS4A1", "IGHG1", "MZB1"],
    "Endothelial":          ["CLDN5", "PECAM1", "VWF", "FLT1", "A2M", "ITM2A"],
    "Pericyte_Mural":       ["RGS5", "PDGFRB", "ACTA2", "NOTCH3", "MYH11"],
    "Proliferating":        ["MKI67", "TOP2A", "CENPF", "CCNB1", "UBE2C"],
    "Tumor_glioma":         ["EGFR", "SOX2", "SOX4", "OLIG2", "CD24", "PTPRZ1",
                             "BCAN", "CHI3L1", "VIM"],
}

# ---------------------------------------------------------- load (full genes)
adata = sc.read_h5ad("GSE273274_processed.h5ad")
full = adata.raw.to_adata()          # log-normalized, all genes
full.obs = adata.obs.copy()          # carry leiden / group / sample_id
full.obsm["X_umap"] = adata.obsm["X_umap"]
print(f"Loaded {full.n_obs} cells x {full.n_vars} genes; {full.obs['leiden'].nunique()} clusters")

# ----------------------------------------------------------------- scoring
score_cols = []
for ct, genes in marker_sets.items():
    present = [g for g in genes if g in full.var_names]
    sc.tl.score_genes(full, present, score_name=f"score_{ct}", ctrl_size=50)
    score_cols.append(f"score_{ct}")
    if len(present) < len(genes):
        print(f"  {ct}: missing {set(genes) - set(present)}")

# mean score per cluster -> assign argmax
mean_scores = full.obs.groupby("leiden", observed=True)[score_cols].mean()
mean_scores.columns = [c.replace("score_", "") for c in mean_scores.columns]
cluster_to_ct = mean_scores.idxmax(axis=1)
mean_scores["assigned"] = cluster_to_ct
mean_scores.to_csv("results/cluster_celltype_scores.csv")
print("\nCluster -> cell type:")
print(cluster_to_ct.to_string())

# ---------------------------------------------- attach labels & save mapping
full.obs["cell_type"] = full.obs["leiden"].map(cluster_to_ct).astype("category")
adata.obs["cell_type"] = full.obs["cell_type"].values

mapping = pd.DataFrame({
    "leiden": cluster_to_ct.index,
    "cell_type": cluster_to_ct.values,
    "n_cells": full.obs["leiden"].value_counts().reindex(cluster_to_ct.index).values,
})
mapping.to_csv("results/cluster_to_celltype.csv", index=False)

# ---------------------------------------- composition: cell type x group
ct_counts = pd.crosstab(full.obs["cell_type"], full.obs["group"])
ct_frac = pd.crosstab(full.obs["cell_type"], full.obs["group"], normalize="columns")
ct_counts.to_csv("results/celltype_counts_by_group.csv")
ct_frac.to_csv("results/celltype_fraction_by_group.csv")
print("\nCell type fraction by group:\n", ct_frac.round(3).to_string())

# ----------------------------------------------------------------- figures
sc.pl.umap(full, color="cell_type", legend_loc="on data", legend_fontsize=7,
           title="Cell type annotation", show=False, save="_celltype.png")
sc.pl.umap(full, color="cell_type", groups=None, show=False, save="_celltype_legend.png")

# dotplot of representative markers per assigned type
dot_markers = {ct: [g for g in gs if g in full.var_names][:4]
               for ct, gs in marker_sets.items()}
sc.pl.dotplot(full, dot_markers, groupby="cell_type", standard_scale="var",
              show=False, save="_celltype_markers.png")

ax = ct_frac.plot(kind="bar", figsize=(11, 5))
ax.set_ylabel("Fraction of cells")
ax.set_title("Cell-type composition by group (GC vs GS)")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("figures/celltype_composition_by_group.png", dpi=120)

# ----------------------------------------------------------------- save
adata.write_h5ad("GSE273274_processed.h5ad")
print("\nDone. Annotated UMAP/dotplot in figures/, tables in results/.")
