"""
GSE273274 single-cell RNA-seq analysis pipeline.

Two conditions, three replicates each:
  GC : 87_R1_GC, 89_R2.GC, 91_R3.GC
  GS : 86_R1.GS, 88_R2.GS, 90_R3.GS

Pipeline: load -> QC -> filter -> normalize -> HVG -> scale -> PCA ->
neighbors -> Leiden clustering -> UMAP -> marker genes -> composition.
Outputs figures to figures/ and a processed .h5ad (gitignored).
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

DATA_DIR = "273274"
samples = {
    "87_R1_GC": {"path": os.path.join(DATA_DIR, "GC", "87_R1_GC"), "group": "GC", "rep": "R1"},
    "89_R2_GC": {"path": os.path.join(DATA_DIR, "GC", "89_R2.GC"), "group": "GC", "rep": "R2"},
    "91_R3_GC": {"path": os.path.join(DATA_DIR, "GC", "91_R3.GC"), "group": "GC", "rep": "R3"},
    "86_R1_GS": {"path": os.path.join(DATA_DIR, "GS", "86_R1.GS"), "group": "GS", "rep": "R1"},
    "88_R2_GS": {"path": os.path.join(DATA_DIR, "GS", "88_R2.GS"), "group": "GS", "rep": "R2"},
    "90_R3_GS": {"path": os.path.join(DATA_DIR, "GS", "90_R3.GS"), "group": "GS", "rep": "R3"},
}

# ---------------------------------------------------------------- load + merge
adatas = []
for sid, info in samples.items():
    print(f"Reading {sid} ...")
    a = sc.read_10x_mtx(info["path"], var_names="gene_symbols", cache=True)
    a.obs_names = [f"{sid}_{bc}" for bc in a.obs_names]
    a.obs["sample_id"] = sid
    a.obs["group"] = info["group"]
    a.obs["replicate"] = info["rep"]
    print(f"  {a.n_obs} cells x {a.n_vars} genes")
    adatas.append(a)

adata = sc.concat(adatas, join="outer", fill_value=0)
adata.var_names_make_unique()
adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
adata.obs["group"] = adata.obs["group"].astype("category")
print(f"\nCombined raw: {adata.n_obs} cells x {adata.n_vars} genes")

# ----------------------------------------------------------------------- QC
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, percent_top=None)

qc_summary = adata.obs.groupby("sample_id")[
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"]
].median()
qc_summary["n_cells"] = adata.obs["sample_id"].value_counts()
qc_summary.to_csv("results/qc_summary_per_sample.csv")
print("\nPer-sample QC medians:\n", qc_summary)

sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
             groupby="sample_id", rotation=45, show=False, save="_qc_pre.png")

# --------------------------------------------------------------- filtering
n0 = adata.n_obs
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[(adata.obs.n_genes_by_counts < 7000) & (adata.obs.pct_counts_mt < 20)].copy()
print(f"\nFiltered: {n0} -> {adata.n_obs} cells ({n0 - adata.n_obs} removed)")

# store raw counts
adata.layers["counts"] = adata.X.copy()

# ----------------------------------------------------------- normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# --------------------------------------------------- HVG / scale / PCA / UMAP
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample_id")
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(adata)
print(f"\nLeiden clusters: {adata.obs['leiden'].nunique()}")

# ---------------------------------------------------------------- figures
sc.pl.umap(adata, color=["leiden"], legend_loc="on data", show=False, save="_leiden.png")
sc.pl.umap(adata, color=["group"], show=False, save="_group.png")
sc.pl.umap(adata, color=["sample_id"], show=False, save="_sample.png")
sc.pl.umap(adata, color=["pct_counts_mt", "n_genes_by_counts"], show=False, save="_qc.png")

# ------------------------------------------------------------ marker genes
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=15, sharey=False, show=False, save="_markers.png")
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv("results/cluster_markers.csv", index=False)
top = (markers.sort_values(["group", "scores"], ascending=[True, False])
       .groupby("group").head(10))
top.to_csv("results/top10_markers_per_cluster.csv", index=False)

# ------------------------------------------------- composition: cluster x group
comp = pd.crosstab(adata.obs["leiden"], adata.obs["group"])
comp_frac = pd.crosstab(adata.obs["leiden"], adata.obs["group"], normalize="columns")
comp.to_csv("results/cluster_counts_by_group.csv")
comp_frac.to_csv("results/cluster_fraction_by_group.csv")

ax = comp_frac.plot(kind="bar", stacked=False, figsize=(10, 5))
ax.set_ylabel("Fraction of cells")
ax.set_xlabel("Leiden cluster")
ax.set_title("Cluster composition by group (GC vs GS)")
plt.tight_layout()
plt.savefig("figures/cluster_composition_by_group.png", dpi=120)

# ----------------------------------------------- GC vs GS pseudobulk-ish DE
sc.tl.rank_genes_groups(adata, "group", groups=["GS"], reference="GC",
                        method="wilcoxon", key_added="de_GS_vs_GC")
de = sc.get.rank_genes_groups_df(adata, group="GS", key="de_GS_vs_GC")
de.to_csv("results/DE_GS_vs_GC.csv", index=False)
print("\nTop GS-vs-GC genes:\n", de.head(15).to_string(index=False))

# ------------------------------------------------------------------- save
adata.write_h5ad("GSE273274_processed.h5ad")
print("\nDone. Figures in figures/, tables in results/.")
