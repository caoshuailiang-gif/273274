"""
Re-annotation: label the neural-progenitor / radial-glia / astrocyte-like
malignant population as Tumor_glioma (this is a glioma dataset; those neural
lineage states are the tumor compartment). Normal microenvironment types are
kept as-is. Regenerates the integrated cell-type figures and tables.
"""

import os
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sc.settings.verbosity = 1
sc.settings.figdir = "figures"
sc.settings.set_figure_params(dpi=120, frameon=False)

# neural malignant states -> Tumor_glioma; everything else stays normal TME
TUMOR_FROM = {"NeuralProgenitor_RadialGlia", "Astrocyte"}

adata = sc.read_h5ad("GSE273274_integrated.h5ad")
adata.obs["cell_type"] = adata.obs["cell_type_int"].astype(str).apply(
    lambda x: "Tumor_glioma" if x in TUMOR_FROM else x).astype("category")
print(adata.obs["cell_type"].value_counts())

# cluster -> final label map
cl_map = (adata.obs.groupby("leiden", observed=True)["cell_type"]
          .agg(lambda s: s.value_counts().index[0]))
pd.DataFrame({"leiden": cl_map.index, "cell_type": cl_map.values,
              "n_cells": adata.obs["leiden"].value_counts().reindex(cl_map.index).values
              }).to_csv("results/harmony_cluster_to_celltype.csv", index=False)

# composition by group
counts = pd.crosstab(adata.obs["cell_type"], adata.obs["group"])
frac = pd.crosstab(adata.obs["cell_type"], adata.obs["group"], normalize="columns")
counts.to_csv("results/harmony_celltype_counts_by_group.csv")
frac.to_csv("results/harmony_celltype_fraction_by_group.csv")
print("\nCell-type fraction by group:\n", frac.round(3).to_string())

# figures (overwrite integrated cell-type views)
sc.pl.umap(adata, color="cell_type", legend_loc="on data", legend_fontsize=7,
           title="Harmony-integrated cell types", show=False, save="_harmony_celltype.png")

ax = frac.plot(kind="bar", figsize=(10, 5))
ax.set_ylabel("Fraction of cells")
ax.set_title("Integrated cell-type composition by group (GC vs GS)")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("figures/harmony_celltype_composition_by_group.png", dpi=120)

adata.write_h5ad("GSE273274_integrated.h5ad")
print("\nDone. Re-annotated as Tumor_glioma; figures/tables updated.")
