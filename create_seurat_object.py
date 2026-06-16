import scanpy as sc
import anndata as ad
import pandas as pd
import os

sc.settings.verbosity = 2

DATA_DIR = "D:/cc-deepseek/GSE273274/273274"

samples = {
    "87_R1_GC": {"path": os.path.join(DATA_DIR, "GC", "87_R1_GC"), "group": "GC", "sample": "R1"},
    "89_R2_GC": {"path": os.path.join(DATA_DIR, "GC", "89_R2.GC"), "group": "GC", "sample": "R2"},
    "91_R3_GC": {"path": os.path.join(DATA_DIR, "GC", "91_R3.GC"), "group": "GC", "sample": "R3"},
    "86_R1_GS": {"path": os.path.join(DATA_DIR, "GS", "86_R1.GS"), "group": "GS", "sample": "R1"},
    "88_R2_GS": {"path": os.path.join(DATA_DIR, "GS", "88_R2.GS"), "group": "GS", "sample": "R2"},
    "90_R3_GS": {"path": os.path.join(DATA_DIR, "GS", "90_R3.GS"), "group": "GS", "sample": "R3"},
}

adatas = []
for sample_id, info in samples.items():
    print(f"Reading {sample_id} ...")
    adata = sc.read_10x_mtx(
        info["path"],
        var_names="gene_symbols",
        cache=False,
    )
    adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]
    adata.obs["sample_id"] = sample_id
    adata.obs["group"] = info["group"]
    adata.obs["replicate"] = info["sample"]
    print(f"  {adata.n_obs} cells x {adata.n_vars} genes")
    adatas.append(adata)

print("\nMerging all samples ...")
adata_combined = ad.concat(adatas, join="outer", fill_value=0)
adata_combined.var_names_make_unique()

print(f"\nCombined object: {adata_combined.n_obs} cells x {adata_combined.n_vars} genes")
print("\nSample breakdown:")
print(adata_combined.obs["sample_id"].value_counts())
print("\nGroup breakdown:")
print(adata_combined.obs["group"].value_counts())

out_path = "D:/cc-deepseek/GSE273274/GSE273274_combined_raw.h5ad"
adata_combined.write_h5ad(out_path)
print(f"\nSaved to: {out_path}")
