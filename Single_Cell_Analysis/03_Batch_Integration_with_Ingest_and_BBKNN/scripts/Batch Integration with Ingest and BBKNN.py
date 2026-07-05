# ============================================================
# Batch Integration with Ingest and BBKNN
# Step 1. Environment setup and load data
# ============================================================

import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Set working directory
project_dir = os.path.expanduser("~/Desktop/Single cell python/Batch Integration")
os.chdir(project_dir)

# Create folder for figures
figures_dir = os.path.join(project_dir, "figures")
os.makedirs(figures_dir, exist_ok=True)

# Scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(
    dpi=120,
    frameon=False,
    figsize=(5, 4),
    facecolor="white"
)

# Check working directory
print("Current working directory:")
print(os.getcwd())

print("\nFiles in project folder:")
print(os.listdir())

# Load datasets
adata_ref = sc.read_h5ad("pbmc3k.h5ad")
adata_query = sc.read_h5ad("pbmc3k_1.h5ad")

print("\nReference dataset:")
print(adata_ref)

print("\nQuery dataset:")
print(adata_query)

# Keep shared genes only
shared_genes = adata_ref.var_names.intersection(adata_query.var_names)

adata_ref = adata_ref[:, shared_genes].copy()
adata_query = adata_query[:, shared_genes].copy()

print("\nNumber of shared genes:", len(shared_genes))