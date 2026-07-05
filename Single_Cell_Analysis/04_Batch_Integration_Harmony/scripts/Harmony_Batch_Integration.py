# ============================================================
# Harmony Batch Integration
# Step 1. Environment setup
# ============================================================

import os
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt

# Set working directory
project_dir = os.path.expanduser("~/Desktop/Single cell python/Batch Integration")
os.chdir(project_dir)

# Create figures folder
figures_dir = os.path.join(project_dir, "figures")
os.makedirs(figures_dir, exist_ok=True)

# Scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(
    dpi=120,
    frameon=False,
    figsize=(5,4),
    facecolor="white"
)

print(os.getcwd())
print(os.listdir())

# ============================================================
# Step 2. Read 10X datasets
# ============================================================

adata_bc2 = sc.read_10x_mtx("BC2/")
adata_bc10 = sc.read_10x_mtx("BC10/")
adata_bc21 = sc.read_10x_mtx("BC21/")

adata_bc2.obs["sample"] = "BC2"
adata_bc10.obs["sample"] = "BC10"
adata_bc21.obs["sample"] = "BC21"

print(adata_bc2)
print(adata_bc10)
print(adata_bc21)

# ============================================================
# Step 3. Merge three datasets
# ============================================================

adata = sc.AnnData.concatenate(
    adata_bc2,
    adata_bc10,
    adata_bc21,
    batch_key="batch",
    batch_categories=["BC2", "BC10", "BC21"]
)

print(adata)

print("\nBatch counts:")
print(adata.obs["batch"].value_counts())

print("\nSample counts:")
print(adata.obs["sample"].value_counts())

# ============================================================
# Step 4. Quality control
# ============================================================

# Filter low-quality cells and genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Identify mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True
)

# Plot QC metrics
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    show=False
)

plt.savefig(
    os.path.join(figures_dir, "05_harmony_qc_violin.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print(adata)

# ============================================================
# Step 5. Filtering, normalization, and highly variable genes
# ============================================================

# Filter cells based on QC metrics
adata = adata[adata.obs["n_genes_by_counts"] > 500, :].copy()
adata = adata[adata.obs["pct_counts_mt"] < 20, :].copy()

print("After QC filtering:")
print(adata)

# Normalize total counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# Log transformation
sc.pp.log1p(adata)

# Identify highly variable genes
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5
)

# Save raw data before subsetting HVGs
adata.raw = adata

# Keep highly variable genes only
adata = adata[:, adata.var["highly_variable"]].copy()

print("After selecting highly variable genes:")
print(adata)

# ============================================================
# Step 6. PCA and Harmony integration
# ============================================================

# Regress out technical effects
sc.pp.regress_out(
    adata,
    ["total_counts", "pct_counts_mt"]
)

# Scale data
sc.pp.scale(
    adata,
    max_value=10
)

# PCA
sc.pp.pca(adata)

# UMAP before Harmony
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.pl.umap(
    adata,
    color=["batch", "sample"],
    title=["Before Harmony: Batch", "Samples"],
    show=False
)

plt.savefig(
    os.path.join(figures_dir, "06_before_harmony.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# Harmony integration using harmonypy directly
import harmonypy as hm

harmony_out = hm.run_harmony(
    adata.obsm["X_pca"],
    adata.obs,
    "sample"
)

print("Z_corr shape:", harmony_out.Z_corr.shape)

adata.obsm["X_pca_harmony"] = harmony_out.Z_corr

print("Original PCA shape:", adata.obsm["X_pca"].shape)
print("Harmony PCA shape:", adata.obsm["X_pca_harmony"].shape)

# Rebuild neighbors using Harmony embedding
sc.pp.neighbors(
    adata,
    use_rep="X_pca_harmony"
)

# UMAP after Harmony
sc.tl.umap(
    adata,
    init_pos="X_pca_harmony"
)

# Leiden clustering
sc.tl.leiden(adata)

# Plot Harmony result
sc.pl.umap(
    adata,
    color=["sample", "leiden"],
    title=["After Harmony: Sample", "Leiden Clusters"],
    show=False
)

plt.savefig(
    os.path.join(figures_dir, "07_after_harmony.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Harmony integration completed.")






























