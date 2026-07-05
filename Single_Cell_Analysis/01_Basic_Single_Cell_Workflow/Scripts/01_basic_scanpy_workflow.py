"""
Basic scRNA-seq Analysis with Scanpy
Author: Yiying Zhang
"""

##################################################
# Step 0. Import packages and set output folders
##################################################

import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

sc.settings.verbosity = 3
sc.set_figure_params(dpi=80, facecolor="white")

os.makedirs("figures", exist_ok=True)
sc.settings.figdir = "figures"

sc.logging.print_header()
print(os.getcwd())


##################################################
# Step 1. Load 10x Genomics scRNA-seq Data
##################################################

adata = sc.read_10x_mtx(
    "/Users/yiyingzhang/Desktop/Single cell python/basic flow/GSE294482_RAW",
    var_names="gene_symbols",
    cache=True,
)

adata.var_names_make_unique()
print(adata)


##################################################
# Step 2. Initial Quality Control
##################################################

sc.pl.highest_expr_genes(
    adata,
    n_top=20,
    save="_highest_expr_genes.png"
)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(adata)


##################################################
# Step 3. Mitochondrial QC Metrics
##################################################

adata.var["mt"] = adata.var_names.str.startswith(("mt-", "Mt-", "MT-"))

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True,
)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save="_qc_violin.png"
)

sc.pl.scatter(
    adata,
    x="total_counts",
    y="pct_counts_mt",
    save="_counts_vs_mt.png"
)

sc.pl.scatter(
    adata,
    x="total_counts",
    y="n_genes_by_counts",
    save="_counts_vs_genes.png"
)


##################################################
# Step 4. Filter Low-quality Cells
##################################################

adata = adata[
    (adata.obs.n_genes_by_counts < 2500)
    & (adata.obs.n_genes_by_counts > 200)
    & (adata.obs.pct_counts_mt < 5),
    :,
].copy()

print(adata)


##################################################
# Step 5. Normalize and Log-transform Data
##################################################

adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print(adata)


##################################################
# Step 6. Highly Variable Genes
##################################################

sc.pp.highly_variable_genes(
    adata,
    layer="counts",
    n_top_genes=2000,
    flavor="seurat_v3",
)

sc.pl.highly_variable_genes(
    adata,
    save="_highly_variable_genes.png"
)

print(adata.var[adata.var["highly_variable"]].head(20))


##################################################
# Step 7. Scale Data and Run PCA
##################################################

adata.layers["scaled"] = adata.X.copy()

sc.pp.regress_out(
    adata,
    ["total_counts", "pct_counts_mt"],
    layer="scaled",
)

sc.pp.scale(
    adata,
    max_value=10,
    layer="scaled",
)

sc.pp.pca(
    adata,
    layer="scaled",
    svd_solver="arpack",
)

sc.pl.pca(
    adata,
    color="total_counts",
    save="_pca_total_counts.png"
)

sc.pl.pca_variance_ratio(
    adata,
    n_pcs=20,
    log=True,
    save="_pca_variance.png"
)


##################################################
# Step 8. Neighbors, Leiden Clustering, and UMAP
##################################################

sc.pp.neighbors(
    adata,
    n_neighbors=10,
    n_pcs=40,
)

sc.tl.leiden(
    adata,
    resolution=0.7,
    random_state=0,
    flavor="igraph",
    n_iterations=2,
    directed=False,
)

sc.tl.umap(adata)

sc.pl.umap(
    adata,
    color="leiden",
    save="_umap_leiden.png"
)


##################################################
# Step 9. Visualize Selected Marker Genes
##################################################

sc.pl.umap(
    adata,
    color=["leiden", "Cd14", "Nkg7"],
    save="_umap_cd14_nkg7.png"
)

sc.pl.umap(
    adata,
    color=["Npas2", "Aff3", "Prim2"],
    save="_umap_selected_genes.png"
)


##################################################
# Step 10. Identify Marker Genes
##################################################

sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon",
    mask_var="highly_variable",
)

sc.pl.rank_genes_groups(
    adata,
    n_genes=25,
    sharey=False,
    save="_rank_genes_groups.png"
)

marker_genes = ["Cd14", "Nkg7", "Npas2", "Aff3", "Prim2"]

sc.pl.dotplot(
    adata,
    marker_genes,
    groupby="leiden",
    save="_marker_dotplot.png"
)

sc.pl.stacked_violin(
    adata,
    marker_genes,
    groupby="leiden",
    save="_marker_stacked_violin.png"
)


##################################################
# Step 11. Save Processed Data
##################################################

adata.write("basic_scanpy_workflow.h5ad")

print("Analysis finished successfully.")





































