# =========================
# Step 1. Load packages and set paths
# =========================

import os
import scvelo as scv
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Project folder
project_dir = "/Users/yiyingzhang/Desktop/Single cell python/scvelo"
os.chdir(project_dir)

# Create output folders
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)

# scVelo settings
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params("scvelo")

print("Current working directory:", os.getcwd())
print("scvelo version:", scv.__version__)
print("scanpy version:", sc.__version__)
print("Step 1 finished.")

# =========================
# Step 2. Load loom data
# =========================

adata = sc.read("data/bc.loom", cache=False)
adata.var_names_make_unique()

print(adata)
print("Layers:", adata.layers.keys())
print("Step 2 finished.")

# =========================
# Step 3. Check spliced / unspliced proportions
# =========================

scv.pl.proportions(
    adata,
    show=False,
    save=None
)

plt.savefig(
    "figures/01_spliced_unspliced_proportions.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Step 3 finished. Figure saved.")

# =========================
# Step 4. Filter and normalize
# =========================

scv.pp.filter_and_normalize(
    adata,
    min_shared_counts=30,
    n_top_genes=2000
)

print(adata)
print("Step 4 finished.")

# =========================
# Step 5. Compute moments
# =========================

scv.pp.moments(
    adata,
    n_pcs=30,
    n_neighbors=30
)

print("Step 5 finished.")

# =========================
# Step 6. Compute RNA velocity
# =========================

scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

print("Step 6 finished.")

# =========================
# Step 7. Clustering and UMAP
# =========================

# Compute neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)

# Leiden clustering
sc.tl.leiden(adata)

# UMAP
sc.tl.umap(adata)

# Plot UMAP
sc.pl.umap(
    adata,
    color="leiden",
    show=False
)

plt.savefig(
    "figures/02_umap_leiden.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Step 7 finished.")

# =========================
# Step 8. Velocity stream plot
# =========================

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="leiden",
    show=False
)

plt.savefig(
    "figures/03_velocity_stream_leiden.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Step 8 finished.")

# =========================
# Step 9. Velocity grid plot
# =========================

scv.pl.velocity_embedding_grid(
    adata,
    basis="umap",
    color="leiden",
    show=False
)

plt.savefig(
    "figures/04_velocity_grid_leiden.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Step 9 finished.")

# =========================
# Step 10. Velocity embedding
# =========================

scv.pl.velocity_embedding(
    adata,
    basis="umap",
    arrow_length=3,
    arrow_size=2,
    dpi=120,
    color="leiden",
    show=False
)

plt.savefig(
    "figures/05_velocity_embedding_leiden.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Step 10 finished.")

# =========================
# Step 11. Single-gene velocity
# =========================

gene = "DOCK4"

if gene in adata.var_names:
    scv.pl.velocity(
        adata,
        [gene],
        dpi=120,
        color="leiden",
        show=False
    )

    plt.savefig(
        f"figures/06_velocity_{gene}.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Step 11 finished. {gene} figure saved.")
else:
    print(f"{gene} is not found in this dataset.")

# =========================
# Step 12. Rank velocity genes
# =========================

scv.tl.rank_velocity_genes(
    adata,
    groupby="leiden",
    min_corr=0.3
)

df_velocity_genes = pd.DataFrame(
    adata.uns["rank_velocity_genes"]["names"]
)

print(df_velocity_genes.head())

df_velocity_genes.to_csv(
    "results/rank_velocity_genes_by_leiden.csv",
    index=False
)

print("Step 12 finished. Results saved.")

# =========================
# Step 13. Cell cycle scores
# =========================

scv.tl.score_genes_cell_cycle(adata)

scv.pl.scatter(
    adata,
    color_gradients=["S_score", "G2M_score"],
    smooth=True,
    perc=[5, 95],
    show=False
)

plt.savefig(
    "figures/07_cell_cycle_scores.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Step 13 finished.")

# =========================
# Step 14. Velocity confidence
# =========================

scv.tl.velocity_confidence(adata)

scv.pl.scatter(
    adata,
    color=["velocity_length", "velocity_confidence"],
    cmap="coolwarm",
    perc=[5, 95],
    show=False
)

plt.savefig(
    "figures/08_velocity_confidence.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Step 14 finished.")

# =========================
# Step 15. Velocity pseudotime
# =========================

scv.tl.velocity_pseudotime(adata)

scv.pl.scatter(
    adata,
    color="velocity_pseudotime",
    cmap="gnuplot",
    show=False
)

plt.savefig(
    "figures/09_velocity_pseudotime.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Step 15 finished.")

# =========================
# Step 16. PAGA
# =========================

sc.tl.paga(adata, groups="leiden")

scv.pl.paga(
    adata,
    basis="umap",
    size=50,
    alpha=0.1,
    min_edge_width=2,
    node_size_scale=1.5,
    show=False
)

plt.savefig(
    "figures/10_paga.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Step 16 finished.")

# =========================
# Step 17. Dynamical model
# =========================

scv.tl.recover_dynamics(adata)

scv.tl.velocity(
    adata,
    mode="dynamical"
)

scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(
    adata,
    basis="umap",
    color="leiden",
    show=False
)

plt.savefig(
    "figures/11_dynamical_velocity_stream.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Step 17 finished.")

# =========================
# Step 18. Latent time
# =========================

scv.tl.latent_time(adata)

scv.pl.scatter(
    adata,
    color="latent_time",
    cmap="gnuplot",
    show=False
)

plt.savefig(
    "figures/12_latent_time.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Step 18 finished.")

# =========================
# Step 19. Dynamical genes heatmap
# =========================

top_genes = adata.var["fit_likelihood"].sort_values(
    ascending=False
).index[:300]

scv.pl.heatmap(
    adata,
    var_names=top_genes,
    tkey="latent_time",
    n_convolve=100,
    col_color="leiden",
    show=False
)

plt.savefig(
    "figures/13_dynamical_genes_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Step 19 finished.")

# =========================
# Step 20. Top-likelihood genes scatter
# =========================

top15 = adata.var["fit_likelihood"].sort_values(
    ascending=False
).index[:15]

scv.pl.scatter(
    adata,
    basis=top15,
    ncols=5,
    frameon=False,
    color="leiden",
    show=False
)

plt.savefig(
    "figures/14_top_likelihood_genes.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Step 20 finished.")

# =========================
# Step 21. Selected genes
# =========================

genes = [
    g for g in ["DOCK4", "COL1A2", "IL32"]
    if g in adata.var_names
]

print("Genes found:", genes)

if len(genes) > 0:

    # Expression on UMAP
    scv.pl.scatter(
        adata,
        color=genes,
        frameon=False,
        show=False
    )

    plt.savefig(
        "figures/15_selected_genes_expression.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    # Expression along latent time
    scv.pl.scatter(
        adata,
        x="latent_time",
        y=genes,
        frameon=False,
        color="leiden",
        show=False
    )

    plt.savefig(
        "figures/16_selected_genes_latent_time.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

print("Step 21 finished.")

# =========================
# Step 22. Rank dynamical genes
# =========================

scv.tl.rank_dynamical_genes(
    adata,
    groupby="leiden"
)

df_dyn = scv.get_df(
    adata,
    "rank_dynamical_genes/names"
)

print(df_dyn.head(5))

df_dyn.to_csv(
    "results/rank_dynamical_genes_by_leiden.csv",
    index=False
)

print("Step 22 finished. Dynamical gene results saved.")





































