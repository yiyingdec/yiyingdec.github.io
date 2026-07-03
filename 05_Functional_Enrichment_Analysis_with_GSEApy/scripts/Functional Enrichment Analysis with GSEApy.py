# ==========================
# P12 Functional Enrichment Analysis with GSEApy
# Step 1. Import packages
# ==========================

import os
import time

import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp

import matplotlib.pyplot as plt
import networkx as nx

print("Scanpy version:", sc.__version__)
print("GSEApy version:", gp.__version__)

# ==========================
# Step 2. Set working directory
# ==========================

os.chdir("/Users/yiyingzhang/Desktop/Single cell python/GSEA/GSEApy-master")

print("Current working directory:")
print(os.getcwd())

# ==========================
# Step 3. Read the AnnData object
# ==========================

adata = sc.read_h5ad("tests/data/ifnb.h5ad")

print(adata)

# ==========================
# Step 4. View metadata
# ==========================

print(adata.obs.head())

# ==========================
# Step 5. Save raw counts and preprocess
# ==========================

adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print(adata)

# ==========================
# Step 6. Set group order and sort cells
# ==========================

adata.obs["stim"] = pd.Categorical(
    adata.obs["stim"],
    categories=["STIM", "CTRL"],
    ordered=True
)

indices = adata.obs.sort_values(["seurat_annotations", "stim"]).index
adata = adata[indices, :].copy()

print(adata.obs[["stim", "seurat_annotations"]].head())

# ==========================
# Step 7. Subset CD14 Mono cells
# ==========================

bdata = adata[adata.obs["seurat_annotations"] == "CD14 Mono"].copy()

print(bdata)
print(bdata.obs["stim"].value_counts())

# ==========================
# Step 8. Run GSEA
# ==========================

t1 = time.time()

res = gp.gsea(
    data=bdata.to_df().T,
    gene_sets="GO_Biological_Process_2021",
    cls=bdata.obs["stim"],
    permutation_num=1000,
    permutation_type="phenotype",
    outdir=None,
    method="s2n",
    threads=16
)

t2 = time.time()

print(f"Running time: {t2 - t1:.2f} seconds")

# ==========================
# Step 9. View enrichment results
# ==========================

print(res.res2d.head(10))

print(res.ranking.shape)

# ==========================
# Step 10. Heatmap of leading-edge genes
# ==========================

i = 7

genes = res.res2d.Lead_genes.iloc[i].split(";")

ax = gp.heatmap(
    df=res.heatmat.loc[genes],
    z_score=None,
    title=res.res2d.Term.iloc[i],
    figsize=(7, 6),
    cmap=plt.cm.viridis,
    xticklabels=False
)

plt.savefig("../figures/heatmap_leading_edge.png", dpi=300, bbox_inches="tight")
plt.show()

# ==========================
# Step 11. GSEA enrichment plot
# ==========================

term = res.res2d.Term

axs = res.plot(terms=term[:5])

plt.savefig(
    "../figures/gsea_enrichment_top5.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# ==========================
# Step 12. Differential expression analysis
# ==========================

sc.tl.rank_genes_groups(
    bdata,
    groupby="stim",
    use_raw=False,
    method="wilcoxon",
    groups=["STIM"],
    reference="CTRL"
)

sc.pl.rank_genes_groups(
    bdata,
    n_genes=25,
    sharey=False,
    show=False
)

plt.savefig(
    "../figures/deg_top25_stim_vs_ctrl.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# ==========================
# Step 13. Extract DE results
# ==========================

deg = sc.get.rank_genes_groups_df(
    bdata,
    group="STIM"
)

print(deg.head())

# ==========================
# Step 14. Create preranked gene list
# ==========================

ranking = deg[["names", "scores"]].dropna()

print(ranking.head())
print(ranking.shape)

# ==========================
# Step 15. Preranked GSEA
# ==========================

pre_res = gp.prerank(
    rnk=ranking,
    gene_sets="KEGG_2016",
    outdir=None
)

print(pre_res.res2d.head(5))

# ==========================
# Step 16. Plot preranked GSEA results
# ==========================

term2 = pre_res.res2d.Term

axes = pre_res.plot(terms=term2[:5])

plt.savefig(
    "../figures/prerank_gsea_kegg_top5.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# ==========================
# Step 17. Select significant DEGs
# ==========================

sig_genes = deg.loc[
    (deg["pvals_adj"] < 0.05) &
    (deg["logfoldchanges"] > 1),
    "names"
].tolist()

print(f"Number of significant genes: {len(sig_genes)}")
print(sig_genes[:10])

# ==========================
# Step 18. Enrichr analysis
# ==========================

enr = gp.enrichr(
    gene_list=sig_genes,
    gene_sets="GO_Biological_Process_2021",
    organism="homo sapiens",
    outdir=None
)

print(enr.results.head())

# ==========================
# Step 19. Enrichr Dotplot
# ==========================

ax = gp.dotplot(
    enr.results,
    column="Adjusted P-value",
    x="Gene_set",
    size=8,
    top_term=10,
    figsize=(8, 6)
)

plt.savefig(
    "../figures/enrichr_dotplot.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# ==========================
# Step 20. Build enrichment map
# ==========================

nodes, edges = gp.enrichment_map(enr.results)

print(nodes.head())
print(edges.head())

# ==========================
# Step 21. Draw enrichment map
# ==========================

G = nx.from_pandas_edgelist(
    edges,
    source="src_idx",
    target="targ_idx"
)

fig, ax = plt.subplots(figsize=(10, 8))

pos = nx.spring_layout(G, seed=42)

nx.draw_networkx_nodes(
    G,
    pos,
    node_size=300,
    node_color="skyblue",
    ax=ax
)

nx.draw_networkx_edges(
    G,
    pos,
    alpha=0.5,
    ax=ax
)
labels = {i: nodes.loc[i, "Term"] for i in nodes.index}

nx.draw_networkx_labels(
    G,
    pos,
    labels=labels,
    font_size=6
)
plt.axis("off")

plt.savefig(
    "../figures/enrichment_map.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# ==========================
# Step 22. Check saved figures
# ==========================

import os

print(os.listdir("../figures"))































