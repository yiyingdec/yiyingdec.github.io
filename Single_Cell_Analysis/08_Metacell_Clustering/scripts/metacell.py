import os
os.getcwd()

# ============================================
# Step 1. Import packages
# ============================================

import os
import shutil
from math import hypot
from typing import *

import anndata as ad
import matplotlib.pyplot as plt
import metacells as mc
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sb

# ============================================
# Step 2. Configure the environment
# ============================================

# Use SVG for high-quality figures
%config InlineBackend.figure_formats = ["svg"]

# Set plotting style
sb.set_style("white")

# Prevent inefficient matrix operations
mc.ut.allow_inefficient_layout(False)

# ============================================
# Step 3. Load a smaller dataset
# ============================================

full = ad.read_h5ad("blobs/blood_aging.clean.h5ad")

mc.ut.top_level(full)
mc.ut.set_name(full, "blood_aging_clean")

print(f"Full: {full.n_obs} cells, {full.n_vars} genes")

# ============================================
# Step 4. Inspect the dataset
# ============================================

print(full)

print(full.obs.head())
print(full.var.head())

print(full.obs.columns)
print(full.var.columns)

# ============================================
# Step 5. Dataset summary
# ============================================

print(full)
print(f"Number of cells : {full.n_obs}")
print(f"Number of genes : {full.n_vars}")

# ============================================
# Step 6. Cleaning the data
# ============================================

PROPERLY_SAMPLED_MIN_CELL_TOTAL = 800
PROPERLY_SAMPLED_MAX_CELL_TOTAL = 20000

# ============================================
# Step 7. Plot UMI distribution
# ============================================

cell_total = np.asarray(full.X.sum(axis=1)).ravel()

print(f"Min UMIs: {cell_total.min()}")
print(f"Max UMIs: {cell_total.max()}")
print(f"Median UMIs: {np.median(cell_total)}")

low_quality = cell_total < PROPERLY_SAMPLED_MIN_CELL_TOTAL
high_quality = cell_total > PROPERLY_SAMPLED_MAX_CELL_TOTAL

print(f"Will exclude {low_quality.sum()} cells with less than {PROPERLY_SAMPLED_MIN_CELL_TOTAL} UMIs")
print(f"Will exclude {high_quality.sum()} cells with more than {PROPERLY_SAMPLED_MAX_CELL_TOTAL} UMIs")

plt.figure(figsize=(6, 4))
plt.hist(cell_total, bins=100)
plt.axvline(PROPERLY_SAMPLED_MIN_CELL_TOTAL, linestyle="--")
plt.axvline(PROPERLY_SAMPLED_MAX_CELL_TOTAL, linestyle="--")
plt.xscale("log")
plt.xlabel("UMIs per cell")
plt.ylabel("Number of cells")
plt.title("UMI distribution before filtering")
plt.show()

# ============================================
# Step 8. Filter cells
# ============================================

keep_cells = (
    (cell_total >= PROPERLY_SAMPLED_MIN_CELL_TOTAL)
    & (cell_total <= PROPERLY_SAMPLED_MAX_CELL_TOTAL)
)

filtered = full[keep_cells].copy()

print(filtered)
print(f"Remaining cells: {filtered.n_obs}")

# ============================================
# Step 9. Save filtered dataset
# ============================================

import os

os.makedirs("results", exist_ok=True)
os.makedirs("figures", exist_ok=True)

filtered.write_h5ad("results/blood_aging_filtered.h5ad")

print("Saved filtered dataset:")
print("results/blood_aging_filtered.h5ad")

# ============================================
# Step 10. Available API
# ============================================

[x for x in dir(mc.pl) if not x.startswith("_")]

# ============================================
# Step 11. Compute direct metacells
# ============================================

# Add missing gene annotations required by the new Metacells API
if "lateral_gene" not in filtered.var.columns:
    filtered.var["lateral_gene"] = False

if "ignored_gene" not in filtered.var.columns:
    filtered.var["ignored_gene"] = False

if "essential_gene" not in filtered.var.columns:
    filtered.var["essential_gene"] = False

# Run direct metacell construction
mc.pl.compute_direct_metacells(
    filtered,
    random_seed=1234
)

print("Metacell assignment finished.")
print(filtered.obs["metacell"].value_counts().head())
print(f"Number of metacells: {filtered.obs['metacell'].nunique()}")

# ============================================
# Step 12. Inspect metacell assignment
# ============================================

print(filtered.obs.head())

print(filtered.obs["metacell"].value_counts().head(10))

print(f"Total metacells: {filtered.obs['metacell'].nunique()}")

print(f"Assigned cells: {(filtered.obs['metacell'] >= 0).sum()}")

import os

os.makedirs("figures", exist_ok=True)

plt.savefig(
    "figures/01_UMI_distribution_before_filtering.png",
    dpi=300,
    bbox_inches="tight"
)

print("Saved: figures/01_UMI_distribution_before_filtering.png")

# ============================================
# Step 13. Save metacell assignment
# ============================================

import os

os.makedirs("results", exist_ok=True)

filtered.obs.to_csv("results/metacell_assignment.csv")

filtered.write_h5ad("results/blood_aging_metacell.h5ad")

print("Saved:")
print("results/metacell_assignment.csv")
print("results/blood_aging_metacell.h5ad")

# ============================================
# Step 14. Plot metacell size distribution
# ============================================

metacell_sizes = filtered.obs["metacell"].value_counts()
metacell_sizes = metacell_sizes[metacell_sizes.index != -1]

print(metacell_sizes.describe())

plt.figure(figsize=(6, 4))
plt.hist(metacell_sizes, bins=50)
plt.xlabel("Cells per metacell")
plt.ylabel("Number of metacells")
plt.title("Metacell size distribution")

plt.savefig(
    "figures/02_metacell_size_distribution.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Saved: figures/02_metacell_size_distribution.png")

# ============================================
# Step 15. Inspect metacell labels
# ============================================

print(filtered.obs[["metacell"]].head(20))

print(filtered.obs["metacell"].nunique())

print(filtered.obs["metacell"].value_counts().head(20))

# ============================================
# Step 16. Save final object
# ============================================

filtered.write_h5ad("results/final_metacell_results.h5ad")

print("Saved:")
print("results/final_metacell_results.h5ad")

# ==========================================
# Step 17. Plot MetaCell size distribution
# ==========================================

import matplotlib.pyplot as plt
import os

os.makedirs("figures", exist_ok=True)

metacell_sizes = filtered.obs["metacell"].value_counts()

# 去掉未分配的 -1
metacell_sizes = metacell_sizes[metacell_sizes.index != -1]

print(metacell_sizes.describe())

plt.figure(figsize=(6,4))
plt.hist(metacell_sizes, bins=40)
plt.xlabel("Cells per MetaCell")
plt.ylabel("Number of MetaCells")
plt.title("MetaCell size distribution")

plt.tight_layout()

plt.savefig(
    "figures/03_metacell_size_distribution.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Saved:")
print("figures/03_metacell_size_distribution.png")

# ============================================
# Step 18. Plot top 20 largest MetaCells
# ============================================

metacell_sizes = filtered.obs["metacell"].value_counts()
metacell_sizes = metacell_sizes[metacell_sizes.index != -1]

top20_metacells = metacell_sizes.head(20)

print(top20_metacells)

plt.figure(figsize=(8, 4))
top20_metacells.plot(kind="bar")

plt.xlabel("MetaCell ID")
plt.ylabel("Number of cells")
plt.title("Top 20 largest MetaCells")

plt.tight_layout()

plt.savefig(
    "figures/04_top20_largest_metacells.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Saved:")
print("figures/04_top20_largest_metacells.png")

# ============================================
# Step 19. Assigned vs Outlier cells
# ============================================

assigned = (filtered.obs["metacell"] != -1).sum()
outlier = (filtered.obs["metacell"] == -1).sum()

print(f"Assigned cells: {assigned}")
print(f"Outlier cells: {outlier}")

plt.figure(figsize=(5,5))

plt.pie(
    [assigned, outlier],
    labels=["Assigned", "Outlier"],
    autopct="%1.1f%%",
    startangle=90
)

plt.title("MetaCell Assignment Summary")

plt.savefig(
    "figures/05_metacell_assignment_summary.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Saved:")
print("figures/05_metacell_assignment_summary.png")

# ============================================
# Step 20. Cell filtering summary
# ============================================

original_cells = 30583
filtered_cells = 28759
assigned_cells = 26009

labels = [
    "Original",
    "After QC",
    "Assigned\nMetaCells"
]

values = [
    original_cells,
    filtered_cells,
    assigned_cells
]

plt.figure(figsize=(6,4))

bars = plt.bar(labels, values)

plt.ylabel("Number of cells")
plt.title("Cell Filtering Summary")

# 在柱子上标数字
for bar in bars:
    height = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width()/2,
        height + 200,
        f"{int(height):,}",
        ha="center",
        fontsize=10
    )

plt.tight_layout()

plt.savefig(
    "figures/06_cell_filtering_summary.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Saved:")
print("figures/06_cell_filtering_summary.png")

# ============================================
# Step 21. Feature gene summary
# ============================================

total_genes = filtered.n_vars

if "feature_gene" in filtered.var.columns:
    feature_genes = filtered.var["feature_gene"].sum()
else:
    feature_genes = 0

non_feature_genes = total_genes - feature_genes

print(f"Total genes: {total_genes}")
print(f"Feature genes: {feature_genes}")

plt.figure(figsize=(5,4))

bars = plt.bar(
    ["Feature genes", "Other genes"],
    [feature_genes, non_feature_genes]
)

plt.ylabel("Number of genes")
plt.title("Feature Gene Summary")

for bar in bars:
    height = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width()/2,
        height,
        f"{int(height):,}",
        ha="center",
        va="bottom"
    )

plt.tight_layout()

plt.savefig(
    "figures/07_feature_gene_summary.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Saved:")
print("figures/07_feature_gene_summary.png")

# ============================================
# Step 22. MetaCell size boxplot
# ============================================

metacell_sizes = filtered.obs["metacell"].value_counts()
metacell_sizes = metacell_sizes[metacell_sizes.index != -1]

plt.figure(figsize=(4,5))

plt.boxplot(metacell_sizes, vert=True)

plt.ylabel("Cells per MetaCell")
plt.title("MetaCell Size Boxplot")

plt.tight_layout()

plt.savefig(
    "figures/08_metacell_size_boxplot.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Saved:")
print("figures/08_metacell_size_boxplot.png")


