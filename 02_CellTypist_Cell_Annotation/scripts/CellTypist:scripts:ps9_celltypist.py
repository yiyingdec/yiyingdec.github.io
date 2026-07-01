# =============================================================================
# CellTypist Tutorial
# Part 1. Project Setup
# =============================================================================

import os
import scanpy as sc
import celltypist
from celltypist import models

# Set project directory
project_dir = "/Users/yiyingzhang/Desktop/Single cell python/CellTypist"
os.chdir(project_dir)

# Create folders
os.makedirs("data", exist_ok=True)
os.makedirs("figures", exist_ok=True)
os.makedirs("models", exist_ok=True)

# Set Scanpy figure parameters
sc.settings.figdir = "figures"
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white"
)

print(f"Working directory: {os.getcwd()}")
print(f"CellTypist version: {celltypist.__version__}")

# =============================================================================
# Part 2. Download Demo Datasets
# =============================================================================

# Read demo dataset with 2,000 cells
adata_2000 = sc.read(
    "data/demo_2000_cells.h5ad",
    backup_url="https://celltypist.cog.sanger.ac.uk/Notebook_demo_data/demo_2000_cells.h5ad"
)

# Read demo dataset with 400 cells
adata_400 = sc.read(
    "data/demo_400_cells.h5ad",
    backup_url="https://celltypist.cog.sanger.ac.uk/Notebook_demo_data/demo_400_cells.h5ad"
)

# Check data dimensions
print("adata_2000 shape:", adata_2000.shape)
print("adata_400 shape:", adata_400.shape)

# Check metadata columns
print("adata_2000 metadata columns:")
print(adata_2000.obs.columns)

print("adata_400 metadata columns:")
print(adata_400.obs.columns)

# =============================================================================
# Part 3. Download CellTypist Models
# =============================================================================

# Download all pretrained models
models.download_models(force_update=True)

# Show model cache directory
print("Model directory:")
print(models.models_path)

# Show available models
print(models.models_description())

# Load immune cell model
model = models.Model.load("Immune_All_Low.pkl")

print(model)

# =============================================================================
# Part 4. Cell Type Annotation
# =============================================================================

# Annotate cells using the pretrained immune model
predictions = celltypist.annotate(
    adata_2000,
    model="Immune_All_Low.pkl",
    majority_voting=True
)

# Convert predictions to AnnData
adata = predictions.to_adata()

# Check prediction results
print(adata.obs.head())

# =============================================================================
# Part 5. UMAP Visualization
# =============================================================================

# Compute UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Plot UMAP
sc.pl.umap(
    adata,
    color=[
        "cell_type",
        "predicted_labels",
        "majority_voting"
    ],
    legend_loc="on data",
    frameon=False,
    show=False
)

# Save figure
import matplotlib.pyplot as plt

plt.savefig(
    "figures/Figure1_CellTypist_UMAP.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# =============================================================================
# Part 6. Dot Plot Comparison
# =============================================================================

# Compare reference labels with CellTypist predictions
celltypist.dotplot(
    predictions,
    use_as_reference="cell_type",
    use_as_prediction="majority_voting"
)

# Save figure
import matplotlib.pyplot as plt

plt.savefig(
    "figures/Figure2_CellTypist_DotPlot.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# =============================================================================
# Part 7. Train a Custom CellTypist Model
# =============================================================================

# Train a custom model using the 2,000-cell dataset
custom_model = celltypist.train(
    adata_2000,
    labels="cell_type",
    n_jobs=4,
    feature_selection=True
)

# Save the custom model
custom_model.write("models/model_from_immune2000.pkl")

# Load the custom model
custom_model = models.Model.load("models/model_from_immune2000.pkl")

print(custom_model)

# =============================================================================
# Part 8. Annotate Query Cells Using the Custom Model
# =============================================================================

# Annotate the 400-cell query dataset using the custom model
custom_predictions = celltypist.annotate(
    adata_400,
    model="models/model_from_immune2000.pkl",
    majority_voting=True
)

# Convert predictions to AnnData
adata_custom = custom_predictions.to_adata()

# Compute UMAP
sc.pp.neighbors(adata_custom)
sc.tl.umap(adata_custom)

# Plot UMAP
sc.pl.umap(
    adata_custom,
    color=[
        "cell_type",
        "predicted_labels",
        "majority_voting"
    ],
    legend_loc="on data",
    frameon=False,
    show=False
)

# Save figure
import matplotlib.pyplot as plt

plt.savefig(
    "figures/Figure3_CustomModel_UMAP.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# =============================================================================
# Part 9. Marker Gene Analysis
# =============================================================================

# Load custom model
custom_model = models.Model.load("models/model_from_immune2000.pkl")

# Display available cell types
print(custom_model.cell_types)

# Extract top marker genes for Mast cells
top3_genes = custom_model.extract_top_markers(
    "Mast cells",
    3
)

print(top3_genes)

# =============================================================================
# Part 10. Marker Gene Validation
# =============================================================================

# Training dataset
sc.pl.violin(
    adata_2000,
    top3_genes,
    groupby="cell_type",
    rotation=90,
    show=False
)

plt.savefig(
    "figures/Figure4_Training_MarkerGenes.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# Query dataset
sc.pl.violin(
    adata_400,
    top3_genes,
    groupby="majority_voting",
    rotation=90,
    show=False
)

plt.savefig(
    "figures/Figure5_Query_MarkerGenes.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()





























