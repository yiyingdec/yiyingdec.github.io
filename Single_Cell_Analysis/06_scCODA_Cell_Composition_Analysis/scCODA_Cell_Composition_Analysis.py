# ============================================
# Step 1. Import packages and set directories
# ============================================

import os
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")

import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import pertpy as pt

# Set project directories
project_dir = Path.home() / "Desktop" / "Single cell python" / "scCODA"

data_dir = project_dir / "data"
figure_dir = project_dir / "figures"
result_dir = project_dir / "results"

figure_dir.mkdir(exist_ok=True)
result_dir.mkdir(exist_ok=True)

os.chdir(project_dir)

print("PertPy version:", pt.__version__)
print("Scanpy version:", sc.__version__)
print("Project directory:", project_dir)


# ============================================
# Step 2. Load Haber 2017 AnnData object
# ============================================

adata = sc.read_h5ad(data_dir / "haber_2017_regions.h5ad")

print(adata)
print("\nObservation columns:")
print(adata.obs.columns.tolist())

print("\nObservation metadata:")
print(adata.obs.head())


# ============================================
# Step 3. Explore cell composition
# ============================================

print("Condition counts:")
print(adata.obs["condition"].value_counts())

print("\nCell type counts:")
print(adata.obs["cell_label"].value_counts())

plt.figure(figsize=(10, 6))

adata.obs["cell_label"].value_counts().plot(kind="bar")

plt.title("Cell Type Distribution")
plt.xlabel("Cell Type")
plt.ylabel("Number of Cells")
plt.xticks(rotation=45, ha="right")

plt.tight_layout()

plt.savefig(
    figure_dir / "01_Cell_Type_Distribution.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ============================================
# Step 4. Create cell composition table
# ============================================

composition_df = pd.crosstab(
    adata.obs["batch"],
    adata.obs["cell_label"]
)

sample_metadata = (
    adata.obs[["batch", "condition"]]
    .drop_duplicates()
    .set_index("batch")
)

composition_df = composition_df.join(sample_metadata)

composition_df.to_csv(result_dir / "cell_composition_table.csv")

print("\nCell composition table:")
print(composition_df)

print("\nComposition table shape:")
print(composition_df.shape)


# ============================================
# Step 5. Prepare data for scCODA
# ============================================

cell_counts = composition_df.drop(columns="condition")
sample_metadata = composition_df[["condition"]]

print("\nCell count matrix:")
print(cell_counts.head())

print("\nSample metadata:")
print(sample_metadata.head())

print("\nCell count matrix shape:", cell_counts.shape)
print("Metadata shape:", sample_metadata.shape)


# ============================================
# Step 6. Create scCODA data object
# ============================================

from pertpy.tools import Sccoda

sccoda = Sccoda()

coda_data = sccoda.load(
    adata,
    type="cell_level",
    generate_sample_level=True,
    cell_type_identifier="cell_label",
    sample_identifier="batch",
    covariate_obs=["condition"],
)

print("\nscCODA data object:")
print(coda_data)

# ============================================
# Step 7. Visualize sample-level composition
# ============================================

# Extract sample-level count table from scCODA object
coda_counts = coda_data["coda"].X
coda_cell_types = coda_data["coda"].var_names
coda_samples = coda_data["coda"].obs_names

coda_counts_df = pd.DataFrame(
    coda_counts,
    index=coda_samples,
    columns=coda_cell_types
)

# Convert counts to proportions
coda_prop_df = coda_counts_df.div(coda_counts_df.sum(axis=1), axis=0)

# Add condition information
coda_prop_df["condition"] = coda_data["coda"].obs["condition"].values

# Save proportion table
coda_prop_df.to_csv(result_dir / "cell_composition_proportions.csv")

# Plot stacked barplot
plot_df = coda_prop_df.drop(columns="condition")

plt.figure(figsize=(10, 6))
plot_df.plot(
    kind="bar",
    stacked=True,
    figsize=(10, 6)
)

plt.title("Cell Composition by Sample")
plt.xlabel("Sample")
plt.ylabel("Proportion")
plt.xticks(rotation=45, ha="right")
plt.legend(title="Cell Type", bbox_to_anchor=(1.05, 1), loc="upper left")

plt.tight_layout()

plt.savefig(
    figure_dir / "02_Cell_Composition_by_Sample.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print(coda_prop_df.head())

# ============================================
# Step 8. Prepare compositional model
# ============================================

sccoda.prepare(
    coda_data,
    formula="condition",
    modality_key="coda",
)

print(coda_data["coda"])
print("\nscCODA model preparation finished.")

# ============================================
# Step 9. Run scCODA model
# ============================================

# Run Hamiltonian Monte Carlo sampling
sccoda.run_nuts(
    coda_data,
    modality_key="coda",
)

print("scCODA model fitting finished.")

# ============================================
# Step 10. Summarize results
# ============================================

results = sccoda.summary(
    coda_data,
    modality_key="coda"
)

print(results)

# Save summary
with open(result_dir / "03_scCODA_summary.txt", "w") as f:
    f.write(str(results))

print("Results saved.")


