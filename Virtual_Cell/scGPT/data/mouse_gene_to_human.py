# ============================================================
# Mouse-to-Human Gene Symbol Conversion
# ============================================================

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


# ------------------------------------------------------------
# File paths
# ------------------------------------------------------------

PROJECT_DIR = Path("/Users/yiyingzhang/Desktop/scGPT")

INPUT_H5AD = PROJECT_DIR / "mouse_paul15_subset.h5ad"
ORTHOLOG_FILE = PROJECT_DIR / "mouse_to_human_orthologs.tsv"
OUTPUT_H5AD = PROJECT_DIR / "human_paul15_subset.h5ad"


# ------------------------------------------------------------
# Check whether the required files exist
# ------------------------------------------------------------

if not INPUT_H5AD.exists():
    raise FileNotFoundError(
        f"Input h5ad file was not found:\n{INPUT_H5AD}"
    )

if not ORTHOLOG_FILE.exists():
    raise FileNotFoundError(
        f"Ortholog mapping file was not found:\n{ORTHOLOG_FILE}"
    )

print("Required files were found successfully.")
print(f"Input dataset: {INPUT_H5AD}")
print(f"Ortholog table: {ORTHOLOG_FILE}")


# ------------------------------------------------------------
# Load the mouse single-cell dataset
# ------------------------------------------------------------

print("\nLoading the mouse single-cell dataset...")

adata = sc.read_h5ad(INPUT_H5AD)

print(f"Original dataset shape: {adata.shape}")
print(f"Number of cells: {adata.n_obs}")
print(f"Number of mouse genes: {adata.n_vars}")

print("\nFirst 10 original gene names:")
print(adata.var_names[:10].tolist())


# ------------------------------------------------------------
# Load the mouse-to-human ortholog table
# ------------------------------------------------------------

print("\nLoading the mouse-to-human ortholog table...")

ortholog = pd.read_csv(
    ORTHOLOG_FILE,
    sep="\t",
    dtype=str
)

print(f"Ortholog table shape: {ortholog.shape}")

print("\nColumns detected in the ortholog table:")
print(ortholog.columns.tolist())

print("\nFirst five rows of the ortholog table:")
print(ortholog.head())


# ------------------------------------------------------------
# Detect the mouse and human gene columns
# ------------------------------------------------------------

mouse_column_candidates = [
    "mouse_gene",
    "mouse_symbol",
    "mouse_gene_symbol",
    "Mouse_gene",
    "Mouse_symbol",
    "Mouse gene",
    "Mouse gene symbol",
    "mouse",
    "Mouse"
]

human_column_candidates = [
    "human_gene",
    "human_symbol",
    "human_gene_symbol",
    "Human_gene",
    "Human_symbol",
    "Human gene",
    "Human gene symbol",
    "human",
    "Human"
]


def detect_column(columns, candidates, species_name):
    """
    Detect the most likely gene-symbol column from a list of candidates.
    """

    normalized_columns = {
        str(column).strip().lower().replace("-", "_").replace(" ", "_"): column
        for column in columns
    }

    for candidate in candidates:
        normalized_candidate = (
            candidate.strip()
            .lower()
            .replace("-", "_")
            .replace(" ", "_")
        )

        if normalized_candidate in normalized_columns:
            return normalized_columns[normalized_candidate]

    for column in columns:
        normalized_column = str(column).strip().lower()

        if (
            species_name.lower() in normalized_column
            and (
                "gene" in normalized_column
                or "symbol" in normalized_column
            )
        ):
            return column

    return None


mouse_column = detect_column(
    ortholog.columns,
    mouse_column_candidates,
    "mouse"
)

human_column = detect_column(
    ortholog.columns,
    human_column_candidates,
    "human"
)


if mouse_column is None or human_column is None:
    raise ValueError(
        "\nThe mouse and human gene columns could not be detected automatically.\n"
        f"Available columns: {ortholog.columns.tolist()}\n"
        "Please inspect the ortholog table column names."
    )

print("\nDetected mapping columns:")
print(f"Mouse gene column: {mouse_column}")
print(f"Human gene column: {human_column}")


# ------------------------------------------------------------
# Clean the ortholog mapping table
# ------------------------------------------------------------

mapping = ortholog[
    [mouse_column, human_column]
].copy()

mapping.columns = [
    "mouse_gene",
    "human_gene"
]

mapping["mouse_gene"] = (
    mapping["mouse_gene"]
    .astype(str)
    .str.strip()
)

mapping["human_gene"] = (
    mapping["human_gene"]
    .astype(str)
    .str.strip()
)

mapping = mapping.replace(
    {
        "": np.nan,
        "nan": np.nan,
        "None": np.nan
    }
)

mapping = mapping.dropna(
    subset=["mouse_gene", "human_gene"]
)

mapping = mapping.drop_duplicates(
    subset=["mouse_gene"],
    keep="first"
)

print(f"\nNumber of valid ortholog pairs: {len(mapping)}")


# ------------------------------------------------------------
# Create the mouse-to-human mapping dictionary
# ------------------------------------------------------------

mouse_to_human = dict(
    zip(
        mapping["mouse_gene"],
        mapping["human_gene"]
    )
)


# ------------------------------------------------------------
# Map mouse genes to human orthologs
# ------------------------------------------------------------

adata.var_names = adata.var_names.astype(str)

adata.var["mouse_gene"] = adata.var_names

adata.var["human_gene"] = (
    adata.var["mouse_gene"]
    .map(mouse_to_human)
)

matched_genes = int(
    adata.var["human_gene"].notna().sum()
)

unmatched_genes = int(
    adata.var["human_gene"].isna().sum()
)

matching_rate = (
    matched_genes / adata.n_vars * 100
    if adata.n_vars > 0
    else 0
)

print("\nGene mapping summary:")
print(f"Matched genes: {matched_genes}")
print(f"Unmatched genes: {unmatched_genes}")
print(f"Matching rate: {matching_rate:.2f}%")


if matched_genes == 0:
    raise ValueError(
        "\nNo mouse genes were matched to human orthologs.\n"
        "The gene identifiers in the h5ad file may not match "
        "the identifiers in the ortholog table."
    )


# ------------------------------------------------------------
# Remove genes without human orthologs
# ------------------------------------------------------------

matched_mask = adata.var["human_gene"].notna().to_numpy()

adata = adata[:, matched_mask].copy()

adata.var_names = (
    adata.var["human_gene"]
    .astype(str)
    .to_numpy()
)

print(
    "\nDataset shape after removing unmatched genes: "
    f"{adata.shape}"
)


# ------------------------------------------------------------
# Merge duplicated human gene symbols
# ------------------------------------------------------------

duplicated_gene_count = int(
    adata.var_names.duplicated().sum()
)

print(
    "\nNumber of duplicated human gene symbols: "
    f"{duplicated_gene_count}"
)


if duplicated_gene_count > 0:

    print("Merging duplicated human gene symbols by summing expression values...")

    human_genes = pd.Index(
        adata.var_names.astype(str)
    )

    unique_human_genes = pd.Index(
        human_genes.unique()
    )

    gene_to_index = {
        gene: index
        for index, gene in enumerate(unique_human_genes)
    }

    new_gene_indices = np.array(
        [
            gene_to_index[gene]
            for gene in human_genes
        ],
        dtype=int
    )

    aggregation_matrix = sp.csr_matrix(
        (
            np.ones(adata.n_vars, dtype=np.float32),
            (
                np.arange(adata.n_vars),
                new_gene_indices
            )
        ),
        shape=(
            adata.n_vars,
            len(unique_human_genes)
        )
    )

    if sp.issparse(adata.X):
        merged_x = adata.X @ aggregation_matrix
    else:
        merged_x = np.asarray(adata.X) @ aggregation_matrix

    merged_adata = sc.AnnData(
        X=merged_x,
        obs=adata.obs.copy()
    )

    merged_adata.var_names = unique_human_genes.astype(str)

    merged_adata.var["human_gene"] = (
        merged_adata.var_names.astype(str)
    )

    adata = merged_adata

    print("Duplicated genes were merged successfully.")

else:
    print("No duplicated human gene symbols were found.")


# ------------------------------------------------------------
# Ensure unique human gene names
# ------------------------------------------------------------

adata.var_names_make_unique()

adata.uns["gene_conversion"] = {
    "source_species": "mouse",
    "target_species": "human",
    "input_file": str(INPUT_H5AD),
    "ortholog_file": str(ORTHOLOG_FILE),
    "matched_mouse_genes": matched_genes,
    "unmatched_mouse_genes": unmatched_genes,
    "matching_rate_percent": float(matching_rate)
}


# ------------------------------------------------------------
# Save the converted human-gene dataset
# ------------------------------------------------------------

print("\nSaving the converted dataset...")

adata.write_h5ad(
    OUTPUT_H5AD,
    compression="gzip"
)

print("\nConversion completed successfully.")
print(f"Final dataset shape: {adata.shape}")
print(f"Output file: {OUTPUT_H5AD}")

print("\nFirst 20 human gene symbols:")
print(adata.var_names[:20].tolist())
