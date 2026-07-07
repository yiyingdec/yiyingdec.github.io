# GenKI: Virtual Gene Perturbation Analysis

This folder documents a complete GenKI workflow for virtual gene perturbation analysis using the mouse Paul15 single-cell RNA-seq dataset.

The objective of this project is to learn how graph neural networks can be applied to model gene regulatory relationships and perform **in silico gene perturbation**, providing candidate genes that may influence cellular states.

---

# Project Overview

GenKI (Gene Knockout Inference) is a graph-based framework for predicting the effect of virtual gene perturbation without performing wet-lab experiments.

In this project, I completed the entire workflow from data preprocessing to model training, prediction, visualization, and result export.

---

# Objectives

The workflow includes:

1. Load single-cell RNA-seq data
2. Perform quality control and preprocessing
3. Construct a Gene Regulatory Network (GRN)
4. Train a Variational Graph AutoEncoder (VGAE)
5. Perform virtual gene perturbation
6. Rank candidate genes
7. Visualize prediction results
8. Export prediction tables

---

# Dataset

Dataset:

```text
mouse_paul15_filtered_detection005.h5ad
```

Species

- Mouse

Data type

- Single-cell RNA sequencing (scRNA-seq)

Platform

- Scanpy / AnnData (.h5ad)

---

# Workflow

```text
Single-cell RNA-seq Data
            │
            ▼
Data preprocessing
            │
            ▼
Gene filtering
            │
            ▼
Gene Regulatory Network (GRN)
            │
            ▼
Variational Graph AutoEncoder (VGAE)
            │
            ▼
Virtual Gene Perturbation
            │
            ▼
Candidate Gene Ranking
            │
            ▼
Visualization & Export
```

---

# Main Analysis Steps

## Step 1. Load data

The mouse Paul15 dataset was loaded into Scanpy.

Input file

```text
mouse_paul15_filtered_detection005.h5ad
```

---

## Step 2. Build Gene Regulatory Network

A Gene Regulatory Network (GRN) was constructed and prepared for GenKI model training.

The GRN serves as the graph structure required by the VGAE model.

---

## Step 3. Initialize GenKI

The GenKI model was initialized using the generated graph and expression matrix.

Model type

```text
VGAE
```

---

## Step 4. Train the model

The VGAE model was trained to learn latent gene representations.

Model parameters included

- Epochs
- Learning rate
- Weight decay
- Random seed

The trained model was saved as

```text
models/vgae_model_state_dict.pt
```

---

## Step 5. Virtual Gene Perturbation

After model training, GenKI predicted the effect of virtual perturbation for each candidate gene.

The prediction score represents the estimated influence of each gene on the cellular system.

---

## Step 6. Candidate Gene Ranking

Genes were ranked according to perturbation scores.

The complete prediction table was exported as

```text
tables/genki_prediction.csv
```

The cleaned Top-20 table was exported as

```text
tables/top20_genki_prediction.csv
```

---

## Step 7. Visualization

The Top 20 predicted genes were visualized using a horizontal bar chart.

Output figure

```text
figures/top20_genki_candidate_genes.png
```

---

# Results

## Top predicted genes

The highest-ranked predicted genes include

| Rank | Gene |
|------:|------|
| 1 | Gnb2l1 |
| 2 | Prkd2 |
| 3 | Rps3 |
| 4 | H2afy |
| 5 | Actb |

These genes showed the highest predicted perturbation scores in this analysis.

---

## Visualization

![Top20](figures/top20_genki_candidate_genes.png)

---

# Project Structure

```text
GenKI/
│
├── README.md
├── GenKI_mouse_paul15.ipynb
├── genki_config.json
│
├── data/
│   └── mouse_paul15_filtered_detection005.h5ad
│
├── models/
│   └── vgae_model_state_dict.pt
│
├── figures/
│   └── top20_genki_candidate_genes.png
│
└── tables/
    ├── genki_prediction.csv
    └── top20_genki_prediction.csv
```

---

# Skills Demonstrated

This project demonstrates experience with

- Single-cell RNA-seq analysis
- Scanpy
- AnnData
- Gene Regulatory Network (GRN)
- Graph Neural Networks
- Variational Graph AutoEncoder (VGAE)
- Virtual Gene Perturbation
- GenKI
- PyTorch model training
- Result visualization
- Data export
- Debugging API version incompatibilities

---

# Technical Notes

During this project, the official tutorial was not fully compatible with **GenKI 0.2.1**.

The workflow was successfully adapted by inspecting the package API and replacing deprecated interfaces with the compatible

```python
genki_model.fit()

genki_model.predict()
```

workflow.

This ensured compatibility with the latest installed GenKI version while preserving the original analysis pipeline.

---

# Future Work

This repository will be continuously expanded with additional virtual cell projects, including

- MetaCell
- CellOracle
- SCENIC
- CellPhoneDB

Together, these projects will build a comprehensive portfolio in virtual cell modeling, gene regulatory network analysis, and AI-driven computational biology.
