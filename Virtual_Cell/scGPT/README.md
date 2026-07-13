# scGPT Foundation Model Workflow

## Overview

This project demonstrates an end-to-end workflow using the pretrained scGPT foundation model for single-cell RNA sequencing analysis.

The project includes four major applications:

- Cell embedding and annotation
- Batch integration
- Perturbation prediction
- Gene embedding and gene-program discovery

All analyses were adapted into a low-memory workflow suitable for CPU execution on an Apple M1 MacBook Air.

---

# Workflow

Mouse dataset

↓

Mouse-to-human gene conversion

↓

Pretrained scGPT model

↓

Cell Embedding

↓

Batch Integration

↓

Perturbation Prediction

↓

Gene Embedding

↓

Gene Programs

---

# Module 1 — Cell Embedding

Dataset

- human_paul15_subset.h5ad

Workflow

Gene expression

↓

Vocabulary matching

↓

Pretrained scGPT

↓

512-dimensional embeddings

↓

Cell annotation

↓

UMAP visualization

Outputs

- Cell embedding UMAP
- Confusion matrix
- Annotated h5ad dataset

---

# Module 2 — Batch Integration

Dataset

Immune_ALL_human.h5ad

Workflow

Multi-batch immune dataset

↓

Pretrained scGPT

↓

Cell embeddings

↓

UMAP

↓

Visualization by batch

↓

Visualization by cell type

Outputs

- Batch Integration UMAP
- Integrated h5ad dataset

---

# Module 3 — Perturbation Prediction

Dataset

Adamson Perturb-seq

Selected perturbations

- BHLHE40
- CREB1

Workflow

Control cells

↓

Pretrained perturbation model

↓

Predicted gene-expression profile

↓

Comparison with experimental data

Outputs

- True vs Predicted Top20 DE genes
- Correlation scatter plot
- Perturbation prediction figures

---

# Module 4 — Gene Embedding

Workflow

Highly variable genes

↓

Extract pretrained gene embeddings

↓

Cosine similarity graph

↓

UMAP

↓

Leiden clustering

↓

Gene programs

Outputs

- Gene embedding UMAP
- Gene program summary
- Gene program CSV

---

# Technical Environment

- Python
- PyTorch
- scGPT
- Scanpy
- AnnData
- NumPy
- pandas
- Matplotlib
- Docker

---

# Memory Optimization

To allow execution on a MacBook Air without GPU, the workflow was optimized by

- CPU inference
- Small batch size
- Cell subsampling
- Explicit memory cleanup
- Garbage collection
- Low-memory data loading

---

# Skills Demonstrated

- Single-cell RNA-seq analysis
- Foundation models
- scGPT inference
- Cell embeddings
- Batch integration
- Cell annotation
- Perturbation prediction
- Gene embeddings
- Gene program discovery
- Low-memory optimization

---

# Repository Structure

scGPT/

├── data/

├── figures/

├── outputs/

├── scGPT_paul15.ipynb

└── README.md

---

# Summary

This project demonstrates how a pretrained single-cell foundation model can support multiple downstream biological tasks using a unified representation.

The same pretrained scGPT model was successfully applied to

- Cell representation learning
- Batch integration
- Perturbation prediction
- Gene representation learning

All analyses were completed on a CPU-only personal laptop using a memory-efficient workflow.
