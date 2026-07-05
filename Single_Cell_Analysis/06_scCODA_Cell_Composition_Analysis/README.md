# Cell Composition Analysis with scCODA

## Overview

This project demonstrates Bayesian differential cell composition analysis using **scCODA** implemented in the **PertPy** framework.

The workflow explores changes in cell-type composition across biological conditions using Bayesian compositional modeling, providing a statistically robust alternative to conventional proportion-based comparisons.

---

## Objectives

- Load and inspect single-cell RNA-seq data
- Explore cell-type composition across samples
- Build a sample-by-cell-type count matrix
- Prepare data for Bayesian compositional analysis
- Perform differential abundance analysis using scCODA
- Summarize and visualize cell composition changes

---

## Dataset

- **Source:** Haber et al. (2017)
- **Species:** Mouse
- **Data type:** Single-cell RNA sequencing
- **Format:** AnnData (.h5ad)

---

## Methods

- Scanpy
- PertPy
- scCODA
- Bayesian compositional analysis
- Hamiltonian Monte Carlo (NUTS)
- Pandas
- Matplotlib

---

## Project Structure

```text
06_scCODA_Cell_Composition_Analysis
├── data/
├── figures/
├── results/
├── scripts/
└── README.md
```

---

## Outputs

- Cell type distribution
- Cell composition by sample
- Cell composition table
- Bayesian model fitting
- Differential abundance summary

---

## Learning Outcomes

- Understand compositional data analysis in single-cell transcriptomics.
- Construct sample-level count matrices for Bayesian modeling.
- Apply scCODA to detect differential cell-type abundance.
- Interpret Bayesian inference results for cell composition analysis.
