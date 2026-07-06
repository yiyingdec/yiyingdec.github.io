# RNA Velocity Analysis with scVelo

## Overview

This project demonstrates a complete RNA velocity analysis workflow using **scVelo** and **Scanpy**.

The workflow includes:

- Data loading (.loom)
- Quality control
- Normalization
- Highly variable gene selection
- RNA velocity estimation
- Velocity graph
- UMAP visualization
- Velocity stream
- Velocity grid
- Velocity embedding
- Velocity confidence
- Velocity pseudotime
- PAGA trajectory
- Dynamical model
- Latent time
- Dynamical genes

---

## Dataset

- Input format: `.loom`
- Example file: `bc.loom`

---

## Software

- Python
- Scanpy
- scVelo
- Pandas
- Matplotlib

---

## Results

### Figures

- Spliced / Unspliced proportions
- UMAP
- Velocity Stream
- Velocity Grid
- Velocity Embedding
- Single Gene Velocity
- Cell Cycle Scores
- Velocity Confidence
- Velocity Pseudotime
- PAGA
- Dynamical Velocity
- Latent Time
- Dynamical Gene Heatmap

### Output files

- `rank_velocity_genes_by_leiden.csv`
- `rank_dynamical_genes_by_leiden.csv`

---

## Project Structure

```

data/
figures/
results/
scripts/

```
