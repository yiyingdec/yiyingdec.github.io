# MetaCell Workflow for Single-Cell RNA-seq Analysis

## Overview

This project demonstrates a complete **MetaCell workflow** for single-cell RNA sequencing (scRNA-seq) analysis using Python.

The workflow includes:

- Data loading
- Quality control (QC)
- UMI filtering
- MetaCell construction
- Result visualization
- Export of processed datasets

The analysis was performed on a blood aging scRNA-seq dataset using the latest MetaCell package.

---

## Workflow

```text
Load Dataset
      ↓
Quality Control (QC)
      ↓
UMI Filtering
      ↓
MetaCell Construction
      ↓
Result Inspection
      ↓
Visualization
      ↓
Export Results
```

---

## Dataset

**Dataset:** Blood Aging scRNA-seq

- Original cells: **30,583**
- Genes: **33,511**

After quality control:

- High-quality cells: **28,759**

MetaCell construction:

- Assigned cells: **26,009**
- Outlier cells: **2,750**
- Final MetaCells: **916**

---

## Analysis Steps

- Import required packages
- Load AnnData object
- Explore dataset structure
- Perform quality control
- Filter cells based on UMI counts
- Construct MetaCells
- Inspect MetaCell assignments
- Save processed datasets
- Generate visualization figures

---

## Figures

```
figures/

01_UMI_distribution_before_filtering.png

02_metacell_size_distribution.png

03_top20_largest_metacells.png

04_metacell_assignment_summary.png

05_cell_filtering_summary.png

06_feature_gene_summary.png

07_metacell_size_boxplot.png
```

---

## Output Files

```
results/

blood_aging_filtered.h5ad

blood_aging_metacell.h5ad

final_metacell_results.h5ad

metacell_assignment.csv
```

---

## Key Results

- Completed quality control for scRNA-seq data.
- Filtered low-quality cells based on UMI counts.
- Constructed **916 MetaCells** from **28,759** high-quality cells.
- Successfully assigned **26,009** cells into MetaCells.
- Exported processed datasets and publication-quality figures.

---

## Software

- Python
- Scanpy
- AnnData
- MetaCell
- NumPy
- Pandas
- Matplotlib

---

## Repository Structure

```
MetaCell/
│
├── blobs/
├── figures/
├── results/
├── metacell.py
└── README.md
```

---

## Notes

This project was implemented using the latest **MetaCell** package. Some APIs differ from older tutorial versions, while the overall analytical workflow remains consistent. The workflow was adapted to ensure compatibility with the current software release.
