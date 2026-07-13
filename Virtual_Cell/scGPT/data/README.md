# Mouse-to-Human Gene Conversion for scGPT

## Overview

This folder contains the preprocessing workflow used to convert mouse gene symbols into their human orthologs before running scGPT.

The original dataset was derived from the mouse Paul15 single-cell RNA-seq dataset. Because the pretrained scGPT whole-human model uses a human gene vocabulary, mouse genes were mapped to their corresponding human orthologs.

The converted dataset was then used for cell embedding, cell-type annotation, and downstream scGPT analyses.

---

## Input Files

### Mouse single-cell dataset

mouse_paul15_subset.h5ad

### Mouse-to-human ortholog table

mouse_to_human_orthologs.tsv

---

## Workflow

Mouse scRNA-seq dataset

↓

Load mouse gene symbols

↓

Match mouse genes to human orthologs

↓

Remove genes without valid mappings

↓

Resolve duplicated human genes

↓

Generate a human-compatible AnnData object

↓

Save converted dataset

---

## Main Steps

1. Load the mouse AnnData object.

2. Load the ortholog mapping table.

3. Match mouse genes to human orthologs.

4. Remove unmapped genes.

5. Handle duplicated human gene symbols.

6. Replace mouse gene names with human gene names.

7. Save the converted AnnData object.

---

## Output

human_paul15_subset.h5ad

The converted dataset contains

- Human gene symbols
- Original cell annotations
- Converted AnnData object
- Gene conversion metadata

---

## Why Conversion Was Needed

The pretrained scGPT foundation model was trained using a human gene vocabulary.

Therefore, mouse genes must first be converted into their corresponding human orthologs before they can be recognized by the pretrained model.

This preprocessing step enables mouse-derived single-cell datasets to be analyzed using the pretrained human scGPT model.

---

## Software

- Python
- Scanpy
- AnnData
- pandas
- NumPy
