# LILRB1 Myeloid Landscape in Human Bone Marrow

Visualization and exploratory analysis workflow for mapping **LILRB1 expression** across healthy and multiple myeloma bone marrow single-cell RNA-seq datasets.

This notebook subsets selected bone marrow cohorts, reclusters them, and generates publication-style dotplots and UMAP visualizations focused on myeloid populations.

---

# Overview

The workflow:

1. Loads a large integrated bone marrow AnnData object (~3M cells)
2. Subsets selected datasets relevant to MM precursor progression
3. Reclusters the subset using Scanpy
4. Generates:
   - UMAPs
   - Global LILRB1 expression maps
   - Dotplots across disease stages
   - Myeloid-restricted visualizations
5. Saves all figures automatically as `.png`

---

# Input Data

Input object:

```python
adata_3p04M_clean_obs_050726.h5ad
```

Expected structure:

- `.obs["dataset"]`
- `.obs["preserved"]`
- `.obs["macro_cell_type_v2"]`
- `.obs["stage_raw"]`
- `.layers["counts"]`

---

# Included Datasets

The analysis subsets the following cohorts:

```python
keep_datasets = [
    "DISCO_ref",
    "GSE124310_MM_immune_RAW",
    "GSE271107_MM_precursors",
    "GSE169396_femoral_27k",
    "GSE194122_BMMC",
    "GSE253355_noanno"
]
```

These include:

- Normal bone marrow
- MGUS
- Smoldering MM
- Multiple myeloma


---

# Disease Stages

Plots are restricted to:

```python
["NBM", "MGUS", "SMM", "MM"]
```

Where:

| Label | Meaning |
|---|---|
| NBM | Normal bone marrow |
| MGUS | Monoclonal gammopathy of undetermined significance |
| SMM | Smoldering multiple myeloma |
| MM | Multiple myeloma |

---

# Reclustering Parameters

The subsetted object is reclustered using Scanpy:

```python
sc.pp.normalize_total(target_sum=1e4)
sc.pp.log1p()

sc.pp.highly_variable_genes(
    n_top_genes=3000,
    flavor="seurat_v3"
)

sc.pp.pca(n_comps=50)

sc.pp.neighbors(
    n_neighbors=15,
    n_pcs=40
)

sc.tl.umap()

sc.tl.leiden(
    resolution=0.5
)
```

---

# Dotplot Logic

For each cell type × disease stage combination:

- Dot size = % cells expressing `LILRB1`
- Dot color = % positive cells
- White → low expression
- Dark red → enriched expression

Expression positivity is defined as:

```python
counts > 0
```

---

# Plot Types

## 1. Global UMAP

Visualizes LILRB1 expression across all cells.

### Parameters

```python
cmap = "bwr"
vmax = 7.5
vmin = -7.5
```

---

## 2. Preserved Cell Type Dotplot

Uses:

```python
celltype_col = "preserved"
```

Top preserved populations ranked by:

```python
mean(% LILRB1+)
```

---

## 3. Macro Cell Type Dotplot

Uses:

```python
celltype_col = "macro_cell_type_v2"
```

Provides broader lineage-level summaries.

---

## 4. Myeloid-Restricted Macro Plot

Restricted to:

```python
[
    "Monocyte/macrophage",
    "Dendritic cell",
    "Neutrophil/granulocyte",
    "Mast cell",
    "Osteoclast",
]
```

Purpose:

- Remove lymphoid populations
- Focus specifically on innate immune remodeling

---

## 5. Myeloid-Restricted Preserved Plot

Regex-based filtering:

```python
[
    "mono",
    "macrophage",
    "dendritic",
    "dc",
    "neutrophil",
    "granulocyte",
    "mast",
    "osteoclast",
]
```

This preserves fine-grained myeloid states while hiding lymphoid compartments.

---

# Output Files

All figures are automatically saved to:

```python
/samurlab1/Joshua/ipynb_store/LILRB1/plots/
```

Example outputs:

```python
LILRB1_all_preserved_stage_raw_dotplot.png
LILRB1_all_macro_stage_raw_dotplot.png
LILRB1_macro_myeloid_stage_raw_dotplot.png
LILRB1_preserved_myeloid_stage_raw_dotplot.png
adata_sub_umap_preserved.png
adata_all_umap_LILRB1.png
```

---

# Key Helper Function

Core plotting function:

```python
lilrb1_dotplot()
```

Supports:

- Preserved labels
- Macro labels
- Myeloid-only filtering
- Adjustable figure sizing
- Automatic PNG export

---

# Dependencies

Main packages:

```python
scanpy
anndata
scvi-tools
numpy
pandas
matplotlib
scipy
```

---

# Notes

- Raw counts are preserved in:

```python
adata.layers["counts"]
```

- Dotplots use raw counts rather than normalized expression.

- Figures are intended for exploratory immune-state visualization and publication-quality refinement.

---

# Author

Joshua Rivera  
Dana-Farber Cancer Institute / Harvard Medical School
