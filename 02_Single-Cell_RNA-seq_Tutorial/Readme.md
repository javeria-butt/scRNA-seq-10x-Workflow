# 🔬 Single-Cell RNA-seq Analysis with Scanpy

**Dataset:** 3k PBMCs from a Healthy Donor (10x Genomics)  
**Tool:** Scanpy (Python-based single-cell analysis)  
**Tutorial source:** scverse/scanpy-tutorials — basic-scrna-tutorial.ipynb  
**Reference:** Wolf et al., *Genome Biology*, 2018  

---

## 📋 Table of Contents
- [Overview](#overview)
- [Pipeline Summary](#pipeline-summary)
- [Dataset](#dataset)
- [Repository Structure](#repository-structure)
- [Step-by-Step Workflow](#step-by-step-workflow)
- [Key Results](#key-results)
- [Installation](#installation)
- [Usage](#usage)
- [Discussion](#discussion)
- [References](#references)

---

## Overview
This repository documents the completion of the Scanpy basic scRNA-seq tutorial using the canonical PBMC3k dataset — 2,700 peripheral blood mononuclear cells sequenced on the 10x Genomics Chromium platform. The analysis covers the full single-cell RNA-seq downstream analysis workflow: from raw count matrix to annotated cell type clusters.

Single-cell RNA sequencing (scRNA-seq) allows the measurement of gene expression in individual cells rather than bulk tissue, revealing cellular heterogeneity that would otherwise be masked. This tutorial demonstrates how to process, cluster, and annotate single-cell data using the Scanpy ecosystem.

---
## Pipeline Summary

Raw count matrix (10x MTX format)  
│  
▼  
[AnnData object]  
│  
▼  
[Quality Control]  
│  
▼  
[Normalisation + Log transform]  
│  
▼  
[Highly Variable Genes]  
│  
▼  
[Regress out confounders]  
│  
▼  
[PCA]  
│  
▼  
[Neighbourhood graph]  
│  
▼  
[UMAP + t-SNE]  
│  
▼  
[Louvain clustering]  
│  
▼  
[Marker genes]  
│  
▼  
[Cell type annotation]  


---

## Dataset

**PBMC3k — 2,700 Peripheral Blood Mononuclear Cells from a Healthy Donor**

| Property | Value |
|----------|------|
| Cells | 2,700 (after QC: ~2,638) |
| Genes | 32,738 |
| Sequencer | Illumina NextSeq 500 |
| Reads per cell | ~69,000 |
| Chemistry | 10x Genomics v1 |
| Reference genome | hg19 (GRCh37) |

### Download the raw data
```bash
mkdir -p data
cd data
curl -O https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
```
##Repository Structure
scRNA-seq-10x-Workflow/
├── README.md
├── environment.yml
├── requirements.txt
├── CITATION.cff
├── Notebook/
│   └── pbmc3k_analysis.ipynb
├── Methods/
│   └── METHODS.md
├── Discussion/
│   └── DISCUSSION.md
├── Results/
│   ├── pbmc3k.h5ad
│   └── figures/
└── Data/

##Step-by-Step Workflow

1. Data Loading
import scanpy as sc
python
```
adata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',
    var_names='gene_symbols',
    cache=True,
)
adata.var_names_make_unique()
```

2. Quality Control
python
```
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
```

3. Normalisation
python
```
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
```
5. Feature Selection
python
```
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
```
5. Dimensionality Reduction
python
```
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.tsne(adata)
```
7. Clustering
python
```
sc.tl.louvain(adata)
```
Result: 8 clusters

7. Marker Gene Identification
python
```
sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')
```

## 8. Cell Type Annotation

| Cluster | Cell type            | Markers            |
|--------|----------------------|--------------------|
| 0      | CD4+ T cells         | IL7R, CCR7         |
| 1      | CD14+ Monocytes      | CD14, LYZ          |
| 2      | CD4+ Memory T        | IL7R, S100A4       |
| 3      | B cells              | MS4A1, CD79A       |
| 4      | CD8+ T cells         | CD8A, CD8B         |
| 5      | NK cells             | GNLY, NKG7         |
| 6      | FCGR3A+ Monocytes    | FCGR3A, MS4A7      |
| 7      | DC / Platelets       | FCER1A, PPBP       |

---

## Key Results

| Step     | Before  | After |
|----------|--------|-------|
| Cells    | 2,700  | 2,638 |
| Genes    | 32,738 | HVGs  |
| PCs      | —      | 40    |
| Clusters | —      | 8     |

##Installation

**Conda**
bash
```
conda env create -f environment.yml
conda activate scanpy-pbmc
```
**pip**
bash
```
pip install -r requirements.txt
```
##Usage##

bash
```
git clone https://github.com/javeria-butt/scRNA-seq-10x-Workflow.git
cd scRNA-seq-10x-Workflow

conda env create -f environment.yml
conda activate scanpy-pbmc

jupyter notebook Notebook/pbmc3k_analysis.ipynb
```
## Discussion

The PBMC3k dataset is a well-characterised benchmark, and the results closely match expected biology. Eight distinct immune populations were identified, including monocyte and T-cell subtypes, demonstrating the resolution of scRNA-seq.

---

## References

- Wolf, F.A. et al. (2018). *Genome Biology*  
- Satija, R. et al. (2015). *Nature Biotechnology*  
- McInnes, L. et al. (2018). *UMAP*  
- van der Maaten, L. & Hinton, G. (2008). *t-SNE*  
- 10x Genomics PBMC3k dataset  
- Scanpy documentation  


