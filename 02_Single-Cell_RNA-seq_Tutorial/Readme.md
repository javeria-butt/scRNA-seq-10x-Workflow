# Single-Cell RNA-seq Analysis — PBMC


## Overview
This repository documents a complete single-cell RNA-seq analysis pipeline applied to the canonical PBMC3k dataset — 2,700 peripheral blood mononuclear cells sequenced on the 10x Genomics Chromium v2 platform. The pipeline follows best practices in the Scanpy ecosystem, covering all major steps from raw count loading through to biologically-annotated cell populations.

8 distinct PBMC cell populations were identified using unsupervised Leiden clustering and validated with canonical marker gene expression.

## Pipeline Summary

Raw Counts (2,700 cells × 32,738 genes)
↓

Quality Control (filter low-quality cells)
↓
Normalisation + Log-Transformation
↓
Highly Variable Gene Selection (~1,800 genes)
↓
Regression + Scaling
↓
PCA (50 components → top 40 selected)
↓
k-Nearest Neighbour Graph (k=10, 40 PCs)
↓
UMAP + t-SNE Embeddings
↓
Leiden Clustering (resolution=0.5)
↓
Marker Gene Identification
↓
Cell Type Annotation (8 populations)


## Step-by-Step Workflow

### 1. Data Loading
The PBMC3k dataset contains 2,700 cells and 32,738 genes in 10x Genomics sparse matrix format.

```python
import scanpy as sc
adata = sc.datasets.pbmc3k()
print(adata.shape)  # (2700, 32738)
```
###2. Quality Control

Three QC metrics are computed per cell and used to filter:
```python
Metric Threshold
n_genes_by_counts (min) > 200
n_genes_by_counts (max) < 2,500
pct_counts_mt (max) < 5%
```

Cells with fewer than 200 genes are empty droplets. Cells with >2,500 genes are likely doublets. High mitochondrial content indicates damaged or dying cells.
```python
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
```

###3. Normalisation and Log-Transformation

Library-size normalisation to 10,000 counts per cell, followed by log1p transformation:
```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # save raw for later use
```

###4. Highly Variable Gene Selection

Genes are filtered by mean expression and dispersion to retain biologically informative features:
```python
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
```

###5. Regression and Scaling

Confounding technical variation is removed, then data are scaled to unit variance:
```python
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
```

###6. PCA

50 principal components are computed; the elbow plot confirms 40 PCs are sufficient:
```python
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
```

###7. kNN Graph + UMAP + t-SNE

A k-nearest neighbour graph (k=10) is constructed in PC space, then used to generate UMAP and t-SNE embeddings:
```python
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.tsne(adata)
```

###8. Leiden Clustering

Graph-based clustering using the Leiden algorithm:
```python
sc.tl.leiden(adata, resolution=0.5)
```

###9. Marker Gene Identification and Annotation

Cell Type Key Marker Genes
CD4+ T (Naive) IL7R, CCR7, SELL, TCF7
CD4+ T (Memory) IL7R, S100A4
CD8+ T Cells CD8A, CD8B
B Cells MS4A1 (CD20), CD79A
CD14+ Monocytes CD14, LYZ, CST3
FCGR3A+ Monocytes FCGR3A, MS4A7
NK Cells GNLY, NKG7
DC / Platelets FCER1A, PPBP

##Key Results and Figures

**Figure 1 — QC Violin Plots**
<img width="2500" height="832" alt="qc_violin" src="https://github.com/user-attachments/assets/bd51ec8b-e067-4f86-b6ce-e8f4ee2673bc" />

Distribution of genes per cell, total UMI counts, and mitochondrial content across all 2,700 cells. The violin plots reveal the spread and density of each QC metric. High MT outliers (damaged cells) and cells at the lower tail (empty droplets) are clearly visible, justifying the thresholds applied.

**Figure 2 — QC Scatter Plots**
<img width="2321" height="793" alt="qc_scatter" src="https://github.com/user-attachments/assets/95816ad7-bd32-4916-85a7-f0376714a837" />

Scatter plots of total UMI counts vs. genes per cell (left) and total UMI counts vs. % MT (right). Red points are cells removed by QC filtering. The plots confirm that filtering selectively removes low-quality cells without major loss of the main cell population.

**Figure 3 — Highly Variable Gene Selection**
<img width="1780" height="793" alt="hvg_plot" src="https://github.com/user-attachments/assets/deb07e18-2254-4d87-ac0c-6253f90f56bd" />

Mean expression vs. normalised dispersion for all genes. Red points are selected HVGs. Annotated marker genes (IL7R, CD14, MS4A1, etc.) fall within the HVG set, confirming that biologically relevant genes are retained for downstream analysis.

**Figure 4 — PCA Elbow Plot**
<img width="2321" height="793" alt="pca_variance" src="https://github.com/user-attachments/assets/37703fd0-51eb-4930-97fe-a35a2be2b184" />

Variance explained per PC (left) and cumulative variance explained (right). The elbow around PC 10–15 and the cumulative curve suggest that 40 PCs capture the majority of biological variation. PCs beyond this point contribute noise rather than signal.

**Figure 5 — UMAP Coloured by Leiden Cluster**
<img width="1600" height="1146" alt="umap_louvain" src="https://github.com/user-attachments/assets/a5922b43-1483-4f34-995a-e8fdcfa9ec1d" />

UMAP embedding coloured by unsupervised Leiden cluster assignment. Clearly separated clusters indicate strong transcriptional differences between cell populations. The number of clusters identified at resolution=0.5 is consistent with the known diversity of PBMC cell types.

**Figure 6 — UMAP Coloured by Annotated Cell Type**
<img width="1781" height="1146" alt="umap_celltype" src="https://github.com/user-attachments/assets/5dec7197-7c09-4b19-a7d4-391e8ab8d8ca" />

UMAP with cells coloured by biological annotation. Eight distinct PBMC populations are clearly resolved. The separation of T cell subsets, monocyte subsets, and lymphocytes is consistent with established PBMC biology, confirming the validity of the clustering approach.

**Figure 7 — t-SNE Embedding**
<img width="1601" height="1146" alt="tsne_louvain" src="https://github.com/user-attachments/assets/d7ecf906-2399-4246-8826-5017d9759a19" />

t-SNE dimensionality reduction coloured by Leiden cluster. t-SNE preserves local neighbourhood structure and provides an alternative view consistent with the UMAP results. The concordance between UMAP and t-SNE embeddings supports the robustness of the identified clusters.

**Figure 8 — UMAP Marker Gene Expression**
<img width="3198" height="1501" alt="umap_marker_expression" src="https://github.com/user-attachments/assets/4a128d0d-59a5-417c-9838-2188d5080695" />


Expression of 8 canonical PBMC marker genes overlaid on UMAP. Each marker is specifically elevated in the expected cluster (e.g., GNLY in NK cells, CD14 in monocytes, MS4A1 in B cells), confirming the accuracy of cell type annotation.

**Figure 9 — Dot Plot**
<img width="2864" height="972" alt="marker_genes_dotplot" src="https://github.com/user-attachments/assets/d3cd2c06-768a-4c4e-9cb4-0fa054ee0ac9" />


Summary of marker gene expression across all annotated cell types. Dot size encodes the percentage of cells in each cluster that express the gene; colour encodes mean log-normalised expression. The clear diagonal enrichment pattern demonstrates high marker specificity.

**Figure 10 — Heatmap**
<img width="2731" height="1146" alt="marker_genes_heatmap" src="https://github.com/user-attachments/assets/72aab987-7207-43b6-89b5-29778a4ff85b" />


Z-score heatmap of key marker genes across 50 representative cells per cell type (sorted by annotation). The block-diagonal structure confirms that cells within each annotated group are transcriptionally coherent and distinct from other populations.


**Figure 11 — Cell Type Proportions**
<img width="2465" height="970" alt="celltype_proportions" src="https://github.com/user-attachments/assets/f0f8993b-118f-4b09-bfce-a52e88a6692d" />

Composition of the PBMC3k sample by annotated cell type. CD4+ T cells (naive + memory combined) form the largest fraction (~43%), consistent with known PBMC composition from healthy donors. CD14+ monocytes represent ~20%, and rare populations (DC/Platelets) are also captured.

