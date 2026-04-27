# 🧬 Single-Cell RNA-seq Master Repository

This repository integrates three major components of single-cell RNA-seq analysis:

1. **01_Preprocessing_of_scRNA_datasets** — Raw data → count matrix  
2. **02_Single-Cell_RNA-seq_Tutorial** — Downstream analysis (Scanpy pipeline)  
3. **03_Anndata_logic_and_memory** — Data structure & memory optimization  

---

# 📁 01_Preprocessing_of_scRNA_datasets

## 01. Pre-processing of scRNA datasets
This repository documents the complete upstream analysis of a 10X Genomics single-cell RNA-seq (scRNA-seq) dataset. The workflow transitions from raw FASTQ files to a filtered, high-quality count matrix suitable for downstream clustering and biological interpretation.

---

## Section 1: Biological Context & Library Architecture

### The Era of 1k PBMCs
The study focuses on 1,000 Peripheral Blood Mononuclear Cells (PBMCs) from a healthy donor. PBMCs are primary cells characterized by low RNA content (~1pg RNA/cell), making them an ideal model for testing the sensitivity of droplet-based sequencing.

### 10X Chromium Technology
We utilized the 10X Chromium system, which isolates single cells into nanoliter droplets (GEMs - Gel Bead-in-Emulsion).

- **Barcoding:** Each droplet contains a unique 10x barcode that indexes all transcripts from a single cell.  
- **Chemistry:** This dataset utilizes Chromium v3.  
- **Read 1 (28bp):** 16bp Cell Barcode + 12bp Unique Molecular Identifier (UMI).  
- **Read 2 (91bp):** Contains the cDNA sequence mapped to the transcriptome.  
- **Whitelist:** Barcodes were validated against the 3M-february-2018.txt whitelist to ensure reads were assigned to known, high-quality gel beads.  

---

## Section 2: Technical Methodology & Pipeline

The analysis was performed on the Galaxy platform using a sub-sampled dataset (approx. 300 cells) for computational efficiency.

### 1. Producing the Count Matrix (STARsolo)
Instead of the standard Cell Ranger pipeline, we implemented RNA STARsolo, a significantly faster "drop-in" solution that provides identical results with more transparent configuration.

- **Mapping:** Reads were aligned to the Human hg19 (GRCh37) reference genome using the Gencode/Ensembl GTF annotation.  
- **UMI Deduplication:** To remove PCR amplification bias, we applied the CellRanger 2-4 algorithm to collapse UMIs, ensuring each mRNA molecule is counted only once.  
- **Strandedness:** The library was processed as "Read strand same as original RNA molecule."  

**UCSC Alignments5_RNA STARSolo**
<img width="1340" height="584" alt="5_RNA STARSolo on 1-6_Alignments(USCS)" src="https://github.com/user-attachments/assets/934bf16e-8c11-48c2-8615-752749a15ddd" />


---

### 2. Quality Control (MultiQC)
Mapping quality was verified by aggregating STAR log files via MultiQC.

- **Key Metric:** We monitored Uniquely Mapped Reads (aiming for >75%) to verify that cDNA reads correctly aligned to the human genome without high rates of multi-mapping or unmapped noise.  

**MultiQC STAR Alignment Plot**
<img width="1600" height="800" alt="MultiQC_star_alignment_plot" src="https://github.com/user-attachments/assets/4febeaca-3a53-45c2-9dad-947576480e49" />


---

### 3. Statistical Cell Calling (DropletUtils)
The raw output from STARsolo contains approximately 5,200 barcodes, the vast majority of which represent empty droplets containing only ambient (background) RNA. We applied two methods within the DropletUtils suite to recover "real" cells:

- **The Rank Barcodes Method:** We generated a Barcode Rank Plot (Log-Total UMI vs Log-Rank). We identified the Knee point (where the UMI count drops significantly) and the Inflection point (the lower boundary of high-quality cells).  
- **The EmptyDrops Method:** We applied a Dirichlet-multinomial model to identify droplets whose RNA profile significantly differs from the ambient "background" pool.  

**Parameters:**
- Lower-bound Threshold: 200 UMIs  
- FDR: 0.01 (limiting false positives to 1%)  

**Barcode Rank PlotDropletUtils**
<img width="460" height="420" alt="DropletUtils_Barcode_Ranks_Plot" src="https://github.com/user-attachments/assets/c396ec7e-e9ea-44a1-9d60-ee426092750f" />


---

## Section 3: Technical Results & Outputs

The successful execution of the upstream pipeline produced the following technical deliverables:

- **Final Cell Count:** ~279-300 validated high-quality cells  

**Detected Cells Plot DropletUtils**
<img width="470" height="439" alt="DropletUtils_Detected_Cells_Plot" src="https://github.com/user-attachments/assets/d399c305-77db-488d-bae7-26cbdf5709b2" />


### Data Structures
- `matrix.mtx`: The digital gene expression (DGE) matrix in sparse format  
- `genes.tsv`: List of Ensembl IDs and Gene Symbols  
- `barcodes.tsv`: List of validated 10X cell barcodes  

---

## Reproducibility & Provenance
For full transparency, the complete Galaxy history, including every tool parameter used (STARsolo, MultiQC, DropletUtils), is available here:
[Galaxy link ](https://galaxy-main.usegalaxy.org/u/javeriabutt/h/scrna-seq-preprocessing)

(Developed as part of the Single Cell Galaxy Portal training.)

---

# 📊 02_Single-Cell_RNA-seq_Tutorial

## Overview
This repository documents a complete single-cell RNA-seq analysis pipeline applied to the canonical PBMC3k dataset — 2,700 peripheral blood mononuclear cells sequenced on the 10x Genomics Chromium v2 platform. The pipeline follows best practices in the Scanpy ecosystem, covering all major steps from raw count loading through to biologically-annotated cell populations.

8 distinct PBMC cell populations were identified using unsupervised Leiden clustering and validated with canonical marker gene expression.
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

---

# 🧠 03_Anndata_logic_and_memory

## AnnData Fundamentals: From Initialization to Advanced Memory Management
This repository documents the comprehensive exploration of the AnnData (Annotated Data) framework, the industry-standard Python structure for handling single-cell RNA-seq (scRNA-seq) data. The documentation follows the established methodologies of Adam Gayoso and Alex Wolf, covering the full lifecycle of a single-cell data object.

---

## 📂 Project Structure

- anndata_initialization.py: Environment setup, sparse matrix creation, and basic indexing.  
- anndata_metadata_annotation.py: Integrating aligned metadata and categorical labeling.  
- anndata_embeddings_and_uns.py: Handling multi-dimensional arrays (UMAP/PCA) and unstructured data.  
- anndata_layers_and_export.py: Managing parallel data states (Raw vs. Normalized) and DataFrame conversion.  
- anndata_disk_io.py: Persistent storage using the HDF5-based .h5ad format.  
- anndata_memory_and_backed_io.py: Advanced memory management, Views, and Backed Mode for large-scale data.  

---

## 🛠 Technical Core Concepts

### 1. The n × d Architecture
AnnData objects represent n observations (cells) and d variables (genes). The framework ensures that the central expression matrix (.X) is always synchronized with:

- `.obs`: Observation-level metadata (Cell type, batch, QC metrics)  
- `.var`: Variable-level metadata (Gene symbols, Ensembl IDs)  

---

### 2. Efficiency through Sparsity
To handle the high dropout rates in scRNA-seq, data is stored in Compressed Sparse Row (CSR) format. This dramatically reduces memory consumption by only recording non-zero values.

---

### 3. Multi-dimensional Metadata
- `.obsm`: Stores coordinates for dimensionality reductions like UMAP or t-SNE  
- `.varm`: Stores gene-related multi-dimensional data like PCA loadings  
- `.uns`: A dictionary for global metadata such as color palettes or analysis parameters  

---

### 4. Memory Optimization: Views vs. Copies
AnnData utilizes a "Lazy Copy" system. Subsetting an object creates a View (a pointer to the original data) to save RAM. Actual copies are only created through explicit calls to `.copy()` or via Implicit Modification when metadata is updated.

---

## 💾 Storage & Large-Scale Handling

### The .h5ad Format
All results are saved as `.h5ad` files, which are language-agnostic HDF5 files. This allows for:

- **Gzip Compression:** Reducing file sizes for distribution  
- **Backed Mode:** Using `ad.read_h5ad(..., backed='r')` to access datasets that exceed the system's physical RAM by keeping an active link to the disk  

---

## 📝 Troubleshooting & Environment Notes

- **Library Version:** Developed using anndata v0.8.0+  
- **Legacy Methods:** Note that `ad.read()` is replaced by `ad.read_h5ad()` in newer versions to ensure format specificity  
- **System Tools:** Inspection of `.h5ad` structures using `!h5ls` requires `hdf5-tools`  

---

## 📚 References
- Gayoso, A., & Wolf, A. (2021). Getting started with anndata.  
- The scverse Community. Exploring the AnnData object structure.  
