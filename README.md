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

UCSC Alignments5_RNA STARSolo on 1-6_Alignments(USCS)

---

### 2. Quality Control (MultiQC)
Mapping quality was verified by aggregating STAR log files via MultiQC.

- **Key Metric:** We monitored Uniquely Mapped Reads (aiming for >75%) to verify that cDNA reads correctly aligned to the human genome without high rates of multi-mapping or unmapped noise.  

MultiQC STAR Alignment PlotMultiQC_star_alignment_plot

---

### 3. Statistical Cell Calling (DropletUtils)
The raw output from STARsolo contains approximately 5,200 barcodes, the vast majority of which represent empty droplets containing only ambient (background) RNA. We applied two methods within the DropletUtils suite to recover "real" cells:

- **The Rank Barcodes Method:** We generated a Barcode Rank Plot (Log-Total UMI vs Log-Rank). We identified the Knee point (where the UMI count drops significantly) and the Inflection point (the lower boundary of high-quality cells).  
- **The EmptyDrops Method:** We applied a Dirichlet-multinomial model to identify droplets whose RNA profile significantly differs from the ambient "background" pool.  

**Parameters:**
- Lower-bound Threshold: 200 UMIs  
- FDR: 0.01 (limiting false positives to 1%)  

Barcode Rank PlotDropletUtils_Barcode_Ranks_Plot

---

## Section 3: Technical Results & Outputs

The successful execution of the upstream pipeline produced the following technical deliverables:

- **Final Cell Count:** ~279-300 validated high-quality cells  

Detected Cells PlotDropletUtils_Detected_Cells_Plot

### Data Structures
- `matrix.mtx`: The digital gene expression (DGE) matrix in sparse format  
- `genes.tsv`: List of Ensembl IDs and Gene Symbols  
- `barcodes.tsv`: List of validated 10X cell barcodes  

---

## Reproducibility & Provenance
For full transparency, the complete Galaxy history, including every tool parameter used (STARsolo, MultiQC, DropletUtils), is available here:

**Link to Shared Galaxy History**

Developed as part of the Single Cell Galaxy Portal training.

---

# 📊 02_Single-Cell_RNA-seq_Tutorial

## Overview
This repository documents a complete single-cell RNA-seq analysis pipeline applied to the canonical PBMC3k dataset — 2,700 peripheral blood mononuclear cells sequenced on the 10x Genomics Chromium v2 platform. The pipeline follows best practices in the Scanpy ecosystem, covering all major steps from raw count loading through to biologically-annotated cell populations.

8 distinct PBMC cell populations were identified using unsupervised Leiden clustering and validated with canonical marker gene expression.

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
