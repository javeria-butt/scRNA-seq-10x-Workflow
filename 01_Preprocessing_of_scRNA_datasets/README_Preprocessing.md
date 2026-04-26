# 01. Pre-processing of 10X Single-Cell RNA Datasets

This sub-repository documents the **upstream processing** of raw single-cell RNA sequencing (scRNA-seq) data. The workflow covers the transition from raw sequencing reads (FASTQ) to a filtered, high-quality digital gene expression matrix.

---

## Section 1: Biological Context & Technical Specifications
The primary objective of this analysis is to resolve cellular heterogeneity within a population of immune cells.

* **Sample Data:** 1k Peripheral Blood Mononuclear Cells (PBMCs) from a healthy donor. 
* **Sequencing Technology:** **10X Genomics Chromium System**. This method uses a droplet-based approach where individual cells are encapsulated in nanoliter droplets with barcoded gel beads.
* **Library Architecture:** **Chromium v3 Chemistry**. 
    * **Read 1:** Contains the 16bp Cell Barcode (CB) and 12bp Unique Molecular Identifier (UMI).
    * **Read 2:** Contains the cDNA sequence (RNA transcript info).
* **Reference Genome:** Human **hg19 (GRCh37)** with associated Ensembl gene annotations.

## Section 2: Technical Methodology & Implementation
The upstream pipeline was executed using the **Galaxy** platform, utilizing high-performance alignment and quantification tools to handle the high-throughput nature of the 10X data.

### 1. Data Organization & QC
* **Tool:** `Galaxy Data Upload` / `MultiQC`
* **Action:** Imported sub-sampled FASTQ files (Lanes L001 and L002). Preliminary quality checks were performed to ensure base-call accuracy and identify the library chemistry (28bp Read 1 confirms v3 chemistry).

### 2. Demultiplexing, Mapping, and Quantification
* **Tool:** `RNA STARsolo`
* **Process:** * **Mapping:** Read 2 sequences were aligned to the hg19 reference.
    * **Demultiplexing:** Read 1 barcodes were matched against the **3M-february-2018 whitelist** to assign reads to specific cells.
    * **Quantification:** UMIs were deduplicated using the **CellRanger2-4 algorithm** to collapse PCR duplicates and provide an accurate count of unique mRNA molecules.

### 3. Barcode Filtering (Empty Droplet Removal)
* **Tool:** `DropletUtils`
* **Process:** A raw matrix typically contains thousands of barcodes that represent empty droplets or ambient RNA. We used the **"Rank Barcodes"** and **"EmptyDrops"** methods to identify the inflection point (knee) in the barcode rank plot, filtering the dataset down to high-quality, "real" cells.

## Section 3: Analysis Results & Data Outputs
The successful execution of the pipeline produced a structured count matrix ready for downstream biological analysis.

### Output Statistics:
* **Raw Barcodes Detected:** Approximately **5,200** initial barcodes identified by STARsolo (labeled as `yesCellBarcodes`).
* **Final Cell Count:** Approximately **300** high-quality cells (post-filtering of the sub-sampled dataset).
* **Data Format:** The results are stored in a **Matrix Market (MTX)** bundle, which includes:
    1.  `matrix.mtx`: The sparse count matrix.
    2.  `barcodes.tsv`: The unique identifiers for each validated cell.
    3.  `genes.tsv`: The list of identified gene symbols and Ensembl IDs.

### Access the Analysis:
For full transparency and reproducibility, the complete Galaxy history, including all tool parameters and intermediate steps, can be accessed at the following link:

**[View Galaxy History: scRNA-seq Preprocessing](https://galaxy-main.usegalaxy.org/u/javeriabutt/h/scrna-seq-preprocessing)**
