# 01. Pre-processing of scRNA datasets

This repository documents the complete upstream analysis of a 10X Genomics single-cell RNA-seq (scRNA-seq) dataset. The workflow transitions from raw FASTQ files to a filtered, high-quality count matrix suitable for downstream clustering and biological interpretation.

---

## Section 1: Biological Context & Library Architecture

### The Era of 1k PBMCs
The study focuses on **1,000 Peripheral Blood Mononuclear Cells (PBMCs)** from a healthy donor. PBMCs are primary cells characterized by low RNA content (~1pg RNA/cell), making them an ideal model for testing the sensitivity of droplet-based sequencing.

### 10X Chromium Technology
We utilized the **10X Chromium system**, which isolates single cells into nanoliter droplets (GEMs - Gel Bead-in-Emulsion).
* **Barcoding:** Each droplet contains a unique 10x barcode that indexes all transcripts from a single cell.
* **Chemistry:** This dataset utilizes **Chromium v3**.
    * **Read 1 (28bp):** 16bp Cell Barcode + 12bp Unique Molecular Identifier (UMI).
    * **Read 2 (91bp):** Contains the cDNA sequence mapped to the transcriptome.
* **Whitelist:** Barcodes were validated against the **3M-february-2018.txt** whitelist to ensure reads were assigned to known, high-quality gel beads.

---

## Section 2: Technical Methodology & Pipeline

The analysis was performed on the **Galaxy** platform using a sub-sampled dataset (approx. 300 cells) for computational efficiency.

### 1. Producing the Count Matrix (STARsolo)
Instead of the standard Cell Ranger pipeline, we implemented **RNA STARsolo**, a significantly faster "drop-in" solution that provides identical results with more transparent configuration.
* **Mapping:** Reads were aligned to the **Human hg19 (GRCh37)** reference genome using the **Gencode/Ensembl GTF** annotation.
* **UMI Deduplication:** To remove PCR amplification bias, we applied the **CellRanger 2-4 algorithm** to collapse UMIs, ensuring each mRNA molecule is counted only once.
* **Strandedness:** The library was processed as "Read strand same as original RNA molecule."

### 2. Quality Control (MultiQC)
Mapping quality was verified by aggregating STAR log files via **MultiQC**. 
* **Key Metric:** We monitored **Uniquely Mapped Reads** (aiming for >75%) to verify that cDNA reads correctly aligned to the human genome without high rates of multi-mapping or unmapped noise.

**[MultiQC STAR Alignment Plot]**
<img width="1600" height="800" alt="MultiQC_star_alignment_plot" src="https://github.com/user-attachments/assets/832a9cd5-8617-430f-b96f-4f65095bc39a" />


### 3. Statistical Cell Calling (DropletUtils)
The raw output from STARsolo contains approximately **5,200 barcodes**, the vast majority of which represent empty droplets containing only ambient (background) RNA. We applied two methods within the **DropletUtils** suite to recover "real" cells:

* **The Rank Barcodes Method:** We generated a **Barcode Rank Plot** (Log-Total UMI vs Log-Rank). We identified the **Knee point** (where the UMI count drops significantly) and the **Inflection point** (the lower boundary of high-quality cells).
* **The EmptyDrops Method:** We applied a Dirichlet-multinomial model to identify droplets whose RNA profile significantly differs from the ambient "background" pool.
    * **Lower-bound Threshold:** 200 UMIs.
    * **FDR:** 0.01 (limiting false positives to 1%).

**[Barcode Rank Plot]**
<img width="460" height="420" alt="DropletUtils_Barcode_Ranks_Plot" src="https://github.com/user-attachments/assets/86411fea-fef0-4cb2-8add-935feaf6cde9" />


---

## Section 3: Technical Results & Outputs

The successful execution of the upstream pipeline produced the following technical deliverables:

* **Final Cell Count:** ~279-300 validated high-quality cells.
**[Detected Cells Plot]**
<img width="470" height="439" alt="DropletUtils_Detected_Cells_Plot" src="https://github.com/user-attachments/assets/56b9cf4a-3951-468d-a355-7d9ea16b2de1" />

* **Data Structures:** * `matrix.mtx`: The digital gene expression (DGE) matrix in sparse format.
    * `genes.tsv`: List of Ensembl IDs and Gene Symbols.
    * `barcodes.tsv`: List of validated 10X cell barcodes.

### Reproducibility & Provenance
For full transparency, the complete Galaxy history, including every tool parameter used (STARsolo, MultiQC, DropletUtils), is available here:

**[Link to Shared Galaxy History](https://galaxy-main.usegalaxy.org/u/javeriabutt/h/scrna-seq-preprocessing)**

---
*Developed as part of the Single Cell Galaxy Portal training.*
