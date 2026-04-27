# Dataset Overview: PBMC3k

## Core Information
* **Provider:** 10x Genomics
* **Sample Type:** ~3,000 Peripheral Blood Mononuclear Cells (Healthy Donor)
* **Technology:** 10x Chromium v2 Chemistry
* **Genomic Mapping:** GRCh38
* **Official Resource:** [Access via 10x Genomics Portal](https://www.10xgenomics.com/datasets/3-k-pb-mc-s-from-a-healthy-donor-1-standard-1-1-0)

---

## Acquisition Methods

## Automated Retrieval (Scanpy)
For users utilizing the **Scanpy** ecosystem, the filtered gene-barcode matrix (approx. 20 MB) can be fetched and cached automatically with the following snippet:

```python
import scanpy as sc

# Fetches and loads the pbmc3k object directly
adata = sc.datasets.pbmc3k()
```

## Manual Download & Format
The native data structure utilizes the standard 10x Genomics sparse matrix convention:

| File Name | Description |
| :--- | :--- |
| **matrix.mtx** | The sparse UMI count matrix |
| **barcodes.tsv** | Unique identifiers for individual cells |
| **genes.tsv** | Mapping for gene names and Ensembl IDs |
