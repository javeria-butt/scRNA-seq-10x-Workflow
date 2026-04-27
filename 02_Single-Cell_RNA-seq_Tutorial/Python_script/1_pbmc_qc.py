import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set Scanpy settings
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=100, facecolor="white")

print("Scanpy version:", sc.__version__)

# ==========================================
# Step 1: Data Loading
# ==========================================
# Load the PBMC3k dataset. Scanpy automatically downloads it on first run.
adata = sc.datasets.pbmc3k()

print(adata)
print(f"\nShape: {adata.n_obs} cells x {adata.n_vars} genes")

# ==========================================
# Step 2: Quality Control
# ==========================================
# Basic filtering: 
# Remove cells with < 200 genes (empty droplets)
# Remove genes found in < 3 cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(f"After basic filter: {adata.shape}")

# Identify mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")

# Calculate QC metrics
sc.pp.calculate_qc_metrics(
    adata, 
    qc_vars=["mt"], 
    percent_top=None, 
    log1p=False, 
    inplace=True
)

print("QC Metrics Head:")
print(adata.obs.head())

# Visualizing QC Metrics
sc.pl.violin(
    adata, 
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4, 
    multi_panel=True, 
    save="_qc_violin.png"
)

sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", save="_scatter_mt.png")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", save="_scatter_genes.png")

# Apply QC filters based on visualization:
# 1. Remove cells with > 2500 genes (likely doublets)
# 2. Remove cells with > 5% mitochondrial counts (damaged/dying cells)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

print(f"After final QC filtering: {adata.shape}")
