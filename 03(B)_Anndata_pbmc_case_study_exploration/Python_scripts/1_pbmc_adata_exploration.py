import anndata
import matplotlib.pyplot as plt
import numpy as np
import pooch
import scanpy as sc

# 1. Retrieve Dataset
datapath = pooch.retrieve(
    path=pooch.os_cache("scverse_tutorials"),
    url="https://exampledata.scverse.org/tutorials/scverse-getting-started-anndata-pbmc3k_processed.h5ad",
    known_hash="md5:b80deb0997f96b45d06f19c694e46243",
)

# 2. Load Object
adata = anndata.read_h5ad(datapath)

# 3. Structural Overview
print("--- AnnData Structure ---")
print(adata)

# 4. Active Data Matrix (.X)
# In this dataset, .X stores normalized and log-transformed counts.
print("\n--- Active Matrix (.X) ---")
print(adata.X)

# Verification of shape and format
print(f"\nShape: {adata.n_obs} cells x {adata.n_vars} genes")
print(f"Matrix Format: {type(adata.X).__name__}")
