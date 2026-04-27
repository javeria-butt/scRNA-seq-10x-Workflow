import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

# Verify library version for reproducibility
print(f"AnnData version: {ad.__version__}")

# Generate random counts using a Poisson distribution and convert to Compressed Sparse Row (CSR) format
counts = csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32)

# Initialize the central AnnData object
adata = ad.AnnData(counts)
print(adata)

# Set unique names for cells (observations) and genes (variables)
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]

# Verify the first 10 cell identifiers
print(adata.obs_names[:10])

# Subset the object using specific observation and variable names
# This returns a 'View' to optimize memory usage
subset_adata = adata[["Cell_1", "Cell_10"], ["Gene_5", "Gene_1900"]]
print(subset_adata)
