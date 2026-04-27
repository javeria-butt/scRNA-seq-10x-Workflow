import scanpy as sc

# ==========================================
# Step 3: Normalisation and Log-Transformation
# ==========================================
# Normalise each cell to 10,000 counts, then log1p-transform to stabilise variance.
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Save the current state (normalised counts) in .raw for later use in plotting/DE
adata.raw = adata

# ==========================================
# Step 4: Highly Variable Gene Selection
# ==========================================
# Retain genes that are both expressed and variable across cells.
sc.pp.highly_variable_genes(
    adata, 
    min_mean=0.0125, 
    max_mean=3, 
    min_disp=0.5
)

# Plot the HVGs and save
sc.pl.highly_variable_genes(adata, save="_hvg.png")

print(f"HVGs selected: {adata.var.highly_variable.sum()}")

# Subset the AnnData object to keep only the highly variable genes
adata = adata[:, adata.var.highly_variable]
