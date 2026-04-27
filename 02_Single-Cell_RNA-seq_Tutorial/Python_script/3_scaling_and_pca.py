import scanpy as sc

# ==========================================
# Step 5: Regression and Scaling
# ==========================================
# Remove technical confounders (total counts and mitochondrial percentage)
# then scale to unit variance with a ceiling of 10.
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)

# ==========================================
# Step 6: PCA
# ==========================================
# Compute 50 principal components using the ARPACK solver.
sc.tl.pca(adata, svd_solver="arpack", n_comps=50)

# Visualize the variance ratio to determine how many PCs to use for downstream analysis.
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save="_pca_variance.png")
