# Adding a synthetic 2D UMAP embedding to observation metadata
adata.obsm["X_umap"] = np.random.normal(0, 1, size=(adata.n_obs, 2))

# Adding 5-dimensional random metadata to the variables (genes)
adata.varm["gene_stuff"] = np.random.normal(0, 1, size=(adata.n_vars, 5))

# Verification of multi-dimensional keys
print(f"Observation matrices: {adata.obsm}")
print(f"Variable matrices: {adata.varm}")
# Storing a random list as general metadata
adata.uns["random"] = [1, 2, 3]

# Displaying the unstructured metadata dictionary
print(f"Unstructured metadata: {adata.uns}")
