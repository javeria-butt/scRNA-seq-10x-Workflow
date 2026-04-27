import scanpy as sc

# ==========================================
# Step 7: Neighbourhood Graph, UMAP, and t-SNE
# ==========================================
# Build the k-nearest neighbor (kNN) graph in PC space.
# We are using 40 principal components as determined by the PCA elbow plot.
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Compute 2D embeddings for visualization.
sc.tl.umap(adata, random_state=42)
sc.tl.tsne(adata, random_state=42)

# ==========================================
# Step 8: Leiden Clustering
# ==========================================
# Cluster the cells using the Leiden algorithm.
# Resolution 0.5 is a standard starting point for ~3,000 cells.
sc.tl.leiden(adata, resolution=0.5, random_state=42)

# Visualize the clusters on UMAP and t-SNE plots and save them.
sc.pl.umap(adata, color=["leiden"], save="_umap_leiden.png")
sc.pl.tsne(adata, color=["leiden"], save="_tsne_leiden.png")
