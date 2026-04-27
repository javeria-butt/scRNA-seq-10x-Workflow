import scanpy as sc

# ==========================================
# Step 9: Marker Gene Expression on UMAP
# ==========================================
# List of known marker genes for common PBMC cell types
marker_genes = [
    "IL7R", "CCR7", "CD14", "LYZ", "MS4A1",
    "CD8A", "GNLY", "NKG7", "FCGR3A", "FCER1A",
    "CST3", "PPBP"
]

# Visualize marker expression to help verify cluster identities
# vmax="p99" clips the top 1% of expression to improve contrast
sc.pl.umap(
    adata, 
    color=marker_genes, 
    ncols=4, 
    vmax="p99", 
    save="_umap_markers.png"
)

# ==========================================
# Step 10: Cell Type Annotation
# ==========================================
# Map leiden clusters to biological cell types based on marker expression observed above
cell_type_map = {
    "0": "CD4+ T (Naive)",
    "1": "CD14+ Monocytes",
    "2": "CD4+ T (Memory)",
    "3": "B Cells",
    "4": "CD8+ T Cells",
    "5": "NK Cells",
    "6": "FCGR3A+ Monocytes",
    "7": "DC / Platelets",
    # Note: If your Leiden run produced more than 8 clusters, 
    # they will be marked as "Unknown" via .fillna()
}

# Assign the biological labels to a new column in the observations (metadata)
adata.obs["cell_type"] = adata.obs["leiden"].map(cell_type_map).fillna("Unknown")

# Final UMAP visualization with biological labels
sc.pl.umap(
    adata, 
    color=["cell_type"], 
    legend_loc="on data",
    legend_fontsize=9, 
    save="_umap_celltypes.png"
)
