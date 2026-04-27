import scanpy as sc

# ==========================================================
# Step 11: Dot Plot and Violin Plot — Marker Genes per Cell Type
# ==========================================================
# Defining a dictionary of specific markers to validate cell type identity
marker_dict = {
    "CD4+ T (Naive)":    ["IL7R", "CCR7", "SELL"],
    "CD14+ Monocytes":   ["CD14", "LYZ", "CST3"],
    "CD4+ T (Memory)":   ["IL7R", "S100A4"],
    "B Cells":           ["MS4A1", "CD79A"],
    "CD8+ T Cells":      ["CD8A", "CD8B"],
    "NK Cells":          ["GNLY", "NKG7"],
    "FCGR3A+ Monocytes": ["FCGR3A", "MS4A7"],
    "DC / Platelets":    ["FCER1A", "PPBP"],
}

# Generate a Dot Plot to show both expression levels and percentage of cells expressing each gene
sc.pl.dotplot(
    adata, 
    marker_dict, 
    groupby="cell_type",
    dendrogram=True, 
    save="_dotplot.png"
)

# Generate a Stacked Violin Plot for a detailed look at the distribution of expression
sc.pl.stacked_violin(
    adata, 
    marker_dict, 
    groupby="cell_type",
    rotation=90, 
    save="_stacked_violin.png"
)
