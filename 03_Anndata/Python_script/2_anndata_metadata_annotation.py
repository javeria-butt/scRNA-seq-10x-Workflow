# Randomly assign cell types to observations
ct = np.random.choice(["B", "T", "Monocyte"], size=(adata.n_obs,))

# Store as a Categorical object for memory efficiency and faster plotting
adata.obs["cell_type"] = pd.Categorical(ct)

# Display the updated observation metadata table
print(adata.obs.head())
# Filter the AnnData object to retain only cells labeled as "B cells"
bdata = adata[adata.obs.cell_type == "B"]

# Verify that the subsetted 'View' maintains all associated metadata
print(bdata)
