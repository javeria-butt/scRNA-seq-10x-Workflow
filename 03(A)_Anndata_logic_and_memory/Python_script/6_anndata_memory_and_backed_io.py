# Create a comprehensive metadata DataFrame aligned with current observation indices
obs_meta = pd.DataFrame({
        'time_yr': np.random.choice([0, 2, 4, 8], adata.n_obs),
        'subject_id': np.random.choice(['subject 1', 'subject 2', 'subject 4', 'subject 8'], adata.n_obs),
        'instrument_type': np.random.choice(['type a', 'type b'], adata.n_obs),
        'site': np.random.choice(['site x', 'site y'], adata.n_obs),
    },
    index=adata.obs.index)

# Re-initialize AnnData to merge the new metadata with existing data
adata = ad.AnnData(adata.X, obs=obs_meta, var=adata.var)
print(adata)
# Slicing creates a 'View' (no new memory used)
adata_view = adata[:5, ['Gene_1', 'Gene_3']]
print(f"Is view: {adata_view.is_view}")

# Using .copy() to create a standalone object in memory
adata_subset = adata[:5, ['Gene_1', 'Gene_3']].copy()
# Attempting to modify a view triggers an ImplicitModificationWarning
# This ensures data integrity by creating a copy before the change is applied
adata_subset = adata[:3, ['Gene_1', 'Gene_2']]
adata_subset.obs['foo'] = range(3) # This triggers the auto-copy
# Reading in 'backed' mode for memory-intensive datasets
adata_backed = ad.read_h5ad('my_results.h5ad', backed='r')

# Verify the object is linked to the disk
print(f"Is backed: {adata_backed.isbacked}")

# Safely closing the file connection
# Note: Newer versions use .file.close() or just close the underlying handle
adata_backed.file.close()
