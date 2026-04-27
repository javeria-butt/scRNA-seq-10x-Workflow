# Apply log(1+x) transformation to the primary matrix and store it as a new layer
# This keeps the original sparse data in .X untouched
adata.layers["log_transformed"] = np.log1p(adata.X)

# Review the updated AnnData structure with the new layer attribute
print(adata)
# Extract the log-transformed layer as a Pandas DataFrame
# Note: AnnData automatically aligns the obs_names and var_names for the DataFrame index
df_log = adata.to_df(layer="log_transformed")

# Display the first 5 rows and columns of the transformed data
print(df_log.iloc[:5, :10])
