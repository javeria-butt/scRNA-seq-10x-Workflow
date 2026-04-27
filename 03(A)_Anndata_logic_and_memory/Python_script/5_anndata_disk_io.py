#Install the tool
!apt-get install hdf5-tools
# Save the object to disk with gzip compression to minimize file size
# AnnData automatically handles the conversion of strings to categoricals during this process
adata.write('my_results.h5ad', compression="gzip")
# Using the shell command 'h5ls' to view the internal directory structure of the .h5ad file
!h5ls 'my_results.h5ad'

