import os
import glob
import numpy as np
import pandas as pd

# Directory containing your estimate CSVs
csv_dir = 'your directory'

# Find all CSVs
lrps_files = glob.glob(os.path.join(csv_dir, '*LRPS*.csv')) #you can change this to switch between method
print(f"Found {len(lrps_files)} LRPS files.")

# Read all as numpy arrays, store in a list
matrices = [pd.read_csv(f, header=0).to_numpy() for f in lrps_files]

# Check all have the same shape
shapes = [m.shape for m in matrices]
assert all(s == shapes[0] for s in shapes), "Not all LRPS matrices have the same shape!"

# Compute element-wise average
G = np.mean(matrices, axis=0)


import numpy as np
from nilearn import plotting, datasets
import copy
import pandas as pd

# Fetch AAL atlas and coordinates for first 113 regions
atlas = datasets.fetch_atlas_aal()
atlas_filename = atlas['maps']
coords = plotting.find_parcellation_cut_coords(atlas_filename)[:116]

# Set a percentile threshold (e.g., keep top 1% of all values)
percentile = 99  # keep connections above 99th percentile

# Work with absolute values for GC strength
G_abs = np.abs(G)
thresh = np.percentile(G_abs, percentile)
G_threshold = np.copy(G)
G_threshold[G_abs < thresh] = 0

print(f"Thresholded at {thresh:.4f}; number of connections kept: {np.count_nonzero(G_threshold)}")

# Plot
plotting.plot_connectome(
    G_threshold,
    coords,
    node_color='auto',
    display_mode='lyrz',
    colorbar=True,
    edge_kwargs={'linewidth': 0.1},
    node_size=13,
    edge_vmin=-0.8,  
    edge_vmax=0.8  
)
import matplotlib.pyplot as plt
plt.savefig("connectome_plot.pdf", format='pdf', bbox_inches='tight')  # Save before showing
plotting.show()
