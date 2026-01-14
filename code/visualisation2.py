import os, glob
import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp
from nilearn import plotting, datasets
import matplotlib.pyplot as plt

csv_dir = 'your directory'
pattern = '*LRPS*.csv'   #can change to different method name       
alpha = 0.01                    

atlas = datasets.fetch_atlas_aal()
atlas_filename = atlas['maps']
coords = plotting.find_parcellation_cut_coords(atlas_filename)[:116]

# --- Load subject matrices ---
files = sorted(glob.glob(os.path.join(csv_dir, pattern)))
print(f'Found {len(files)} subject matrices.')
mats = [pd.read_csv(f, header=0).to_numpy() for f in files]
X = np.stack(mats, axis=0)     

# --- One-sample t-test vs 0 ---
tvals, pvals = ttest_1samp(X, popmean=0.0, axis=0)

# Mask with significance
T_plot = np.where(pvals < alpha, tvals, 0.0)

# Clean up diagonals (optional)
np.fill_diagonal(T_plot, 0.0)

plotting.plot_connectome(
    adjacency_matrix=T_plot,
    node_coords=coords,        
    node_color='auto',
    edge_cmap='RdBu_r',        
    node_size=15,
    edge_kwargs={'linewidth': 0.3},
    colorbar=True,
    edge_vmin=-9, 
    edge_vmax= 9
)
import matplotlib.pyplot as plt
plt.savefig("lrps_tplot.pdf", format='pdf', bbox_inches='tight')  # Save before showing
plotting.show()

