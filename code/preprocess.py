import os
import glob
import nibabel as nib
import numpy as np
import pandas as pd
from nilearn.input_data import NiftiLabelsMasker
from nilearn import datasets

# Input and output directories
in_dir = 'your directory'
out_dir = 'your directory'

# Create output directory if it doesn't exist
os.makedirs(out_dir, exist_ok=True)

# Fetch AAL atlas
atlas = datasets.fetch_atlas_aal()
masker = NiftiLabelsMasker(labels_img=atlas['maps'], standardize=True)

# Find all *_bold.nii.gz files
fmri_files = glob.glob(os.path.join(in_dir, '*_bold.nii.gz'))

for fmri_file in fmri_files:
    print(f"Processing {os.path.basename(fmri_file)}...")
    try:
        # Extract regional time series
        region_ts = masker.fit_transform(fmri_file)
        # Convert to DataFrame with region labels as columns
        df = pd.DataFrame(region_ts, columns=atlas['labels'])
        # Build output file name
        out_file = os.path.join(
            out_dir,
            os.path.basename(fmri_file).replace('_bold.nii.gz', '_regionts.csv')
        )
        # Save to CSV
        df.to_csv(out_file, index=False)
        print(f"Saved region time series to {out_file}")
    except Exception as e:
        print(f"Error processing {fmri_file}: {e}")
