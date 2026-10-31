"""Paths to the input data and the results folder.

Defaults are relative to the python_code folder; override with environment variables:
  LRPS_DATA     folder holding the processed ROI time series, one sub-folder per dataset
                (Axcpt/, Cuedts/, resting/); default data/, which ships with the package
  LRPS_EVENTS   folder holding the DMCC55B run-1 *_events.tsv and *_bold.json files, one sub-folder per
                task (Axcpt/, Cuedts/), from OpenNeuro ds003465; default data/events/
  LRPS_ATLAS    folder of the AAL (SPM12) atlas: atlas/AAL.nii and ROI_MNI_V4.txt (figures only); default data/atlas/aal/
  LRPS_RESULTS  output folder
"""

from __future__ import annotations

import os
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
DATA_DIR = Path(os.environ.get("LRPS_DATA", HERE / "data"))
EVENTS_DIR = Path(os.environ.get("LRPS_EVENTS", HERE / "data" / "events"))
ATLAS_DIR = Path(os.environ.get("LRPS_ATLAS", HERE / "data" / "atlas" / "aal"))
RESULTS_DIR = Path(os.environ.get("LRPS_RESULTS", HERE / "results"))
