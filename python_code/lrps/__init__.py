"""Low-rank pre-smoothing (LRPS) for VAR models of fMRI ROI time series.

Replication code for the simulation study (Section 3) and the fMRI experiments (Section 4).
"""

from .config import ATLAS_DIR, DATA_DIR, EVENTS_DIR, RESULTS_DIR  # noqa: F401
from .data import load_dataset  # noqa: F401
from .models import fit_lrps, fit_ols, fit_svar, lasso_path_columns, lrps_path  # noqa: F401
