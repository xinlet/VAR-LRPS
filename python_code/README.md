# Replication code: low-rank pre-smoothing (LRPS) for VAR models of fMRI data

This folder holds the Python code for the simulation study (Section 3) and the fMRI experiments (Section 4) of the
paper. It compares three VAR(1) estimators:

- **VAR (OLS)**: ordinary least squares.
- **LRPS**: OLS regression of the rank-k SVD approximation of the response on the lagged data. The estimate is
  B_k = B_ols V_k V_k', where V_k holds the leading k right singular vectors of the training response matrix.
- **sVAR**: the sparse VAR. Each region is regressed separately on the lagged data with the lasso, using
  scikit-learn's `lasso_path` (`fit_intercept=False`, no standardisation). The objective is
  (1/2n)‖y − Xb‖² + λ‖b‖₁, which is the glmnet scaling.

All data needed to run the code are included in `data/`, so nothing has to be downloaded.

## 1. Data

Sources, licences and the references to cite are in [`data/README.md`](data/README.md).

| Folder | Contents | Source |
|---|---|---|
| `data/Axcpt/`, `data/Cuedts/` | Processed AAL-116 ROI time series, one CSV per participant | DMCC55B, OpenNeuro [ds003465](https://doi.org/10.18112/openneuro.ds003465.v1.0.7) (CC0); Etzel et al., *Scientific Data* 2022 |
| `data/resting/` | Processed AAL-116 ROI time series, one CSV per participant | OpenNeuro [ds005003](https://doi.org/10.18112/openneuro.ds005003.v2.0.0) (CC0); Leuthardt et al. |
| `data/events/Axcpt/`, `data/events/Cuedts/` | Run-1 `*_events.tsv` and `*_bold.json` (task-block timing, used by scripts 02, 04 and 05) | OpenNeuro ds003465, unmodified |
| `data/atlas/aal/` | AAL atlas for SPM12: `atlas/AAL.nii`, `ROI_MNI_V4.txt` (used for the glass-brain figures) | [GIN, Bordeaux](https://www.gin.cnrs.fr/en/tools/aal/), GNU GPL; Tzourio-Mazoyer et al., *NeuroImage* 2002 |

The ROI time series were extracted with the AAL atlas from the preprocessed fMRI data, as described in Section 4.1
of the paper.

The ROI files are named `sub-<id>_<...>_regionts.csv` and have one column per AAL region, in `ROI_MNI_V4.txt`
order. Each participant's series is z-scored per region. Participants are sorted by file name; the first 10 are
the tuning set and the rest are the evaluation set.

| Dataset | Participants | Tuning / evaluation | Time points × ROIs |
|---|---|---|---|
| Axcpt | 33 | 10 / 23 | 610 × 116 |
| Cuedts | 35 | 10 / 25 | 650 × 116 |
| Resting | 35 | 10 / 25 | 160 × 116 |

To keep any of these folders elsewhere, set these environment variables:

| Variable | Folder it points to |
|---|---|
| `LRPS_DATA` | The folder with `Axcpt/`, `Cuedts/` and `resting/` |
| `LRPS_EVENTS` | The folder with `Axcpt/` and `Cuedts/` holding the events files |
| `LRPS_ATLAS` | The AAL folder, containing `atlas/AAL.nii` and `ROI_MNI_V4.txt` |
| `LRPS_RESULTS` | The output folder (default `results/`) |

## 2. Setup

```bash
pip install -r requirements.txt
```

The code was run with Python 3.14, numpy 2.5, scipy 1.18, pandas 3.0, scikit-learn 1.9, matplotlib 3.11,
nibabel 5.4 and nilearn 0.14. Figures use the Arial font.

## 3. Running

Run the scripts from this folder, in order. `--jobs` sets the number of worker processes. Each analysis script
draws its figures when it finishes (PDF, SVG and 300 dpi TIFF in `results/figures/`), so there is no separate
plotting step; add `--no-figures` to skip them.

```bash
python 01_simulations.py --jobs 8         # Section 3:   Figs 1-3, S1, S2
python 02_block_cv.py --jobs 8            # Section 4.2: Figs 4-5, Tables 2-4
python 03_resting_expanding.py --jobs 8   # Section 4.3: Figs 6-7
python 04_connectivity.py --jobs 8        # Section 4.4: Fig 8, S3 Fig
python 05_clusters.py                     # Section 4.4: Figs 9-10
```

| Script | Approximate run time on 8 cores |
|---|---|
| `01_simulations.py` | 30 min |
| `02_block_cv.py` | 15 min |
| `03_resting_expanding.py` | 40 min (mostly the sVAR penalty path at every forecast origin, where n is close to p) |
| `04_connectivity.py`, `05_clusters.py` | 5 min each |

Scripts 01–03 are independent of each other. Script 04 uses the sVAR penalties reported in the paper; add
`--lambda-from-results` to use those selected by scripts 02 and 03.

To redraw figures from saved results without rerunning an analysis, use `python make_figures.py` (all figures),
`--step 04` (one script's figures) or `--only fig8,fig9`.

## 4. What reproduces what

| Script | Manuscript item | Main output in `results/` |
|---|---|---|
| `01_simulations.py` | Figs 1–3, S1 Fig and S2 Fig, and the results in Sections 3.2–3.3 | `simulations/summary_estimation.csv`, `summary_prediction.csv`, `selection_comparison.csv` |
| `02_block_cv.py` | Fig 4 (block timing), Fig 5, Tables 2–4 | `block_cv/table_tuning.csv` (Table 2), `table_performance.csv` (Table 3), `table_paired.csv` (Table 4) |
| `03_resting_expanding.py` | Figs 6–7, and the resting-state results in Section 4.3 | `resting/summary.json`, `table_performance.csv`, `table_paired.csv` |
| `04_connectivity.py` | Fig 8 and S3 Fig (edge counts) | `connectivity/summary_significance.csv`, `<dataset>/estimates_<method>.npz` |
| `05_clusters.py` | Figs 9–10 | `clusters/<dataset>/cluster_summary_allm.csv`, `cluster_roi_coefficients_allm.csv`, `drive_varshare_by_subject.csv`, `matching.csv` |

Every figure is written to `results/figures/` as `fig<N>_<name>` / `figS<N>_<name>`.

## 5. Code layout

| Module | Contents |
|---|---|
| `lrps/config.py` | Data and results paths (environment variables above) |
| `lrps/data.py` | Loading and z-scoring participants, the tuning/evaluation split, task blocks from `events.tsv` (trials more than 20 s apart start a new block), and within-block lag pairs |
| `lrps/models.py` | `fit_ols`, `lrps_path` / `fit_lrps`, `lrps_factors`, `lasso_path_columns` / `fit_svar` (scikit-learn), and the VAR(1) simulator |
| `lrps/simulation.py` | The true coefficient matrices and grids of Section 3 |
| `lrps/stats.py` | Paired-t comparisons, Benjamini–Hochberg, and across-participant edge tests with BH and Bonferroni corrections |
| `lrps/clusters.py` | Matched spherical k-means, Procrustes alignment, cluster summaries, variance shares, and cross-dataset matching |
| `lrps/glassbrain.py` | AAL centroids and glass-brain node and edge plots (nilearn outlines) |
| `lrps/figures.py` | All manuscript figures; `make_figures(step)` draws the figures of one analysis script |
| `lrps/utils.py` | Parallel map and the consensus (mode) rule |

**Conventions.** Lag pairs are X = (y_1, …, y_{n−1})' and Y = (y_2, …, y_n)', and the model is Y = X B + E, with
no intercept. So B = A' and B[j, i] is the directed edge from region j at time t − 1 to region i at time t.

## 6. Notes on the lasso

The lasso is solved with scikit-learn's coordinate descent (`sklearn.linear_model.lasso_path`, with tol = 1e-8).
Our original analysis code used a hand-written LARS path. For penalties λ ≥ 0.005 the two give the same
coefficients (to about 1e-9). Below that, scikit-learn reaches a slightly lower lasso objective. sVAR results
that depend on very small penalties can therefore differ marginally from the numbers in the paper. These are the
smallest values on the simulation grid (Fig 2) and the low end of the CV curves. The OLS VAR, LRPS and cluster
results do not use the lasso and reproduce exactly.
