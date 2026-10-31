# Prediction of fMRI activity using vector autoregressive models: a comparison of sparse and low-rank approaches

## Abstract

Vector autoregressive (VAR) models have a history of being used to examine functional connectivity in the brain, as captured by functional MRI studies. Such models allow for an estimation of Granger-causal relationships between regions of interest across the brain. Unfortunately, since the number of parameters in the VAR model scales as the square of the number of regions, and this is typically large compared to the number of temporal observations, these parameter estimates will exhibit high variance. To address this challenge, we introduce a low-rank pre-smoothing method that applies a low-rank approximation to the observations before fitting a VAR model. We estimate these models using individual subject data from both task-based and resting-state conditions, tuning hyperparameters at the population level. Our low-rank approach is directly compared against sparse and unconstrained estimation methods. Evaluations of predictive performance and model structure reveal that our pre-smoothing technique enables robust individual-level parameter estimation and significantly reduces prediction error, a finding further validated by synthetic experiments where the ground-truth parameters are known.

## Authors

[Xinle Tian](https://xinlet.github.io/), [Alex Gibberd](https://sites.google.com/view/gibberd/), [Sandipan Roy](https://researchportal.bath.ac.uk/en/persons/sandipan-roy/), [Matthew Nunes](https://mattnunes.github.io/) 

## BioRxiv link
Preprint paper link can be found at [https://www.biorxiv.org/content/10.64898/2026.06.11.731556v1.abstract]<br />

## Datasets
DMCC55B dataset can be found at [https://openneuro.org/datasets/ds003465/versions/1.0.7]<br />
Resting-state dataset can be found at [https://https://openneuro.org/datasets/ds005003/versions/2.0.0]<br />

## Explore the project

| File or folder | Purpose |
|---|---|
| [README.md](README.md) | Detailed setup, implementation notes, script arguments and manuscript-to-output mapping. |
| `01_simulations.py` | Synthetic coefficient-estimation and prediction experiments. |
| `02_block_cv.py` | Task-based forecasting with one or two training blocks. |
| `03_resting_expanding.py` | Resting-state forecasting as the training sample grows. |
| `04_connectivity.py` | Across-participant coefficient summaries and adjusted connection tests. |
| `05_clusters.py` | LRPS network clustering and regional driving-pattern analysis. |
| `lrps/` | Shared model, data, statistics and plotting functions. |
| `results/` | Saved numerical results and manuscript figures. |

To inspect the findings without rerunning the analyses, start with `results/figures/`, the task tables in `results/block_cv/`, and the resting-state tables in `results/resting/`.

## Getting started

Install the dependencies from this directory:

```bash
python -m pip install -r requirements.txt
```

To redraw the figures from the supplied results:

```bash
python make_figures.py
```

The existing [technical README](README.md) explains how to run individual experiments, configure paths and choose the number of worker processes. Analyses and figure generation write to `results/` by default. Set `LRPS_RESULTS` to a separate output directory when running new experiments if you want to preserve the supplied results.

## Attribution

Please cite the accompanying manuscript and the original datasets and atlas when using this project. Dataset and atlas attribution details are collected in [data/README.md](data/README.md). Those third-party terms are separate from any licence chosen for the project source code.

     
