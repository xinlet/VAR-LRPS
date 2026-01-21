# Robust vector autoregressive models for functional MRI activity: A low-rank pre-smoothing approach

## Abstract

Vector autoregressive (VAR) models have a history of being used to examine functional connectivity in the brain, as captured by functional MRI studies. Such VAR models allow for an estimation of Granger-causal relationships between regions-of-interest (ROI) across the brain. Unfortunately, since the number of parameters in the VAR model scales as the square of the number of ROI, and this is typically large compared to the number of temporal observations, these parameter estimates will exhibit high variance. A common approach is to regularise the VAR parameters, e.g. one may assume a sparse-VAR model where the number of interacting regions is assumed to be relatively small. In this paper, we propose an alternative low-rank pre-smoothing (LRPS) scheme to ensure the robustness of parameter estimates which performs a low-rank approximation of the observations, prior to fitting the VAR model. We fit models to individuals within a population spanning a range of tasks (including resting state data), where we demonstrate that the LRPS estimates preserve causal functional connectivity structure at the population level, whilst allowing for individual differences. This contrasts with the sparse-VAR approach, which allows individual causal relationships, but with structure that is not preserved when averaging across the population. The proposed LRPS scheme thus provides a new method to enable robust estimation of VAR parameters at the individual level, greatly reducing predictive error, whilst also enabling some interpretation at the population level. Synthetic experiments further demonstrate the effectiveness of the pre-smoothing technique.

## Authors

[Xinle Tian](https://xinlet.github.io/), [Alex Gibberd](https://sites.google.com/view/gibberd/) , [Matthew Nunes](https://people.bath.ac.uk/man54/homepage.html), [Sandipan Roy](https://sites.google.com/view/sandipanroy)

## ArXiv link
Preprint paper link can be found at [TBD]<br />

## Datasets
DMCC55B dataset can be found at [https://openneuro.org/datasets/ds003465/versions/1.0.7]<br />
Movie-watching dataset can be found at [https://openneuro.org/datasets/ds000228/versions/1.1.1]<br />
Resting-state dataset can be found at [https://https://openneuro.org/datasets/ds005003/versions/2.0.0]<br />

## Workflow
This work follows a three-stage pipeline combining Python and R:

1. **Pre-processing (Python)**
   - Raw data (`*_bold.nii.gz`) are processed in Python.
   - The AAL atlas is applied to extract region-level measurements.
   - Subjects with missing variables required for the analysis are excluded
     during this stage.
   - The preprocess code can be found at `code/preprocess.py`.
   - Outputs are saved as CSV files with predefined filenames in `data/processed/`.

2. **Estimation (R)**
   - The processed data are loaded into R.
   - Statistical estimation and model fitting are performed.
   - Demonstration code for estimation can be found at `code/estimation.R`.
   - It is worth noting that the demonstration code is organised into three main sections:
      (i) generation of VAR coefficient matrices for each method;
      (ii) evaluation of predictive performance;
      (iii) generation of cross-validation (CV) plots.
     
3. **Visualisation (Python)**
   - Results from the R estimation step are loaded in Python.
   - Figures and plots used in the paper are generated.
   - Demonstration code for visualisation can be found at `code/visualisation1.py` and `code/visualisation2.py`.
   - `code/visualisation1.py` provides multi-subject Granger causality visualisation。
   - `code/visualisation2.py` provides the group-level connectome based on one-sample t-tests of estimated VAR matrices.

## Other Visualisation

- **Rank Choice and MSPE (Simulated Data):**  
  `code/sim_mspe.R` implemented the visualisations for cross-validation rank selection and MSPE performance based on simulated datasets.

- **Principal Angle:**  
  `code/principal_angle.R` implemented the visualisation of principal angles computed between VAR matrices for each dataset.

     
