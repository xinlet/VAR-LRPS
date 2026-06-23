# Prediction of fMRI activity using vector autoregressive models: a comparison of sparse and low-rank approaches

## Abstract

Vector autoregressive (VAR) models have a history of being used to examine functional connectivity in the brain, as captured by functional MRI studies. Such models allow for an estimation of Granger-causal relationships between regions of interest across the brain. Unfortunately, since the number of parameters in the VAR model scales as the square of the number of regions, and this is typically large compared to the number of temporal observations, these parameter estimates will exhibit high variance. To address this challenge, we introduce a low-rank pre-smoothing method that applies a low-rank approximation to the observations before fitting a VAR model. We estimate these models using individual subject data from both task-based and resting-state conditions, tuning hyperparameters at the population level. Our low-rank approach is directly compared against sparse and unconstrained estimation methods. Evaluations of predictive performance and model structure reveal that our pre-smoothing technique enables robust individual-level parameter estimation and significantly reduces prediction error, a finding further validated by synthetic experiments where the ground-truth parameters are known.

## Authors

[Xinle Tian](https://xinlet.github.io/), [Alex Gibberd](https://sites.google.com/view/gibberd/), [Sandipan Roy](https://researchportal.bath.ac.uk/en/persons/sandipan-roy/), [Matthew Nunes](https://people.bath.ac.uk/man54/homepage.html) 

## BioRxiv link
Preprint paper link can be found at [https://www.biorxiv.org/content/10.64898/2026.06.11.731556v1.abstract]<br />

## Datasets
DMCC55B dataset can be found at [https://openneuro.org/datasets/ds003465/versions/1.0.7]<br />
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

     
