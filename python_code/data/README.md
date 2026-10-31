# Data: sources and attribution

Everything the code needs is in this folder. Please cite the original sources below when you use these data.

| Folder | Contents | Source | Licence |
|---|---|---|---|
| `Axcpt/`, `Cuedts/` | Processed ROI time series (AAL, 116 regions), one CSV per participant: 33 Axcpt and 35 Cuedts participants, baseline session, run 1 | DMCC55B, OpenNeuro ds003465 | CC0 |
| `resting/` | Processed ROI time series (AAL, 116 regions), one CSV per participant: 35 participants | Retrospective Task/Rest fMRI Data from Tumor Patients, OpenNeuro ds005003 | CC0 |
| `events/Axcpt/`, `events/Cuedts/` | Unmodified `*_events.tsv` and `*_bold.json` files for the same runs, used for the task-block timing | DMCC55B, OpenNeuro ds003465 | CC0 |
| `atlas/aal/` | AAL atlas for SPM12, version 4 (25/08/2015): `atlas/AAL.nii`, `atlas/AAL.xml`, `ROI_MNI_V4.txt` and the distribution's readme, all unmodified | Groupe d'Imagerie Neurofonctionnelle (GIN), Bordeaux | GNU General Public License (see `atlas/aal/readme_aal_for_SPM12.txt`) |

The ROI time series were extracted from the preprocessed fMRI data with the AAL atlas, as described in Section 4.1
of the paper. Each CSV has one column per AAL region, in `ROI_MNI_V4.txt` order.

## DMCC55B (task data: `Axcpt/`, `Cuedts/`, `events/`)

- Dataset: Braver TS, Kizhner A, Tang R, Freund MC, Etzel JA. *DMCC55B*. OpenNeuro, ds003465, version 1.0.7.
  <https://openneuro.org/datasets/ds003465/versions/1.0.7>, doi:[10.18112/openneuro.ds003465.v1.0.7](https://doi.org/10.18112/openneuro.ds003465.v1.0.7)
- The events and bold.json files were downloaded from OpenNeuro's public S3 copy of ds003465
  (`https://s3.amazonaws.com/openneuro.org/ds003465/<sub>/ses-wave1bas/func/<file>`).
- The dataset asks users to cite:
  - Etzel JA, Brough RE, Freund MC, Kizhner A, Lin Y, Singh MF, et al. The Dual Mechanisms of Cognitive Control
    dataset, a theoretically-guided within-subject task fMRI battery. *Scientific Data* 2022; 9: 114.
    Preprint: doi:[10.1101/2021.05.28.446178](https://doi.org/10.1101/2021.05.28.446178).
  - doi:[10.1162/jocn_a_01768](https://doi.org/10.1162/jocn_a_01768).

## Resting state (`resting/`)

- Dataset: Leuthardt EC, Shimony JS, Snyder AZ, Dierker D, Park KY. *Retrospective Task/Rest fMRI Data from Tumor
  Patients*. OpenNeuro, ds005003, version 2.0.0.
  <https://openneuro.org/datasets/ds005003/versions/2.0.0>, doi:[10.18112/openneuro.ds005003.v2.0.0](https://doi.org/10.18112/openneuro.ds005003.v2.0.0)
- The dataset asks users to cite doi:[10.1016/j.nic.2017.06.011](https://doi.org/10.1016/j.nic.2017.06.011)
  and doi:[10.1371/journal.pone.0236423](https://doi.org/10.1371/journal.pone.0236423). Its suggested
  acknowledgement is: "Data was sourced from the OpenNeuro database, accession number ds005003".

## AAL atlas (`atlas/aal/`)

- Source: Anatomical Automatic Labeling (AAL) for SPM12, Groupe d'Imagerie Neurofonctionnelle (GIN, UMR5296,
  Bordeaux). <https://www.gin.cnrs.fr/en/tools/aal/>, archive `aal_for_SPM12.tar.gz`.
- It is distributed as copyright freeware under the GNU General Public License.
- Reference: Tzourio-Mazoyer N, Landeau B, Papathanassiou D, Crivello F, Etard O, Delcroix N, Mazoyer B, Joliot M.
  Automated anatomical labeling of activations in SPM using a macroscopic anatomical parcellation of the MNI MRI
  single-subject brain. *NeuroImage* 2002; 15: 273–289.
- The cerebellar regions follow Schmahmann JD et al., *NeuroImage* 1999; 10: 233–260.
