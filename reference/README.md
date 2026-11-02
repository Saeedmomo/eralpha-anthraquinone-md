# Reference data and V9 ensemble scope

The original six reference files are retained unchanged. The antagonist files
without an ensemble-size suffix contain the original balanced comparison: three
trajectories for each of eight ligands, 24 trajectories in total. This comparison
underlies manuscript Table 7 and SI Table S8 and Figure S5. The two agonist files
contain nine trajectories, three for each of three ligands.

| Files | Ensemble | Residue rows |
| --- | --- | --- |
| `fp_project_by_residue.csv` | Balanced 24-run antagonist | 482 |
| `fp_github_by_residue.csv` | Balanced 24-run antagonist | 555 |
| `fp_project_by_type.csv`, `fp_github_by_type.csv` | Balanced 24-run antagonist; one row per trajectory | Not applicable |
| `fp_agonist_project_by_residue.csv` | Nine-run agonist with TIF2 peptide | 142 |
| `fp_agonist_github_by_residue.csv` | Nine-run agonist with TIF2 peptide | 176 |
| `fp_antagonist25_project_by_residue.csv` | Updated 25-run antagonist | 501 |
| `fp_antagonist25_github_by_residue.csv` | Updated 25-run antagonist | 577 |

The updated antagonist files append the original extra-replicate outputs for
`anthraflavic_r5` to the balanced files, preserving every input field value. That
source directory is reported as the retained 2,6-DHAQ replicate R4 in V9.
`anthraflavic_r2`, the retained run with low direct bridging, remains included.
The excluded original `anthraflavic_r4` is absent from all these matrices.
`retained_runs.csv` is an unchanged copy of the reproducibility deposit's final
34-trajectory inventory and records the source names and velocity seeds.

## Reconstructing the V9 fingerprint PCA

The project-criteria PCA in main Figure 10 uses the updated 25-run antagonist
file and the nine-run agonist project file, as separate ensembles. For each
ensemble, retain residue labels whose maximum `score` across its trajectories
is at least 0.1. Pivot `score` into a trajectory-by-residue matrix, filling absent
contacts with zero. Centre each residue column without scaling its variance,
then apply PCA. This is the procedure in the manuscript figure-generation code.

The resulting project matrices have 25 trajectories by 32 residue features and
nine trajectories by 17 residue features. The 32 and 17 counts describe the
filtered PCA features, not all residue labels in the long-format contact files.
The first two variance percentages are 57.75758374 and 25.55001541 for the
antagonist ensemble, and 74.60795806 and 14.15660019 for the agonist ensemble.
The retained historical 24-run antagonist matrix gives 58.08306526 and
25.46707790 percent with the same procedure.

The staged source `fp_analysis/fingerprint_matrix.csv` is the earlier 24-run
matrix. It has not been relabelled as a 25-run source matrix: the updated input
was assembled from the balanced file and the extra-replicate file during figure
generation. The files added here expose the complete updated long-format data.

`validation_v9.json` records the unchanged original file hashes, exact field
preservation, ensemble sizes and recomputed variance percentages. It reports
data-integrity checks and PCA reconstruction; it does not assert a new run of
the trajectory-dependent regression tests.
