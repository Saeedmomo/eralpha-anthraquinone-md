# ERalpha anthraquinone molecular dynamics

Author and maintainer: [Said Moshawih](https://github.com/Saeedmomo).

This repository contains the interaction-fingerprint implementation and reference data for *Molecular dynamics guides anthraquinone drug design by resolving estrogen receptor α recognition geometry*, manuscript V11.

## V11 derived analyses and unchanged simulation inputs

The [V11 reproducibility release](https://github.com/Saeedmomo/eralpha-anthraquinone-md/releases/tag/v11-reproducibility-2026-11-03) provides a compact analysis supplement and a tagged copy of this repository. It adds the following records to the unchanged simulation-input archives in the [V9 release](https://github.com/Saeedmomo/eralpha-anthraquinone-md/releases/tag/v9-reproducibility-2026-11-01).

| Record | Contents |
|---|---|
| `data/v11_geometry/` | 11,022 frame-level observations from twenty-two two-hydroxyl antagonist trajectories; instantaneous hydroxyl span, cavity width and anchor distances; original-coordinate reconstruction and standard-library verification |
| `analysis/verify_v11_equations.py` | Derived ensemble regimes, allowance scan, engaged anchor distances and the 482-row weighted-score identity |
| `figures/v11_figure1/` | Updated Figure 1, original artwork, four-isomer first-run geometry, hydroxyl records, panel-F plotting source and integrity report |
| `reference/` | Balanced 24-run and updated 25-run antagonist fingerprints, nine-run agonist fingerprints and the retained-run inventory |

The instantaneous condition is `w(t) <= s(t) + 2d0`, with `d0 = 3.5 Å`. None of the 4,552 saved frames exceeding this bound showed simultaneous GLU353/HIS524 engagement. These are correlated trajectory frames, not independent replicates. The mean-width reach margin uses the separately prespecified allowance `a = 2.8 Å`. The two conditions are distinguished in the data and verification scripts.

The original four-isomer comparison uses three replicates per isomer. The expanded antagonist ensemble retains a fourth 2,6-DHAQ trajectory, giving 25 antagonist runs. Figure 1F uses first-replicate spans and within-trajectory standard deviations; its 8.62 and 7.22 Å guides illustrate equations (3) and (4) at the rounded first-2,6 mean width of 14.22 Å. They are reference-width examples, not universal cutoffs.

The unchanged V9 archives remain the input source:

| Archive | Contents |
|---|---|
| `eralpha_anthraquinone_reproducibility_V9.zip` | Inputs and topologies for 34 retained conventional runs, 45 umbrella windows and twelve AWH walkers; seed, extension and analysis records; pull series; fingerprints and locally dated predictions |
| `legacy53_processed_reproducibility_V9.zip` | Processed historical screening data and reconstruction scripts |
| `V9_supporting_figures.zip` | Earlier supporting-series figure corrections and their source tables and validation records |

Umbrella records cover 708 ns of nominal sampling and AWH records cover 87.504 ns. Coordinate trajectories are excluded. The V11 frame table permits verification of the published bound without coordinates; regenerating that table requires the saved coordinate and contact records. Repeating simulations requires the V9 inputs and the stated software. Local prediction timestamps document project provenance, not external registration.

The historical Schrödinger series comprises eleven 100 ns and forty-two 50 ns trajectories. Comparative analyses use the common first 50 ns. All 53 records have energy data; 51 have stored structural/contact observations. The hydroxyl-count–engagement association is Pearson r = 0.673 for n = 51, r = 0.671 after preparation exclusions (n = 49), and r = 0.723 in the energy-stable subset (n = 18). Structural fields for 2OUZ and 3ERT remain unavailable. Marginal contacts to both anchors are distinct from simultaneous bridging.

## Fingerprint definitions

The implementation reads `prod.tpr` and `prod.xtc`, reports weighted contact scores and per-type frame occupancies, canonicalises AMBER residue names, applies minimum-image distances and retains receptor residue identity. Use `--criteria project` for the V11 fingerprint calculations.

| Term | Project implementation |
|---|---|
| Hydrogen bond | Heavy-atom distance below 3.5 Å and bonded D–H···A angle above 135°; candidate atoms are N, O or S |
| Hydrophobic | Ligand carbon within 4.0 Å of a side-chain carbon or sulfur in a listed hydrophobic residue |
| Ionic-labelled proximity | Ligand atom within 3.7 Å of ARG, LYS, ASP or GLU; no formal-charge assignment is performed |
| Cation–π-labelled proximity | Ligand nitrogen within 4.5 Å of a protein aromatic-ring centroid |
| Polar proximity, stored as `occ_water_bridge` | Direct ligand–protein N/O/S distance below 2.8 Å; no water molecule is required |

The last column name is retained for compatibility with the archived matrices. Explicit single-water bridges in the manuscript were calculated in a separate trajectory analysis. The fingerprint hydrogen-bond atom selection also differs from the independent cpptraj N/O-only inventory; the SI documents the distinction.

Weights are 5.0, 1.5, 4.0, 3.5 and 2.0, respectively, in the table order above. Scores are contact tallies in arbitrary units. They are not binding energies or affinities. A frame can contribute to several categories, so occupancies should be compared by their stated definitions.

`--criteria github` retains the earlier implementation thresholds, including a 2.5 Å heavy-atom hydrogen-bond cutoff. This remains the API default for compatibility. Select `--criteria project` explicitly to reproduce the manuscript fingerprints.

## Installation

Python 3.11 or later and MDAnalysis 2.11.0 or later are required for the GROMACS 2026 TPR inputs used here. MDAnalysis 2.11 added GROMACS 2026 support; see the [official release notes](https://www.mdanalysis.org/2026/10/21/release-2.11.0/).

```bash
python -m venv venv
source venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows PowerShell, activate the environment with `venv\Scripts\Activate.ps1`. Alternatively, create the supplied environment with `mamba env create -f environment.yml` and activate `era-fingerprint`.

## Run the fingerprint analysis

```bash
python interaction_fingerprint.py \
    --root /path/to/md_s4 \
    --runs dihydroxy23 alizarin_r3 estradiol_r1 \
    --criteria project \
    --out-prefix fp_project
```

The script writes `<prefix>_by_type.csv` and `<prefix>_by_residue.csv`. Each named run requires its production TPR and XTC. `--stride` controls frame sampling, and `--offset` maps topology residue numbers to ERalpha numbering; the project default is 305.

For the agonist systems containing the TIF2 peptide, restrict receptor numbering to topology residues 1–245:

```bash
python interaction_fingerprint.py \
    --root /path/to/md_s4 \
    --runs agopep_26dhaq_r1 \
    --criteria project \
    --receptor-resid-range 1 245 \
    --out-prefix fp_agonist
```

Residues outside that interval retain their topology numbers, have `chain=other`, and receive no ERalpha residue number.

## Reference data and verification

The original antagonist reference files contain the balanced 24-run ensemble used for main Table 2, SI Table S8 and SI Figure S10. The updated files add the fourth retained 2,6-isomer replicate for the 25-run antagonist analysis. Agonist references contain nine runs. See [the reference inventory](reference/README.md) for file roles and the PCA selection.

```bash
export ERA_MD_ROOT=/path/to/md_s4
pytest -q
```

The existing trajectory regression tests compare three antagonist runs under both criteria sets with their archived values. Separate checks cover receptor numbering and configurable residue offsets. Trajectory-dependent tests skip when `ERA_MD_ROOT` is unset; such a run does not verify trajectory reproduction.

## Verify the V11 derived analyses

The following checks require only Python 3's standard library and the files in this repository:

```bash
python data/v11_geometry/verify_v11_frame_bound.py
python analysis/verify_v11_equations.py
```

The first verifies the instantaneous framewise bound. The second reconstructs the mean-width regimes, allowance sensitivity, anchor-distance calibration and weighted fingerprint scores. Coordinate extraction and figure regeneration have separate scientific-package requirements documented beside their sources. These checks do not claim to rerun the simulations or the trajectory-dependent fingerprint tests.

## Citation and provenance

Moshawih, S. (2026). *ERalpha anthraquinone molecular dynamics inputs and analysis data V11* [Data set]. GitHub, release v11-reproducibility-2026-11-03. [Versioned release](https://github.com/Saeedmomo/eralpha-anthraquinone-md/releases/tag/v11-reproducibility-2026-11-03).

Repository citation metadata are provided in [CITATION.cff](CITATION.cff). The contact scheme derives from the earlier [ligand–protein interaction implementation](https://github.com/Saeedmomo/Binding-Affinity-Prediction-Using-Molecular-Dynamics). The archive license applies to newly generated numerical data and records; software and force-field materials retain their documented original terms.
