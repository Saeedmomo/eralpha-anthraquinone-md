# ERalpha anthraquinone molecular dynamics

Author and maintainer: [Said Moshawih](https://github.com/Saeedmomo).

This repository contains the interaction-fingerprint implementation and reference data for *Hydroxyl availability and spacing govern anthraquinone recognition by estrogen receptor alpha*, manuscript and Supporting Information V9.

## Simulation inputs and analysis data

The [V9 reproducibility release](https://github.com/Saeedmomo/eralpha-anthraquinone-md/releases/tag/v9-reproducibility-2026-11-01) contains the following archives. SHA-256 files and validation records accompany the downloads.

| Archive | Contents |
|---|---|
| `eralpha_anthraquinone_reproducibility_V9.zip` | Inputs and topologies for 34 retained conventional runs, 45 umbrella windows and 12 AWH walkers; seed and extension records; pull series; WHAM records; analysis scripts; fingerprints; locally dated prediction files |
| `legacy53_processed_reproducibility_V9.zip` | Processed historical screening data and reconstruction scripts for the corrected supporting-series analysis |
| `V9_supporting_figures.zip` | Corrected main Figure 11 and SI Figure S4, their plotting scripts, source tables and validation records |

Umbrella records cover 708 ns of nominal sampling, and AWH records extend to 87.504 ns. Coordinate trajectories are excluded. The archives provide simulation inputs and processed analysis records; repeating contact calculations from coordinates requires the original trajectories. Local prediction files preserve project records and timestamps; no external registration is claimed.

The historical series comprises 53 energy records and 51 complexes with stored structural/contact data. The primary hydroxyl-count–engagement association is Pearson r = 0.673 for n = 51. Preparation exclusions give n = 49 and r = 0.671; the energy-stable subset gives n = 18 and r = 0.723. Missing structural observations for 2OUZ and 3ERT are marked unavailable. Contacts to both anchors across snapshots are distinguished from simultaneous bridging.

## Fingerprint definitions

The implementation reads `prod.tpr` and `prod.xtc`, reports weighted contact scores and per-type frame occupancies, canonicalises AMBER residue names, applies minimum-image distances and retains receptor residue identity. Use `--criteria project` for the V9 fingerprint calculations.

| Term | Project implementation |
|---|---|
| Hydrogen bond | Heavy-atom distance below 3.5 Å and bonded D–H···A angle above 135°; candidate atoms are N, O or S |
| Hydrophobic | Ligand carbon within 4.0 Å of a side-chain carbon or sulfur in a listed hydrophobic residue |
| Ionic-labelled proximity | Ligand atom within 3.7 Å of ARG, LYS, ASP or GLU; no formal-charge assignment is performed |
| Cation–π-labelled proximity | Ligand nitrogen within 4.5 Å of a protein aromatic-ring centroid |
| Polar proximity, stored as `occ_water_bridge` | Direct ligand–protein N/O/S distance below 2.8 Å; no water molecule is required |

The last column name is retained for compatibility with the archived matrices. Explicit single-water bridges in V9 were calculated in a separate trajectory analysis. The fingerprint hydrogen-bond atom selection also differs from the independent cpptraj N/O-only inventory; the SI documents the distinction.

Weights are 5.0, 1.5, 4.0, 3.5 and 2.0, respectively, in the table order above. Scores are contact tallies in arbitrary units. They are not binding energies or affinities. A frame can contribute to several categories, so occupancies should be compared by their stated definitions.

`--criteria github` retains the earlier implementation thresholds, including a 2.5 Å heavy-atom hydrogen-bond cutoff. This remains the API default for compatibility. Select `--criteria project` explicitly to reproduce V9.

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

The original antagonist reference files contain the balanced 24-run ensemble used for main Table 7, SI Table S8 and SI Figure S5. The updated files add the fourth retained 2,6-isomer replicate for the 25-run antagonist analysis. Agonist references contain nine runs. See [the reference inventory](reference/README.md) for file roles and the PCA selection.

```bash
export ERA_MD_ROOT=/path/to/md_s4
pytest -q
```

The existing trajectory regression tests compare three antagonist runs under both criteria sets with their archived values. Separate checks cover receptor numbering and configurable residue offsets. Trajectory-dependent tests skip when `ERA_MD_ROOT` is unset; such a run does not verify trajectory reproduction.

## Citation and provenance

Moshawih, S. (2026). *ERalpha anthraquinone recognition simulation inputs and analysis data V9* [Data set]. GitHub, release v9-reproducibility-2026-11-01. [Versioned release](https://github.com/Saeedmomo/eralpha-anthraquinone-md/releases/tag/v9-reproducibility-2026-11-01).

Repository citation metadata are provided in [CITATION.cff](CITATION.cff). The contact scheme derives from the earlier [ligand–protein interaction implementation](https://github.com/Saeedmomo/Binding-Affinity-Prediction-Using-Molecular-Dynamics). The archive license applies to newly generated numerical data and records; software and force-field materials retain their documented original terms.
