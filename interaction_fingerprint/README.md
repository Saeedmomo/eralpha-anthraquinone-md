# Per-residue ligand–protein interaction fingerprints

Adapted from [`ligand_protein_interaction_forces.py`](https://github.com/Saeedmomo/Binding-Affinity-Prediction-Using-Molecular-Dynamics),
which classifies every ligand–residue contact as a hydrogen bond, cation–π, hydrophobic,
ionic or water-mediated contact, weights each by a per-type constant and averages over
frames. This version reads GROMACS trajectories, keeps residue identity, and runs on
both the antagonist- and agonist-state estrogen receptor α systems used in this study.

## What it is not

**The per-type constants are labelled kcal/mol in the original source. They are fixed
weights applied to contact counts, not computed energies.** The `score` column is a
weighted contact tally. It is not an interaction energy, not a binding affinity, and a
larger value does not mean a better ligand. One physical contact can also enter the sum
more than once under different thresholds, and the ionic term applies no charge test, so
absolute magnitudes carry no meaning. Compare within one residue column, never across
columns and never against an affinity. The occupancy columns are fractions of frames and
are the quantities to use for anything physical.

## Install

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```

**MDAnalysis 2.11.0 or newer is required.** GROMACS 2024 and later write `.tpr` files at
tpx version 138, and earlier MDAnalysis releases refuse to parse them with
`Your tpx version is 138, which this parser does not support`. The script checks the
version at startup and tells you what to install.

## Use

```bash
python interaction_fingerprint.py \
    --root /path/to/runs \
    --runs system_r1 system_r2 system_r3 \
    --criteria project \
    --out-prefix fp_project
```

Each run directory must contain `prod.tpr` and `prod.xtc`. Two CSVs are written:
`<prefix>_by_type.csv`, one row per run aggregated over the twenty amino-acid types, and
`<prefix>_by_residue.csv`, one row per run and residue.

### Options that matter

| option | meaning |
|---|---|
| `--criteria project` | this work's definitions: hydrogen bond 3.5 Å with a real donor–H···acceptor angle above 135°, hydrophobic contact 4.0 Å |
| `--criteria github` | the original repository's thresholds: hydrogen bond 2.5 Å heavy-atom, hydrophobic 3.6 Å. On these systems this returns **zero** hydrogen bonds in every trajectory, because 2.5 Å between heavy atoms is shorter than a hydrogen bond |
| `--offset N` | added to the topology resid to give receptor numbering. Default 305, which is correct for the ERα systems here. **Set it for your own system.** |
| `--receptor-resid-range LO HI` | topology resid range of the receptor chain. Residues outside it keep their raw numbering and are marked `chain=other`. Use this whenever the system carries a second chain |
| `--stride N` | analyse every Nth frame |

### A system with more than one chain

The agonist-state systems here carry a TIF2 coactivator peptide at topology resid
246–256. Without `--receptor-resid-range` every residue is offset, so the peptide would
be reported as ERα 551–561 — residue numbers that do not exist:

```bash
python interaction_fingerprint.py --root /path/to/runs \
    --runs agopep_26dhaq_r1 --criteria project \
    --receptor-resid-range 1 245 --out-prefix fp_agonist
```

Rows outside the range then carry `chain=other`, a blank `era`, and a label of the form
`LEU248@other`.

## Four changes from the published descriptor

Each changes the numbers rather than their presentation, so each is recorded.

1. **GROMACS input.** The original reads a single PDB. These are `prod.tpr` + `prod.xtc`,
   34 941 atoms (35 328 with the peptide), 501 frames at 100 ps.
2. **AMBER residue names.** The original matches the twenty standard names. These systems
   are AMBER-built and contain HID and HIE and no residue named HIS, so under the
   unmodified matcher every histidine in the pocket — HIS524 included, the residue this
   project turns on — scored exactly zero. Names are canonicalised before matching.
3. **Periodic images.** Distances and angles are evaluated under the minimum-image
   convention, so a ligand that has drifted across a box face is still measured against
   the protein and not against empty space.
4. **Residue identity is kept.** The published descriptor sums every ALA in the protein
   into one number. The same score is also emitted per residue, in receptor numbering.

Under `--criteria project` the hydrogen bond is evaluated against the **actual bonded
hydrogen**, taken from the topology, rather than the original's nearest-neighbour
stand-in. This is why the `.tpr` is read rather than a `.gro`: a `.gro` carries no bonds.

## Verifying an installation

```bash
export ERA_MD_ROOT=/path/to/runs
pytest -q
```

The tests assert that the descriptor reproduces the reference matrices in `reference/`
with a **maximum absolute difference of zero** across the score and all five occupancy
columns, that a second chain is never given receptor numbering, and that `--offset`
shifts the numbering without changing any score. Without `ERA_MD_ROOT` the
trajectory-dependent tests skip and the rest still run.

The reference CSVs were produced from the published antagonist-state trajectories. They
are the contract: if a change to this script, to MDAnalysis, or to the toolchain alters
them, the tests fail rather than the numbers changing quietly.

## Reproducing the published analyses

Antagonist state, 24 trajectories:

```bash
python interaction_fingerprint.py --root $ERA_MD_ROOT --criteria project \
  --runs alizarin alizarin_r2 alizarin_r3 alizarin_anion_r1 alizarin_anion_r2 \
         alizarin_anion_r3 anthraflavic anthraflavic_r2 anthraflavic_r3 \
         danthron danthron_r2 danthron_r3 dihydroxy23 dihydroxy23_r2 dihydroxy23_r3 \
         dihydroxy27_r1 dihydroxy27_r2 dihydroxy27_r3 estradiol_r1 estradiol_r2 \
         estradiol_r3 nativeE4D_r1 nativeE4D_r2 nativeE4D_r3 \
  --out-prefix fp_project
```

Agonist state, 9 trajectories, with the peptide excluded from renumbering:

```bash
python interaction_fingerprint.py --root $ERA_MD_ROOT --criteria project \
  --receptor-resid-range 1 245 \
  --runs agopep_26dhaq_r1 agopep_26dhaq_r2 agopep_26dhaq_r3 \
         agopep_27dhaq_r1 agopep_27dhaq_r2 agopep_27dhaq_r3 \
         agopep_estradiol_r1 agopep_estradiol_r2 agopep_estradiol_r3 \
  --out-prefix fp_agonist
```

Note that the antagonist 2,6-dihydroxyanthraquinone ensemble reported in the manuscript
is four replicates, `anthraflavic`, `_r2`, `_r3` and `_r5`; `anthraflavic_r4` is
quarantined on a trajectory-quality criterion and is not part of any analysis. The
reference matrix above predates that decision and contains the first three only, which
is why it is a reference for reproducibility rather than the manuscript ensemble.
