# Per-residue ligand–protein interaction fingerprints

A per-residue interaction descriptor for GROMACS molecular dynamics trajectories. It
classifies every ligand–residue contact as a hydrogen bond, a cation–π contact, a
hydrophobic contact, an ionic contact or a water-mediated contact, and returns both an
aggregated twenty-amino-acid vector and a residue-resolved profile in the receptor's own
numbering.

It was developed for a study of how dihydroxyanthraquinones are recognised by the
estrogen receptor α ligand-binding domain, across both the antagonist and the agonist
receptor conformation, and it runs on either without modification.

## What makes it different from a contact count

**Real bonded hydrogens.** Under `--criteria project` a hydrogen bond requires a
heavy-atom donor–acceptor pair within 3.5 Å *and* a genuine D–H···A angle above 135°,
measured against the hydrogen actually bonded to the donor in the topology rather than a
nearest-neighbour stand-in. This is why the `.tpr` is read rather than a `.gro`: a `.gro`
carries no bonds.

**AMBER protonation states are handled.** AMBER-built systems contain `HID` and `HIE` and
no residue named `HIS`. A matcher written against the twenty standard names scores every
histidine as zero — including the catalytically relevant ones. Residue names are
canonicalised before matching.

**Periodic images are respected.** Distances and angles are evaluated under the
minimum-image convention, so a ligand that has drifted across a box face is still
measured against the protein rather than against empty space, without unwrapping tens of
thousands of atoms per frame.

**Residue identity is preserved.** Aggregating every alanine in a protein into one number
discards the only thing worth knowing when working against a single binding site, so the
score is also emitted per residue.

**Multi-chain systems are not silently misnumbered.** See `--receptor-resid-range` below.

## What the numbers are, and are not

The per-type constants are weights applied to contact counts. **The `score` column is a
weighted contact tally, not an energy.** It is not an interaction energy, not a binding
affinity, and a larger value does not mean a better ligand. One physical contact can enter
the sum more than once under different thresholds, and the ionic term applies no charge
test, so absolute magnitudes carry no meaning. Compare within one residue column, never
across columns and never against an affinity.

The `occ_*` columns are fractions of frames. **Those are the quantities to use for
anything physical.**

## Install

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```

or with conda/mamba:

```bash
mamba env create -f environment.yml && mamba activate era-fingerprint
```

**MDAnalysis 2.11.0 or newer is required.** GROMACS 2024 and later write `.tpr` files at
tpx version 138, and earlier MDAnalysis releases refuse to parse them with
`Your tpx version is 138, which this parser does not support`. The script checks at
startup and tells you what to install rather than failing with a parser traceback.

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

| option | meaning |
|---|---|
| `--criteria project` | hydrogen bond 3.5 Å with a real donor–H···acceptor angle above 135°, hydrophobic contact 4.0 Å |
| `--criteria github` | 2.5 Å heavy-atom hydrogen bond, hydrophobic 3.6 Å. Note that a 2.5 Å heavy-atom cut-off is shorter than a hydrogen bond and returns **zero** hydrogen bonds on these systems |
| `--offset N` | added to the topology resid to give receptor numbering (default 305, which suits the ERα systems this was written for — **set it for your own**) |
| `--receptor-resid-range LO HI` | topology resid range of the receptor chain; residues outside it keep their raw numbering and are marked `chain=other` |
| `--stride N` | analyse every Nth frame |

### Systems with more than one chain

If a system carries a second chain — a coactivator peptide, a second protomer — then
offsetting every residue maps that chain onto receptor numbers that do not exist. In the
agonist-state systems here the TIF2 peptide occupies topology resid 246–256 and would be
reported as ERα 551–561. Give the receptor's range instead:

```bash
python interaction_fingerprint.py --root /path/to/runs \
    --runs agopep_26dhaq_r1 --criteria project \
    --receptor-resid-range 1 245 --out-prefix fp_agonist
```

Rows outside the range then carry `chain=other`, a blank `era` column and a label of the
form `LEU248@other`, so a contact with a second chain is visible as such.

## Verifying an installation

```bash
export ERA_MD_ROOT=/path/to/runs
pytest -q
```

The tests assert that the descriptor reproduces the reference matrices in `reference/`
with a **maximum absolute difference of zero** across the score and all five occupancy
columns, that a second chain is never given receptor numbering, and that `--offset` shifts
the numbering without changing any score. Without `ERA_MD_ROOT` the trajectory-dependent
tests skip and the rest still run, so `pytest` is useful on a fresh clone with no data.

The reference CSVs are the contract: if a change to this script, to MDAnalysis or to the
toolchain alters the numbers, the tests fail rather than the results changing quietly.

## Trajectories

The trajectories are large and are not held in this repository. `--root` and the
`ERA_MD_ROOT` test variable both point at the directory containing the run folders.

## Provenance

The contact-classification scheme derives from
[`ligand_protein_interaction_forces.py`](https://github.com/Saeedmomo/Binding-Affinity-Prediction-Using-Molecular-Dynamics),
from earlier work in this line. Everything described under *What makes it different from a
contact count* above is specific to this implementation.
