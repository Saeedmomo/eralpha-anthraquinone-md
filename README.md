# eralpha-anthraquinone-md

Analysis code and reference data for a molecular dynamics study of how
dihydroxyanthraquinones are recognised by the estrogen receptor α ligand-binding domain,
in both the antagonist and the agonist receptor conformation.

## Contents

| path | what it is |
|---|---|
| [`interaction_fingerprint/`](interaction_fingerprint) | the adapted per-residue ligand–protein interaction descriptor, its reference matrices, and tests that assert it reproduces them exactly |

## A note on what the descriptor reports

The per-type constants in the interaction descriptor are labelled kcal/mol in the source
it derives from. **They are fixed weights applied to contact counts, not computed
energies.** The `score` column is a weighted contact tally: it is not an interaction
energy, not a binding affinity, and a larger value does not indicate a better ligand. The
occupancy columns are fractions of frames and are the quantities to use for anything
physical. See [`interaction_fingerprint/README.md`](interaction_fingerprint/README.md).

## Requirements

MDAnalysis 2.11.0 or newer. GROMACS 2024 and later write `.tpr` files at tpx version 138,
which earlier MDAnalysis releases refuse to parse.

## Trajectories

The trajectories themselves are large and are not held in this repository. The tests take
their location from the `ERA_MD_ROOT` environment variable and skip cleanly without it.
