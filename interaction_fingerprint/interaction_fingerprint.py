"""Per-residue ligand-protein interaction fingerprints for the Stage A trajectories.

Derived from ligand_protein_interaction_forces.py in
Saeedmomo/Binding-Affinity-Prediction-Using-Molecular-Dynamics, which classifies each
ligand-residue contact as hydrogen bond, cation-pi, hydrophobic, ionic or water-bridged,
weights it by a per-type constant and averages over frames. The original reads one
structure file and aggregates by residue NAME. Four changes were needed here, all
recorded so the difference from the published descriptor is explicit:

  1. GROMACS input. The original takes a single PDB. These runs are prod.tpr + prod.xtc,
     34 941 atoms (35 328 with the coactivator peptide), 501 frames at 100 ps. Reading a
     .tpr written by GROMACS 2024 or later requires MDAnalysis 2.11 or newer, because
     earlier releases do not support tpx file version 138; this is checked at startup.

  2. AMBER residue names. The original matches against the 20 standard names, but
     AMBER writes protonation states: this system has 6 HID and 7 HIE and no HIS at
     all, so every histidine - including HIS524, the residue the whole project turns
     on - would have scored exactly zero. Names are canonicalised before matching.

  3. Periodic images. Distances and angles are computed under the minimum-image
     convention rather than on raw coordinates, so a ligand that has drifted across a
     box face is still measured against the protein and not against empty space. This
     avoids unwrapping 34 941 atoms per frame.

  4. Residue identity is kept. The published descriptor sums every ALA in the protein
     into one number. For ligand design against one pocket that discards the only thing
     worth knowing, so the same score is also emitted per residue, in ERalpha numbering
     (prmtop resid + 305, the project's OFF).

--criteria github reproduces the repo's thresholds exactly, subject to 1-3 above.
--criteria project substitutes this project's established definitions: a hydrogen bond
is a heavy-atom donor-acceptor pair within 3.5 A with a real D-H...A angle above 135
degrees, measured with the actual bonded hydrogen rather than the original's
nearest-neighbour stand-in, and hydrophobic contact is a ligand carbon within 4.0 A of
a protein side-chain carbon or sulfur.

NOTE ON UNITS. The per-type constants are labelled kcal/mol in the source. They are
fixed weights applied to contact counts, not computed energies, and the column totals
are weighted contact tallies. They are not interaction energies and must never be
reported as such.
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

MIN_MDANALYSIS = (2, 11, 0)


def _require_mdanalysis(mda):
    """GROMACS 2026 writes tpx 138, which MDAnalysis below 2.11 refuses to parse. The
    failure is a bare parser traceback, so say what to do about it instead."""
    got = tuple(int(x) for x in mda.__version__.split('.')[:3]
                if x.isdigit()) or (0,)
    if got < MIN_MDANALYSIS:
        raise SystemExit(
            'MDAnalysis %s is installed; this script needs >= %s to read the .tpr files '
            'written by GROMACS 2024 and later (tpx 138). Install it with\n'
            '    pip install "MDAnalysis>=%s"'
            % (mda.__version__, '.'.join(map(str, MIN_MDANALYSIS)),
               '.'.join(map(str, MIN_MDANALYSIS))))


AMINO = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
         'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
HYDROPHOBIC = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR'}
CHARGED = {'ARG', 'LYS', 'ASP', 'GLU'}
AROMATIC = {'PHE', 'TYR', 'TRP', 'HIS'}

# AMBER protonation and disulfide variants map onto the standard twenty.
CANON = {'HID': 'HIS', 'HIE': 'HIS', 'HIP': 'HIS', 'HISE': 'HIS', 'HISD': 'HIS',
         'CYX': 'CYS', 'CYM': 'CYS', 'ASH': 'ASP', 'GLH': 'GLU', 'LYN': 'LYS',
         'ARN': 'ARG', 'TYM': 'TYR'}

CUTOFF = {'hbond': 2.5, 'cation_pi': 4.5, 'hydrophobic': 3.6,
          'ionic': 3.7, 'water_bridge': 2.8}
ENERGY = {'hbond': 5.0, 'cation_pi': 3.5, 'hydrophobic': 1.5,
          'ionic': 4.0, 'water_bridge': 2.0}
ANGLE = {'hbond_donor': 120.0, 'hbond_acceptor': 90.0}

PROJECT_CUTOFF = dict(CUTOFF, hbond=3.5, hydrophobic=4.0)
DONOR_ACCEPTOR = ('N', 'O', 'S')
OFF = 305
RING = {'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
        'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2']}


def mic(delta, box):
    """Minimum-image displacement for an orthorhombic box."""
    return delta - box * np.round(delta / box)


def mic_angle(a, b, c, box):
    """Angle a-b-c in degrees with b at the vertex, under minimum image."""
    v1, v2 = mic(a - b, box), mic(c - b, box)
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 < 1e-8 or n2 < 1e-8:
        return float('nan')
    return float(np.degrees(np.arccos(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))))


def analyse(tpr, xtc, criteria='github', stride=1, ligand_sel='resname LIG',
            offset=OFF, receptor_range=None):
    import MDAnalysis as mda
    from MDAnalysis.analysis.distances import distance_array

    _require_mdanalysis(mda)
    u = mda.Universe(tpr, xtc)
    ligand = u.select_atoms(ligand_sel)
    protein = u.select_atoms('protein')
    if ligand.n_atoms == 0:
        raise ValueError(f'{ligand_sel!r} matched no atoms')

    cutoff = PROJECT_CUTOFF if criteria == 'project' else CUTOFF

    lig_el = np.array([n[0].upper() for n in ligand.names])
    prot_el = np.array([n[0].upper() for n in protein.names])
    prot_raw = np.array([a.residue.resname for a in protein])
    prot_resname = np.array([CANON.get(r, r) for r in prot_raw])
    prot_resid = np.array([a.residue.resid for a in protein])
    # side chain = everything that is not backbone N, CA, C, O
    prot_side = np.array([n not in ('N', 'CA', 'C', 'O') for n in protein.names])

    # Real donor hydrogens, for the project criterion. The tpr carries bonds, so each
    # polar hydrogen can be attached to the heavy atom it actually belongs to.
    lig_local = {a.index: k for k, a in enumerate(ligand)}
    prot_local = {a.index: k for k, a in enumerate(protein)}
    lig_h, prot_h = {}, {}
    for i, at in enumerate(ligand):
        if at.name[0].upper() in DONOR_ACCEPTOR:
            lig_h[i] = [lig_local[b.index] for b in at.bonded_atoms
                        if b.name[0].upper() == 'H' and b.index in lig_local]
    for j, at in enumerate(protein):
        if at.name[0].upper() in DONOR_ACCEPTOR:
            prot_h[j] = [prot_local[b.index] for b in at.bonded_atoms
                         if b.name[0].upper() == 'H' and b.index in prot_local]

    # aromatic ring atom indices, precomputed once
    rings = []
    for res in protein.residues:
        rn = CANON.get(res.resname, res.resname)
        if rn in RING:
            idx = [k for k, a in enumerate(protein) if a.residue.resindex == res.resindex
                   and a.name in RING[rn]]
            if len(idx) >= 3:
                rings.append((rn, res.resid, np.array(idx)))

    by_type = defaultdict(float)
    by_res = defaultdict(float)
    counts = defaultdict(int)
    per_type_res = defaultdict(lambda: defaultdict(float))
    n_frames = 0

    for ts in u.trajectory[::stride]:
        n_frames += 1
        box = ts.dimensions[:3]
        lp, pp = ligand.positions, protein.positions
        d = distance_array(lp, pp, box=ts.dimensions)
        frame = defaultdict(set)              # kind -> {(resname, resid)}

        # ---- hydrogen bonds
        li, pi = np.where(d < cutoff['hbond'])
        for i, j in zip(li, pi):
            if lig_el[i] not in DONOR_ACCEPTOR or prot_el[j] not in DONOR_ACCEPTOR:
                continue
            ok = False
            if criteria == 'project':
                # genuine D-H...A on either direction of donation: the ligand heavy
                # atom donating to the protein one, or the reverse
                for hs, hpos, dn, ac in ((lig_h.get(i, []), lp, lp[i], pp[j]),
                                         (prot_h.get(j, []), pp, pp[j], lp[i])):
                    for h in hs:
                        if mic_angle(dn, hpos[h], ac, box) > 135.0:
                            ok = True
                            break
                    if ok:
                        break
            else:
                # the repo's geometry: the nearest heavy neighbour stands in for the
                # covalent partner that defines the bond direction
                dl = np.linalg.norm(mic(lp - lp[i], box), axis=1)
                dl[i] = np.inf
                k = int(np.argmin(dl))
                dp = np.linalg.norm(mic(pp - pp[j], box), axis=1)
                dp[j] = np.inf
                m = int(np.argmin(dp))
                donor = mic_angle(lp[k], lp[i], pp[j], box)
                acceptor = mic_angle(lp[i], pp[j], pp[m], box)
                ok = (np.isfinite(donor) and np.isfinite(acceptor)
                      and donor >= ANGLE['hbond_donor']
                      and acceptor >= ANGLE['hbond_acceptor'])
            if ok:
                frame['hbond'].add((prot_resname[j], prot_resid[j]))

        # ---- hydrophobic
        li, pi = np.where(d < cutoff['hydrophobic'])
        for i, j in zip(li, pi):
            if prot_resname[j] not in HYDROPHOBIC:
                continue
            if criteria == 'project' and not (
                    lig_el[i] == 'C' and prot_el[j] in ('C', 'S') and prot_side[j]):
                continue
            frame['hydrophobic'].add((prot_resname[j], prot_resid[j]))

        # ---- ionic
        li, pi = np.where(d < cutoff['ionic'])
        for i, j in zip(li, pi):
            if prot_resname[j] in CHARGED:
                frame['ionic'].add((prot_resname[j], prot_resid[j]))

        # ---- water bridge, as the repo defines it (no water is actually involved;
        #      it is a second, wider polar-contact term). Kept for comparability.
        li, pi = np.where(d < cutoff['water_bridge'])
        for i, j in zip(li, pi):
            if lig_el[i] in DONOR_ACCEPTOR and prot_el[j] in DONOR_ACCEPTOR:
                frame['water_bridge'].add((prot_resname[j], prot_resid[j]))

        # ---- cation-pi
        cations = np.where(lig_el == 'N')[0]
        if len(cations):
            for rn, rid, idx in rings:
                cen = pp[idx].mean(0)
                if np.min(np.linalg.norm(mic(lp[cations] - cen, box), axis=1)) \
                        < cutoff['cation_pi']:
                    frame['cation_pi'].add((rn, rid))

        for kind, keys in frame.items():
            for rn, rid in keys:
                by_type[rn] += ENERGY[kind]
                by_res[(rn, rid)] += ENERGY[kind]
                per_type_res[kind][(rn, rid)] += 1.0
            counts[kind] += len(keys)

    if n_frames == 0:
        raise ValueError('no frames')

    types = {r: by_type.get(r, 0.0) / n_frames for r in AMINO}
    types['_n_frames'] = n_frames
    for k in CUTOFF:
        types[f'_n_{k}'] = counts.get(k, 0)

    res_rows = []
    for (rn, rid), v in sorted(by_res.items(), key=lambda t: -t[1]):
        # Only residues inside the receptor chain carry the renumbering. Anything else -
        # a coactivator peptide, a second protomer - keeps its raw resid and is flagged,
        # because offsetting it would invent a receptor residue number.
        inside = receptor_range is None or receptor_range[0] <= rid <= receptor_range[1]
        row = {'resname': rn, 'resid': rid,
               'era': (rid + offset) if inside else '',
               'chain': 'receptor' if inside else 'other',
               'label': f'{rn}{rid + offset}' if inside else f'{rn}{rid}@other',
               'score': v / n_frames}
        for k in CUTOFF:
            row[f'occ_{k}'] = per_type_res[k].get((rn, rid), 0.0) / n_frames
        res_rows.append(row)
    return types, pd.DataFrame(res_rows), n_frames


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', default=os.path.expanduser('~/md_s4'))
    ap.add_argument('--runs', nargs='+', required=True)
    ap.add_argument('--criteria', choices=('github', 'project'), default='github')
    ap.add_argument('--offset', type=int, default=OFF,
                    help='added to the topology resid to give the receptor numbering '
                         '(default %(default)s, correct for this project\'s ERalpha '
                         'systems; set it for your own)')
    ap.add_argument('--receptor-resid-range', type=int, nargs=2, metavar=('LO', 'HI'),
                    default=None,
                    help='topology resid range of the receptor chain. Residues outside '
                         'it keep their raw numbering and are marked chain=other, which '
                         'is what you want when the system carries a second chain such '
                         'as a coactivator peptide. Default: treat everything as '
                         'receptor.')
    ap.add_argument('--stride', type=int, default=1)
    ap.add_argument('--out-prefix', required=True)
    a = ap.parse_args()

    rows, per_res = [], []
    for run in a.runs:
        tpr = os.path.join(a.root, run, 'prod.tpr')
        xtc = os.path.join(a.root, run, 'prod.xtc')
        if not (os.path.exists(tpr) and os.path.exists(xtc)):
            print(f'  [skip] {run}: missing tpr or xtc', file=sys.stderr)
            continue
        try:
            types, df, nf = analyse(tpr, xtc, a.criteria, a.stride,
                                    offset=a.offset,
                                    receptor_range=tuple(a.receptor_resid_range)
                                    if a.receptor_resid_range else None)
            rows.append({'run': run, **types})
            df.insert(0, 'run', run)
            per_res.append(df)
            top = ', '.join(df.head(3)['label'] + '=' + df.head(3)['score'].round(1).astype(str))
            print(f'  [ok] {run:20s} {nf:4d} frames  hb {types["_n_hbond"]:5d} '
                  f'phob {types["_n_hydrophobic"]:6d} ionic {types["_n_ionic"]:5d}  | {top}')
        except Exception as e:
            print(f'  [FAIL] {run}: {e}', file=sys.stderr)

    if not rows:
        return 1
    pd.DataFrame(rows).to_csv(f'{a.out_prefix}_by_type.csv', index=False)
    pd.concat(per_res, ignore_index=True).to_csv(f'{a.out_prefix}_by_residue.csv',
                                                 index=False)
    print(f'\nwrote {a.out_prefix}_by_type.csv  ({len(rows)} runs x 20 amino acids)')
    print(f'wrote {a.out_prefix}_by_residue.csv')
    return 0


if __name__ == '__main__':
    sys.exit(main())
