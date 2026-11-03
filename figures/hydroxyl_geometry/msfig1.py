"""Figure 1 - the four regioisomers, their measured hydroxyl state, and the two descriptors."""
import io

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from msfig_common import (C, DOUBLE, NAME, ROW, SHORT, SMILES, SYS, panel, rd,
                          save, style)

FREE = '#009E73'      # hydroxyl free of the peri chelate
CHEL = '#9E9E9E'      # hydroxyl consumed by the intramolecular chelate
QUIN = '#D55E00'      # quinone carbonyl oxygen


def rgb(h):
    return tuple(int(h[i:i + 2], 16) / 255 for i in (1, 3, 5))


def peri_hydroxyls(mol):
    """(peri, free) hydroxyl oxygen indices, decided from connectivity.

    A hydroxyl is peri when the aromatic carbon carrying it is adjacent to a
    ring-fusion carbon that is itself bonded to a quinone carbonyl carbon. That is
    exactly the arrangement that lets O-H fold back onto C=O, and it reproduces the
    classification measured from the trajectories in s4_an_hydroxyl.csv.
    """
    carbonylC = {a.GetIdx() for a in mol.GetAtoms()
                 if a.GetSymbol() == 'C' and
                 any(b.GetBondTypeAsDouble() == 2 and
                     b.GetOtherAtom(a).GetSymbol() == 'O' for b in a.GetBonds())}
    peri, free = [], []
    for a in mol.GetAtoms():
        if a.GetSymbol() != 'O' or a.GetTotalNumHs() != 1:
            continue
        ca = a.GetNeighbors()[0]
        hit = any(m.GetIdx() in carbonylC
                  for n in ca.GetNeighbors() if n.GetIdx() != a.GetIdx()
                  for m in n.GetNeighbors() if m.GetIdx() != ca.GetIdx())
        (peri if hit else free).append(a.GetIdx())
    return peri, free


def trim(im, bg=255, pad=6):
    a = np.asarray(im.convert('RGB'))
    mask = (a < bg - 8).any(axis=2)
    ys, xs = np.where(mask)
    if len(ys) == 0:
        return im
    return im.crop((max(xs.min() - pad, 0), max(ys.min() - pad, 0),
                    min(xs.max() + pad, a.shape[1]), min(ys.max() + pad, a.shape[0])))


def draw(s, px=1000):
    mol = Chem.MolFromSmiles(SMILES[s])
    AllChem.Compute2DCoords(mol)
    peri, free = peri_hydroxyls(mol)
    quin = [a.GetIdx() for a in mol.GetAtoms()
            if a.GetSymbol() == 'O' and a.GetTotalNumHs() == 0]
    col = {}
    col.update({i: rgb(CHEL) for i in peri})
    col.update({i: rgb(FREE) for i in free})
    col.update({i: rgb(QUIN) for i in quin})
    d = rdMolDraw2D.MolDraw2DCairo(px, int(px * 0.70))
    o = d.drawOptions()
    o.bondLineWidth = 3
    o.scaleBondWidth = False
    o.fillHighlights = True
    o.highlightRadius = 0.30
    o.minFontSize = 34
    o.maxFontSize = 34
    o.padding = 0.07
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol, highlightAtoms=list(col),
                                       highlightAtomColors=col)
    d.FinishDrawing()
    return trim(Image.open(io.BytesIO(d.GetDrawingText()))), len(free)


def main():
    style()
    hyd = rd('s4_an_hydroxyl.csv')
    geo = {r['system']: r for r in rd('s4_an_geometry.csv')}

    fig = plt.figure(figsize=(DOUBLE, 3.45))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.02, 1.00],
                          hspace=0.20, wspace=0.52,
                          left=0.085, right=0.985, top=0.90, bottom=0.14)

    for k, s in enumerate(SYS):
        ax = fig.add_subplot(gs[0, k])
        img, nfree = draw(s)
        ax.imshow(np.asarray(img))
        ax.axis('off')
        ax.set_title('%s\n(%s), %d free OH'
                     % (NAME[s].replace(' acid', ''), SHORT[s], nfree),
                     fontsize=7.0, pad=2, linespacing=1.3)
        ax.text(-0.08, 1.46, 'ABCD'[k], transform=ax.transAxes, fontsize=9,
                fontweight='bold', va='top', ha='left')

    # ---- E: measured chelate occupancy, one marker per hydroxyl oxygen
    axE = fig.add_subplot(gs[1, :2])
    for i, s in enumerate(SYS):
        rows = [r for r in hyd if r['system'] == ROW[s]]
        off = np.linspace(-0.16, 0.16, len(rows))
        for r, dx in zip(rows, off):
            occ = float(r['occ_chelated'])
            cl = 'alpha' in r['classification']
            axE.plot([i + dx, i + dx], [0, occ], color=CHEL if cl else FREE, lw=1.2,
                     zorder=1)
            axE.plot(i + dx, occ, 'o', ms=5.5, mfc=CHEL if cl else FREE,
                     mec='black', mew=0.5, zorder=3)
            axE.text(i + dx, occ + 0.075, r'$\alpha$' if cl else r'$\beta$',
                     ha='center', va='bottom', fontsize=7.5)
    axE.axhline(0, color='black', lw=0.6)
    axE.set_ylim(-0.06, 1.28)
    axE.set_yticks([0, 0.5, 1.0])
    axE.set_ylabel('chelate occupancy\n(fraction of 501 frames)')
    axE.set_xticks(range(4))
    axE.set_xticklabels([SHORT[s] for s in SYS], fontsize=7)
    axE.set_xlim(-0.5, 3.5)
    panel(axE, 'E', x=-0.155, y=1.02)

    # ---- F: measured O...O separation of the two hydroxyl oxygens
    axF = fig.add_subplot(gs[1, 2:])
    sp = [float(geo[ROW[s]]['span_mean_A']) for s in SYS]
    sd = [float(geo[ROW[s]]['span_sd_A']) for s in SYS]
    axF.bar(range(4), sp, yerr=sd, width=0.60, capsize=2.0,
            color=[C[s] for s in SYS], edgecolor='black', linewidth=0.4,
            error_kw=dict(lw=0.7))
    for i, (v, e) in enumerate(zip(sp, sd)):
        axF.text(i, v + e + 0.25, '%.2f' % v, ha='center', va='bottom', fontsize=6.8)
    axF.axhline(8.6, color='black', ls=(0, (4, 2)), lw=0.8)
    axF.text(-0.42, 8.9, 'reach needed to bridge\nboth ends (8.6 Å)',
             ha='left', va='bottom', fontsize=6.3, linespacing=1.15)
    axF.set_xticks(range(4))
    axF.set_xticklabels([SHORT[s] for s in SYS], fontsize=7)
    axF.set_ylim(0, 13.2)
    axF.set_yticks([0, 2.5, 5, 7.5, 10, 12.5])
    axF.set_ylabel('hydroxyl O···O separation (Å)')
    panel(axF, 'F', x=-0.21, y=1.02)

    handles = [plt.Line2D([], [], marker='o', ls='', ms=4.5, mfc=c, mec='black',
                          mew=0.4, label=t)
               for c, t in ((FREE, 'free hydroxyl O'), (CHEL, 'chelated hydroxyl O'),
                            (QUIN, 'quinone carbonyl O'))]
    fig.legend(handles=handles, loc='lower center', ncol=3,
               bbox_to_anchor=(0.5, -0.025), handletextpad=0.3, columnspacing=1.6)
    save(fig, 'Figure1')
    plt.close(fig)


if __name__ == '__main__':
    main()
