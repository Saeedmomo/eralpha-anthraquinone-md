"""Shared style and data access for the manuscript figures.

Everything is read from the delivered Stage-4 CSV tables and from the per-frame
.dat/.npz files copied out of the simulation directories. Nothing is typed in by
hand, so a figure cannot drift away from the tables it is supposed to show.
"""
import csv
import json
import os

import matplotlib as mpl
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
TAB = r'D:\Chua files\Third_article_MDS\Stage4_package\analysis_tables'
OUT = r'D:\Chua files\Third_article_MDS\Manuscript\figures'
OUTSI = r'D:\Chua files\Third_article_MDS\Manuscript\figures_SI'

SYS = ['danthron', 'alizarin', 'dihydroxy23', 'anthraflavic']
# the row label used in every delivered CSV
ROW = {'danthron': 'Danthron (1,8-diOH)',
       'alizarin': 'Alizarin (1,2-diOH)',
       'dihydroxy23': '2,3-dihydroxyanthraquinone',
       'anthraflavic': 'Anthraflavic acid (2,6-diOH)'}
NAME = {'danthron': 'Danthron', 'alizarin': 'Alizarin',
        'dihydroxy23': '2,3-DHAQ', 'anthraflavic': 'Anthraflavic acid'}
SHORT = {'danthron': '1,8-diOH', 'alizarin': '1,2-diOH',
         'dihydroxy23': '2,3-diOH', 'anthraflavic': '2,6-diOH'}
SMILES = {'danthron': 'O=C1c2cccc(O)c2C(=O)c2c(O)cccc21',
          'alizarin': 'O=C1c2ccccc2C(=O)c2c1ccc(O)c2O',
          'dihydroxy23': 'O=C1c2ccccc2C(=O)c2cc(O)c(O)cc21',
          'anthraflavic': 'O=C1c2ccc(O)cc2C(=O)c2ccc(O)cc21'}

# Okabe-Ito; distinguishable in colour, in greyscale and to colour-blind readers
C = {'danthron': '#0072B2', 'alizarin': '#E69F00',
     'dihydroxy23': '#009E73', 'anthraflavic': '#D55E00'}
TRIAD_C = {'GLU353': '#0072B2', 'ARG394': '#009E73', 'HIS524': '#CC79A7'}
GREY = '#4D4D4D'

SINGLE, DOUBLE = 3.33, 7.0        # ACS column widths, inches


def style():
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,
        'axes.labelsize': 8,
        'axes.titlesize': 8,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 7,
        'legend.frameon': False,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 2.5,
        'ytick.major.size': 2.5,
        'lines.linewidth': 0.9,
        'patch.linewidth': 0.5,
        'savefig.dpi': 600,
        'figure.dpi': 110,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })


def panel(ax, letter, x=-0.16, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9,
            fontweight='bold', va='bottom', ha='left')


def rd(fn):
    with open(os.path.join(TAB, fn), encoding='utf-8') as fh:
        return list(csv.DictReader(fh))


def dat(s, fn):
    """two-column cpptraj-style .dat -> (x, y) float arrays"""
    a = np.loadtxt(os.path.join(DATA, s, fn), comments='#')
    return a[:, 0], a[:, 1]


def js(s, fn):
    with open(os.path.join(DATA, s, fn), encoding='utf-8') as fh:
        return json.load(fh)


def npz(s, fn):
    return np.load(os.path.join(DATA, s, fn))


def save(fig, name, si=False):
    d = OUTSI if si else OUT
    os.makedirs(d, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(d, '%s.%s' % (name, ext)),
                    bbox_inches='tight', pad_inches=0.02)
    print('   wrote %s.png / .pdf  (%s)' % (name, d))


def repel(xy, labels, ax, iters=400, pad=0.018, fixed=None, fixed_r=0.030):
    """Push text anchors apart so labels in a crowded 2-D map do not overlap.

    Returns a label position for each point. The original F12 pocket map drew the
    residue name at the point itself, which collided whenever two side chains
    projected close together; this does the obvious simple thing instead - treat
    the label boxes as repelling particles and let them settle.
    """
    xy = np.asarray(xy, float)
    pos = xy.copy()
    # estimated half-size of a label box in axes-fraction units
    w = np.array([0.5 * (0.011 * len(t.split('\n')[0]) + 0.012) for t in labels])
    h = np.full(len(labels), 0.030)
    for _ in range(iters):
        moved = False
        for i in range(len(pos)):
            for j in range(i + 1, len(pos)):
                dx = pos[j, 0] - pos[i, 0]
                dy = pos[j, 1] - pos[i, 1]
                ox = (w[i] + w[j] + pad) - abs(dx)
                oy = (h[i] + h[j] + pad) - abs(dy)
                if ox > 0 and oy > 0:
                    moved = True
                    if ox / (w[i] + w[j]) < oy / (h[i] + h[j]):
                        s = 0.5 * ox * (1 if dx >= 0 else -1)
                        pos[i, 0] -= s
                        pos[j, 0] += s
                    else:
                        s = 0.5 * oy * (1 if dy >= 0 else -1)
                        pos[i, 1] -= s
                        pos[j, 1] += s
        # the drawn molecule is an obstacle too: a residue label that lands on top of
        # the ring skeleton is just as unreadable as one that lands on another label
        if fixed is not None and len(fixed):
            for i in range(len(pos)):
                for fx, fy in fixed:
                    dx, dy = pos[i, 0] - fx, pos[i, 1] - fy
                    ox = (w[i] + fixed_r) - abs(dx)
                    oy = (h[i] + fixed_r) - abs(dy)
                    if ox > 0 and oy > 0:
                        moved = True
                        if ox / (w[i] + fixed_r) < oy / (h[i] + fixed_r):
                            pos[i, 0] += ox * (1 if dx >= 0 else -1)
                        else:
                            pos[i, 1] += oy * (1 if dy >= 0 else -1)
        # weak spring back to the true position so labels stay near their point
        pos += 0.05 * (xy - pos)
        pos = np.clip(pos, 0.035, 0.965)
        if not moved:
            break
    return pos
