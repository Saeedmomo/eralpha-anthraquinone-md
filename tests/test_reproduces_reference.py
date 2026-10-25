# -*- coding: utf-8 -*-
"""Verification tests for interaction_fingerprint.py.

The reference CSVs in ../reference were produced from the published antagonist-state
trajectories. The point of these tests is that a change to the script, the MDAnalysis
version or the toolchain cannot silently alter the numbers: the descriptor must return
the reference matrix with a maximum absolute difference of zero.

They need the trajectories, which are large and are not in this repository. Point the
tests at them with

    export ERA_MD_ROOT=/path/to/md_s4
    pytest -q

Without that variable the trajectory-dependent tests skip and the pure-logic tests
still run, so `pytest` is useful even on a clone with no data.
"""
import os
import subprocess
import sys

import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
REF = os.path.join(PKG, 'reference')
SCRIPT = os.path.join(PKG, 'interaction_fingerprint.py')
ROOT = os.environ.get('ERA_MD_ROOT')

COLS = ['score', 'occ_hbond', 'occ_cation_pi', 'occ_hydrophobic',
        'occ_ionic', 'occ_water_bridge']
# three antagonist runs that exercise all five contact types
RUNS = ['dihydroxy23', 'alizarin_r3', 'estradiol_r1']

needs_data = pytest.mark.skipif(
    not ROOT or not os.path.isdir(ROOT),
    reason='set ERA_MD_ROOT to the directory holding the run folders')


def test_mdanalysis_is_new_enough():
    """GROMACS 2024+ writes tpx 138; MDAnalysis below 2.11 cannot read it."""
    mda = pytest.importorskip('MDAnalysis')
    got = tuple(int(x) for x in mda.__version__.split('.')[:3] if x.isdigit())
    assert got >= (2, 11, 0), (
        'MDAnalysis %s cannot parse tpx 138; install MDAnalysis>=2.11.0' % mda.__version__)


def test_script_imports_and_shows_help():
    out = subprocess.run([sys.executable, SCRIPT, '--help'],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    for flag in ('--criteria', '--offset', '--receptor-resid-range', '--out-prefix'):
        assert flag in out.stdout, 'missing %s in --help' % flag


@needs_data
@pytest.mark.parametrize('criteria', ['project', 'github'])
def test_reproduces_reference_matrix(tmp_path, criteria):
    """The descriptor must return the published antagonist matrix exactly."""
    prefix = str(tmp_path / ('out_%s' % criteria))
    out = subprocess.run(
        [sys.executable, SCRIPT, '--root', ROOT, '--runs'] + RUNS +
        ['--criteria', criteria, '--out-prefix', prefix],
        capture_output=True, text=True)
    assert out.returncode == 0, out.stderr

    new = pd.read_csv(prefix + '_by_residue.csv')
    old = pd.read_csv(os.path.join(REF, 'fp_%s_by_residue.csv' % criteria))
    old = old[old['run'].isin(RUNS)]
    assert len(new) == len(old), 'row count %d vs reference %d' % (len(new), len(old))

    key = ['run', 'resid', 'label']
    m = old.merge(new, on=key, how='outer', suffixes=('_old', '_new'), indicator=True)
    assert (m['_merge'] == 'both').all(), 'rows present in only one matrix'
    for c in COLS:
        worst = (m[c + '_old'] - m[c + '_new']).abs().max()
        assert worst == 0.0, '%s differs from the reference by %.3e' % (c, worst)


@needs_data
def test_second_chain_is_not_renumbered(tmp_path):
    """A system with a coactivator peptide must not map it onto receptor numbering.

    The agonist-state systems carry the peptide at topology resid 246-256. Without
    --receptor-resid-range those residues would be reported as ERalpha 551-561, which
    do not exist.
    """
    run = 'agopep_26dhaq_r1'
    if not os.path.isdir(os.path.join(ROOT, run)):
        pytest.skip('%s not present' % run)
    prefix = str(tmp_path / 'ago')
    out = subprocess.run(
        [sys.executable, SCRIPT, '--root', ROOT, '--runs', run,
         '--criteria', 'project', '--receptor-resid-range', '1', '245',
         '--out-prefix', prefix],
        capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    df = pd.read_csv(prefix + '_by_residue.csv')
    assert 'chain' in df.columns
    receptor = df[df['chain'] == 'receptor']
    assert receptor['era'].notna().all()
    assert receptor['resid'].between(1, 245).all()
    other = df[df['chain'] == 'other']
    # whether or not the peptide is contacted, nothing outside the receptor may be
    # given a receptor residue number
    assert other['era'].isna().all() or (other['era'].astype(str) == '').all()
    assert not ((df['resid'] > 245) & (df['chain'] == 'receptor')).any()


@needs_data
def test_offset_is_configurable(tmp_path):
    """Changing --offset must shift the reported numbering by exactly that much."""
    run = RUNS[0]
    a, b = str(tmp_path / 'a'), str(tmp_path / 'b')
    for prefix, off in ((a, 305), (b, 0)):
        out = subprocess.run(
            [sys.executable, SCRIPT, '--root', ROOT, '--runs', run,
             '--criteria', 'project', '--offset', str(off), '--out-prefix', prefix],
            capture_output=True, text=True)
        assert out.returncode == 0, out.stderr
    da = pd.read_csv(a + '_by_residue.csv').sort_values('resid').reset_index(drop=True)
    db = pd.read_csv(b + '_by_residue.csv').sort_values('resid').reset_index(drop=True)
    assert (da['era'] - db['era'] == 305).all()
    assert (da['score'].values == db['score'].values).all(), 'offset changed the scores'
