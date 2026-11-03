"""Read-only verification of V11's framewise bridging-bound statement.

No simulation input or saved trajectory is altered. Output is an audit JSON only.
Run in the existing WSL /home/saeed/mdaenv environment.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import numpy as np
from scipy.io import netcdf_file

parser = argparse.ArgumentParser()
parser.add_argument('--data-root', type=Path, default=Path('/home/saeed/md_s4'))
parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent)
args = parser.parse_args()
S = args.data_root
args.output_dir.mkdir(parents=True, exist_ok=True)
GROUPS = [
    ('1,8-DHAQ', ['danthron', 'danthron_r2', 'danthron_r3']),
    ('1,2-DHAQ', ['alizarin', 'alizarin_r2', 'alizarin_r3']),
    ('2,3-DHAQ', ['dihydroxy23', 'dihydroxy23_r2', 'dihydroxy23_r3']),
    ('17b-estradiol', ['estradiol_r1', 'estradiol_r2', 'estradiol_r3']),
    ('2,6-DHAQ', ['anthraflavic', 'anthraflavic_r2', 'anthraflavic_r3', 'anthraflavic_r5']),
    ('2,7-DHAQ', ['dihydroxy27_r1', 'dihydroxy27_r2', 'dihydroxy27_r3']),
    ('benzoxathiin', ['nativeE4D_r1', 'nativeE4D_r2', 'nativeE4D_r3']),
]
def numbers(path):
    return np.loadtxt(path, comments='#', ndmin=2)

totals = {k: {'frames': 0, 'bridged': 0} for k in ['system_mean_span', 'run_mean_span', 'instantaneous_span']}
rows = []
frame_rows = []
for name, runs in GROUPS:
    records = {r: json.loads((S / r / 'an/structure.json').read_text()) for r in runs}
    group_mean = float(np.mean([records[r]['hydroxyl_span_max'] for r in runs]))
    for run in runs:
        an = S / run / 'an'
        rec = records[run]
        pair = max(rec['hydroxyl_O_spans'], key=lambda p: p['mean'])
        atoms = [int(pair[k][1:]) - 1 for k in ['o1', 'o2']]
        path = an / 'dry_nofit.nc'
        with netcdf_file(path, 'r', mmap=True) as data:
            coords = data.variables['coordinates']
            assert coords.units in [b'angstrom', b'angstroms']
            first = np.array(coords[:, atoms[0], :], dtype=np.float64)
            second = np.array(coords[:, atoms[1], :], dtype=np.float64)
            del coords
        spans = np.linalg.norm(first - second, axis=1)
        assert abs(float(spans.mean()) - pair['mean']) < 1e-6
        assert abs(float(spans.std()) - pair['sd']) < 1e-6
        widths = numbers(S / 'stageA/span_check' / run / 'span.dat')[:, 1:5].min(axis=1)
        glu = numbers(an / 'triad_GLU353.dat')[:, 1]
        his = numbers(an / 'triad_HIS524.dat')[:, 1]
        assert len(spans) == len(widths) == len(glu) == len(his) == 501
        bridged = np.maximum(glu, his) < 3.5
        for i in range(len(spans)):
            frame_rows.append({'compound': name, 'run': run, 'frame_index': i,
                'cavity_width_A': repr(float(widths[i])),
                'instantaneous_hydroxyl_span_A': repr(float(spans[i])),
                'run_mean_hydroxyl_span_A': repr(float(pair['mean'])),
                'system_mean_hydroxyl_span_A': repr(group_mean),
                'GLU353_min_ligand_O_distance_A': repr(float(glu[i])),
                'HIS524_min_ligand_O_distance_A': repr(float(his[i]))})
        per_run = {}
        for kind, span in [('system_mean_span', group_mean), ('run_mean_span', pair['mean']), ('instantaneous_span', spans)]:
            beyond = widths > span + 7.0
            per_run[kind] = {'frames': int(beyond.sum()), 'bridged': int((beyond & bridged).sum())}
            totals[kind]['frames'] += per_run[kind]['frames']
            totals[kind]['bridged'] += per_run[kind]['bridged']
        rows.append({'compound': name, 'run': run, 'atom_pair_one_based': pair['o1'] + ' ' + pair['o2'],
                     'reconstructed_span_mean': float(spans.mean()), 'reconstructed_span_sd': float(spans.std()),
                     'system_mean_span': group_mean, 'counts': per_run,
                     'frames_equal_to_contact_cutoff': int(((glu == 3.5) | (his == 3.5)).sum())})

csv_path = args.output_dir / 'v11_frame_bound_per_frame.csv'
with csv_path.open('w', newline='', encoding='utf-8') as fh:
    writer = csv.DictWriter(fh, fieldnames=list(frame_rows[0]))
    writer.writeheader()
    writer.writerows(frame_rows)
report = {
    'method': 'Frame widths from the four-pair span.dat minimum; direct bridging from archived minimum ligand-oxygen GLU353/HIS524 distances <3.5 Å. Instantaneous hydroxyl spans reconstructed from the two atom indices in structure.json and saved whole-molecule dry_nofit.nc coordinates. Reconstructed means/SDs match archived structure.json within1e-6 Å. Same22-run set as check_v11.py.',
    'no_input_or_trajectory_mutation': True,
    'total_runs': len(rows), 'total_frames': 501 * len(rows),
    'totals': totals, 'run_records': rows,
    'derived_per_frame_csv': csv_path.name,
    'derived_per_frame_csv_sha256': hashlib.sha256(csv_path.read_bytes()).hexdigest(),
    'extraction_script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
out = args.output_dir / 'frame_bound_audit.json'
out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({k: v for k, v in report.items() if k != 'run_records'}, ensure_ascii=False, indent=2))
