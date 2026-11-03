"""Verify the study geometric-bound counts from the released per-frame CSV.

Python 3 standard library only. No coordinate trajectory or scientific package is
needed for this derived-data verification. For coordinate-level reconstruction,
see frame_bound_audit.py and its recorded original trajectory/atom-pair sources.
"""
from pathlib import Path
from collections import defaultdict
import argparse
import csv
import hashlib
import json
import math

parser = argparse.ArgumentParser()
parser.add_argument('csv_file', nargs='?', type=Path,
                    default=Path(__file__).with_name('frame_bound_per_frame.csv'))
args = parser.parse_args()
with args.csv_file.open(newline='', encoding='utf-8') as fh:
    rows = list(csv.DictReader(fh))
assert len(rows) == 11022, 'Unexpected frame count'
by_run = defaultdict(list)
totals = {kind: {'frames': 0, 'bridged': 0} for kind in
          ['system_mean_span', 'run_mean_span', 'instantaneous_span']}
span_fields = {'system_mean_span': 'system_mean_hydroxyl_span_A',
               'run_mean_span': 'run_mean_hydroxyl_span_A',
               'instantaneous_span': 'instantaneous_hydroxyl_span_A'}
for row in rows:
    by_run[row['run']].append(row)
    width = float(row['cavity_width_A'])
    glu = float(row['GLU353_min_ligand_O_distance_A'])
    his = float(row['HIS524_min_ligand_O_distance_A'])
    assert all(math.isfinite(x) and x >= 0 for x in [width, glu, his])
    bridged = max(glu, his) < 3.5
    for kind, field in span_fields.items():
        span = float(row[field])
        assert math.isfinite(span) and span > 0
        if width > span + 7.0:
            totals[kind]['frames'] += 1
            totals[kind]['bridged'] += int(bridged)
assert len(by_run) == 22
for run, subset in by_run.items():
    assert sorted(int(row['frame_index']) for row in subset) == list(range(501)), run
    expected = float(subset[0]['run_mean_hydroxyl_span_A'])
    observed = math.fsum(float(row['instantaneous_hydroxyl_span_A']) for row in subset) / 501
    assert abs(expected - observed) < 1e-6, run
assert totals == {'system_mean_span': {'frames': 4552, 'bridged': 0},
                  'run_mean_span': {'frames': 4554, 'bridged': 0},
                  'instantaneous_span': {'frames': 4552, 'bridged': 0}}, totals
print(json.dumps({'status': 'passed', 'runs': len(by_run), 'frames': len(rows),
                  'csv_sha256': hashlib.sha256(args.csv_file.read_bytes()).hexdigest(),
                  'counts': totals}, indent=2))
