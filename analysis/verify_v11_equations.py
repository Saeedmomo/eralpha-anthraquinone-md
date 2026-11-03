"""Reproduce V11's equation-derived numerical summaries from released records.

This verification uses Python 3's standard library. It does not require a saved
coordinate trajectory and does not rerun molecular dynamics.
"""
from collections import defaultdict
from pathlib import Path
import csv
import json
import math
import statistics

ROOT = Path(__file__).resolve().parents[1]
with (ROOT/'data/v11_geometry/v11_frame_bound_per_frame.csv').open(newline='', encoding='utf-8') as fh:
    frames = list(csv.DictReader(fh))
by_run = defaultdict(list)
for row in frames:
    by_run[row['run']].append(row)
assert len(by_run) == 22 and len(frames) == 11022
A, D0 = 2.8, 3.5
summaries = {}
for run, rows in by_run.items():
    assert len(rows) == 501
    width = statistics.fmean(float(r['cavity_width_A']) for r in rows)
    span = float(rows[0]['run_mean_hydroxyl_span_A'])
    bridge = statistics.fmean(max(float(r['GLU353_min_ligand_O_distance_A']),
                                 float(r['HIS524_min_ligand_O_distance_A'])) < D0 for r in rows)
    summaries[run] = {'mean_width_A': width, 'mean_span_A': span,
                      'reach_margin_A': span + 2*A - width,
                      'mean_contact_margin_A': span + 2*D0 - width,
                      'direct_bridge_fraction': bridge}

forbidden = [r for r, v in summaries.items() if v['mean_contact_margin_A'] < 0]
intermediate = [r for r, v in summaries.items() if v['mean_contact_margin_A'] >= 0 > v['reach_margin_A']]
inside = [r for r, v in summaries.items() if v['reach_margin_A'] >= 0]
assert len(forbidden) == 9 and all(summaries[r]['direct_bridge_fraction'] == 0 for r in forbidden)
assert len(intermediate) == 1
assert abs(summaries[intermediate[0]]['reach_margin_A'] + 0.52) < 0.005
assert abs(summaries[intermediate[0]]['direct_bridge_fraction'] - 0.076) < 0.001
assert len(inside) == 12
inside_range = [min(summaries[r]['direct_bridge_fraction'] for r in inside),
                max(summaries[r]['direct_bridge_fraction'] for r in inside)]
assert abs(inside_range[0] - 0.489) < 0.001 and abs(inside_range[1] - 0.998) < 0.001

critical = {r: (v['mean_width_A'] - v['mean_span_A'])/2 for r, v in summaries.items()}
engaged = {r: v['direct_bridge_fraction'] >= 0.4 for r, v in summaries.items()}
lower = max(critical[r] for r in summaries if engaged[r])
upper = min(critical[r] for r in summaries if not engaged[r])
assert abs(lower - 2.502) < 0.002 and abs(upper - 3.060) < 0.002
equivalent = [a for a in (2.40 + 0.01*i for i in range(81))
              if all((v['mean_span_A'] + 2*a - v['mean_width_A'] >= 0) == engaged[r]
                     for r, v in summaries.items())]
assert abs(min(equivalent) - 2.51) < 1e-9 and abs(max(equivalent) - 3.05) < 1e-9

groups = defaultdict(list)
for row in frames:
    groups[row['compound']].append(row)
glu_means, his_means = [], []
for compound in ['17b-estradiol', '2,6-DHAQ', '2,7-DHAQ', 'benzoxathiin']:
    rows = groups[compound]
    distances = []
    for field in ['GLU353_min_ligand_O_distance_A', 'HIS524_min_ligand_O_distance_A']:
        subset = [float(row[field]) for row in rows if float(row[field]) <= D0]
        assert subset
        distances.append(statistics.fmean(subset))
    glu_means.append(distances[0]); his_means.append(distances[1])
glu, his = statistics.fmean(glu_means), statistics.fmean(his_means)
assert abs(glu - 2.665) < 0.002 and abs(his - 2.913) < 0.002
assert abs(glu + his - 5.578) < 0.002

weights = {'hbond': 5.0, 'cation_pi': 3.5, 'hydrophobic': 1.5, 'ionic': 4.0, 'water_bridge': 2.0}
with (ROOT/'reference/fp_project_by_residue.csv').open(newline='', encoding='utf-8') as fh:
    scores = list(csv.DictReader(fh))
assert len(scores) == 482
errors = [abs(float(row['score']) - sum(weight*float(row['occ_'+kind])
          for kind, weight in weights.items())) for row in scores]
assert max(errors) < 2e-15
glu23 = [float(r['score']) for r in scores if r['label']=='GLU353' and r['run'].startswith('dihydroxy23')]
asp_control = [float(r['score']) for r in scores if r['label']=='ASP351' and r['run'].startswith('nativeE4D')]
assert len(glu23) == len(asp_control) == 3
assert abs(statistics.fmean(glu23) - 11.00) < 0.005
assert abs(statistics.fmean(asp_control) - 10.59) < 0.005

report = {
    'status': 'passed', 'runs': len(summaries),
    'mean_width_regimes': {'beyond_contact_bound': len(forbidden),
                          'within_contact_bound_negative_reach_margin': len(intermediate),
                          'positive_reach_margin': len(inside)},
    'positive_margin_bridge_range': inside_range,
    'critical_allowance_edges_A': [lower, upper],
    'equivalent_scanned_allowance_interval_A': [min(equivalent), max(equivalent)],
    'engaged_anchor_distance_means_A': {'GLU353': glu, 'HIS524': his, 'sum': glu+his},
    'fingerprint_score_rows': len(scores), 'maximum_score_reconstruction_error': max(errors),
    'GLU353_23_score_mean': statistics.fmean(glu23),
    'ASP351_benzoxathiin_score_mean': statistics.fmean(asp_control),
    'note': 'Mean-width regimes summarize trajectories; the separate frame-bound verifier uses instantaneous spans.',
}
print(json.dumps(report, indent=2))
