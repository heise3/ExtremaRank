"""Validate R-only benchmark records and attach external macOS peak RSS."""
import argparse
import json
from pathlib import Path
import re

parser = argparse.ArgumentParser()
parser.add_argument('--timing-dir',type=Path,required=True)
parser.add_argument('--results',type=Path,default=Path('results/v06/mac'))
args = parser.parse_args()
records = {}
for path in sorted(args.results.glob('*.json')):
    value = json.loads(path.read_text())
    if 'method' not in value:
        continue
    assert value['exact_match'] is True
    source = args.timing_dir / (path.stem+'.time')
    text = source.read_text()
    matched = re.search(r'(\d+)\s+maximum resident set size',text)
    assert matched, f'Missing OS peak measurement: {source}'
    value['peak_rss_bytes'] = int(matched[1])
    value['memory_measurement'] = 'macOS /usr/bin/time -l maximum resident set size, bytes; one fresh process per case/method'
    metrics = value.get('computation') or {}
    if value['method'].startswith('metal'):
        assert 0 < metrics['gpu_buffer_bytes'] <= value['gpu_budget_bytes']
        assert metrics['gpu_commands'] > 0
    path.write_text(json.dumps(value,indent=2)+'\n')
    (args.results / (path.stem+'.time.txt')).write_text(text)
    records[(value['case'],value['method'])] = value

comparison = []
for name in ('sparse-100k','sparse-250k','kang-real'):
    group = [z for (case,method),z in records.items() if case==name]
    for field in ('cells','features','units','stored_values','input_sha256','expected_hash','expected_pipeline_hash'):
        assert len({z[field] for z in group}) == 1, f'Mismatched {name} {field}'
    old,new = [records[(name,z)] for z in ('released-0.5','native')]
    whole,bounded = [records[(name,z)] for z in ('metal-whole','metal-64MiB' if name=='kang-real' else 'metal-8MiB')]
    comparison.append(dict(case=name,
        pipeline_speedup_vs_released_v05=old['pipeline_median_seconds']/new['pipeline_median_seconds'],
        peak_rss_reduction_vs_released_v05=1-new['peak_rss_bytes']/old['peak_rss_bytes'],
        metal_buffer_reduction_vs_whole=1-bounded['computation']['gpu_buffer_bytes']/whole['computation']['gpu_buffer_bytes']))
(args.results / 'comparisons.json').write_text(json.dumps(comparison,indent=2)+'\n')
print(json.dumps(comparison,indent=2))
