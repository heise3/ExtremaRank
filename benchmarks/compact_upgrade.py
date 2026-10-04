"""Compact new public inputs/results with verified, deterministic gzip bytes."""
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
records = []
paths = [ROOT/'data/prepared/Golub/matrix.csv', ROOT/'data/prepared/Golub/matrix.binary64',
         ROOT/'data/prepared/Golub/limma/refits.tsv']
paths += list((ROOT/'results/upgrade').rglob('*.csv'))
paths += list((ROOT/'results/upgrade').rglob('*.tsv'))
for path in sorted(paths):
    if not path.exists():
        continue
    raw = path.read_bytes()
    target = path.with_name(path.name+'.gz')
    packed = gzip.compress(raw, compresslevel=9, mtime=0)
    assert gzip.decompress(packed) == raw
    target.write_bytes(packed)
    records.append(dict(original_path=str(path.relative_to(ROOT)),archived_path=str(target.relative_to(ROOT)),
                        original_sha256=hashlib.sha256(raw).hexdigest(),archive_sha256=hashlib.sha256(packed).hexdigest(),
                        original_bytes=len(raw),archive_bytes=len(packed),lossless_verified=True))
    path.unlink()
    audit_path = path.parent/'audit.json'
    if audit_path.exists():
        audit = json.loads(audit_path.read_text())
        for key,value in list(audit.get('outputs',{}).items()):
            if value == path.name:
                audit['outputs'][key] = target.name
        if audit.get('preparation',{}).get('prepared_path') == path.name:
            audit['preparation']['prepared_path'] = target.name
        audit['distribution_note'] = 'Large CSV/TSV files compressed losslessly after execution; reruns write ordinary tables.'
        audit_path.write_text(json.dumps(audit,indent=2)+'\n')
    preparation_path = path.parent/'preparation.json'
    if preparation_path.exists():
        prep = json.loads(preparation_path.read_text())
        if prep.get('prepared_path') == path.name:
            prep['prepared_path'] = target.name
        preparation_path.write_text(json.dumps(prep,indent=2)+'\n')
manifest_path = ROOT/'results/upgrade_archive_manifest.json'
old_records = json.loads(manifest_path.read_text())['files'] if manifest_path.exists() else []
merged = {r['original_path']:r for r in old_records}
merged.update({r['original_path']:r for r in records})
manifest_path.write_text(json.dumps(dict(files=list(merged.values())),indent=2)+'\n')
for record in merged.values():
    packed = (ROOT/record['archived_path']).read_bytes()
    assert hashlib.sha256(packed).hexdigest() == record['archive_sha256']
    assert hashlib.sha256(gzip.decompress(packed)).hexdigest() == record['original_sha256']
print(json.dumps(dict(new_files=len(records),verified_files=len(merged))))
