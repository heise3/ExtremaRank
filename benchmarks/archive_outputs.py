"""Compress executed TSVs for distribution; verify lossless round-trip first."""
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    records = []
    for path in sorted((ROOT/'results').rglob('*.tsv')):
        raw = path.read_bytes()
        target = path.with_suffix('.tsv.gz')
        with target.open('wb') as file:
            with gzip.GzipFile(filename='',mode='wb',fileobj=file,mtime=0) as compressed:
                compressed.write(raw)
        restored = gzip.decompress(target.read_bytes())
        assert restored == raw
        records.append({'original_path':str(path.relative_to(ROOT)),
                        'archived_path':str(target.relative_to(ROOT)),
                        'original_bytes':len(raw),'archived_bytes':target.stat().st_size,
                        'original_sha256':hashlib.sha256(raw).hexdigest(),
                        'archive_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
                        'lossless_verified':True})
        path.unlink()
        if path.parent.name.endswith('_cli'):
            audit_path = path.parent/'audit.json'
            audit = json.loads(audit_path.read_text())
            for key,value in list(audit['outputs'].items()):
                if value == path.name:
                    audit['outputs'][key] = target.name
            audit['distribution_note'] = 'Executed TSVs compressed losslessly after the run; reproduction writes ordinary TSV.'
            audit_path.write_text(json.dumps(audit,indent=2)+'\n')
    output = ROOT/'results/archive_manifest.json'
    if records:
        output.write_text(json.dumps({'kind':'lossless distribution compression',
                                    'files':records},indent=2)+'\n')
    print(f'Compressed {len(records)} executed tables losslessly')


if __name__=='__main__':
    main()
