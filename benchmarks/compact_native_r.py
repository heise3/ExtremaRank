"""Losslessly compress generated v0.4 tables and deduplicate frozen inputs."""
import gzip
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]/'results/v04'
for dataset in ('golub','kang-b','cptac','nutrimouse'):
    original=ROOT/f'{dataset}-up/frozen_hex.tsv'
    raw=original.read_bytes() if original.exists() else gzip.decompress(original.with_name(original.name+'.gz').read_bytes())
    for direction in ('down','absolute'):
        duplicate=ROOT/f'{dataset}-{direction}/frozen_hex.tsv'
        if duplicate.exists():
            assert duplicate.read_bytes()==raw
            duplicate.unlink()
for path in sorted(ROOT.rglob('*.tsv')):
    if path.stat().st_size < 250000:
        continue
    raw=path.read_bytes(); compressed=path.with_name(path.name+'.gz')
    content=gzip.compress(raw,compresslevel=9,mtime=0)
    assert gzip.decompress(content)==raw
    if compressed.exists():assert gzip.decompress(compressed.read_bytes())==raw
    compressed.write_bytes(content);path.unlink()
print('Generated tables compressed losslessly; each dataset retains one frozen hex matrix.')
