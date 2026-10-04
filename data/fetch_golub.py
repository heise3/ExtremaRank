"""Restore pinned public multtest files; never install or execute downloaded code."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT/'data/source/Golub'
parser = argparse.ArgumentParser()
parser.add_argument('--offline', action='store_true')
args = parser.parse_args()
provenance = json.loads((DIRECTORY/'provenance.json').read_text())
missing = []
for item in provenance['source_files']:
    path = DIRECTORY/item['file']
    if path.exists():
        if hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
            raise RuntimeError('Pinned source hash mismatch; refusing drift')
    else:
        missing.append(item)
if missing:
    if args.offline:
        raise RuntimeError('Missing pinned public source files in offline mode')
    with urlopen(provenance['package_url'], timeout=60) as response:
        raw = response.read()
    if hashlib.sha256(raw).hexdigest() != provenance['package_sha256']:
        raise RuntimeError('Public package hash drift; refusing substitution')
    with tarfile.open(fileobj=io.BytesIO(raw),mode='r:gz') as archive:
        for item in missing:
            member = archive.extractfile(item['package_member'])
            if member is None:
                raise RuntimeError('Expected source member missing')
            data = member.read()
            if hashlib.sha256(data).hexdigest() != item['sha256']:
                raise RuntimeError('Extracted source hash mismatch')
            (DIRECTORY/item['file']).write_bytes(data)
print(f"Verified {len(provenance['source_files'])} pinned public Golub source files")
