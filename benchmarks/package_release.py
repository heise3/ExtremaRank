"""Create a deterministic source archive and per-file provenance manifest."""
import argparse
import hashlib
import json
import re
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
VERSION = re.search(r'__version__ = "([^"]+)"', (ROOT/'src/extremarank/__init__.py').read_text()).group(1)


def selected_files():
    for path in sorted(ROOT.rglob('*')):
        rel = path.relative_to(ROOT)
        if not path.is_file() or any(p in {'__pycache__','build','.git','dist','.venv'} or
                                     p.endswith('.egg-info') for p in rel.parts):
            continue
        if path.name == 'RELEASE_MANIFEST.json' or path.suffix in {'.pyc','.pyo'}:
            continue
        yield path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT.parent/f'ExtremaRank-v{VERSION}-source.zip')
    args = parser.parse_args()
    paths = list(selected_files())
    records = [dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,
                    sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths]
    manifest = ROOT/'RELEASE_MANIFEST.json'
    manifest.write_text(json.dumps(dict(version=VERSION,files=records,
        note='Manifest excludes itself; archive hash is stored alongside archive.'),indent=2)+'\n')
    paths.append(manifest)
    with zipfile.ZipFile(args.output,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for path in sorted(paths):
            info=zipfile.ZipInfo('extremarank/'+str(path.relative_to(ROOT)),(2026,10,3,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=0o644 << 16
            archive.writestr(info,path.read_bytes())
    with zipfile.ZipFile(args.output) as archive:
        assert archive.testzip() is None
        for r in records:
            assert hashlib.sha256(archive.read('extremarank/'+r['path'])).hexdigest()==r['sha256']
    checksum=hashlib.sha256(args.output.read_bytes()).hexdigest()
    args.output.with_suffix('.zip.sha256').write_text(checksum+'  '+args.output.name+'\n')
    print(json.dumps(dict(archive=str(args.output),files=len(paths),bytes=args.output.stat().st_size,
                         sha256=checksum,verified=True)))


if __name__=='__main__':
    main()
