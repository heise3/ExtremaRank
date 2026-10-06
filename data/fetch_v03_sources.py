"""Verify bundled fixed source bytes; retrieve missing public data on request."""
import argparse, hashlib, json, tarfile, tempfile, urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parent/'source/v03'

def verify(path,sha):
 assert hashlib.sha256(path.read_bytes()).hexdigest()==sha, 'Source hash mismatch: '+str(path)

def main():
 p=argparse.ArgumentParser();p.add_argument('--offline',action='store_true');a=p.parse_args();checks=0
 records=json.loads((ROOT/'download_manifest.json').read_text())['sources']
 kang=json.loads((ROOT/'Kang/manifest.json').read_text());records.append({'url':kang['url'],'path':'Kang/EH2259.rds','sha256':kang['sha256']})
 for r in records:
  if r.get('bundled') is False:continue
  dest=ROOT/r['path']
  if not dest.exists():
   if a.offline:raise FileNotFoundError(dest)
   dest.parent.mkdir(parents=True,exist_ok=True);urllib.request.urlretrieve(r['url'],dest)
  verify(dest,r['sha256']);checks+=1
 for directory,archive_record in [('mixOmics',json.loads((ROOT/'mixOmics/manifest.json').read_text())),('muscat',json.loads((ROOT/'muscat/manifest.json').read_text()))]:
  missing=[r for r in archive_record['files'] if not (ROOT/directory/r['path']).exists()]
  if missing:
   if a.offline:raise FileNotFoundError(missing)
   with tempfile.TemporaryDirectory() as tmp:
    arc=Path(tmp)/'source.tar.gz';urllib.request.urlretrieve(archive_record['source_url'],arc);verify(arc,archive_record['package_archive_sha256'])
    with tarfile.open(arc) as source:
     for r in missing:
      (ROOT/directory/r['path']).write_bytes(source.extractfile(r['member']).read())
  for r in archive_record['files']:verify(ROOT/directory/r['path'],r['sha256']);checks+=1
 print(json.dumps({'verified_source_files':checks,'offline':a.offline,'hash_mismatches':0}))
if __name__=='__main__':main()
