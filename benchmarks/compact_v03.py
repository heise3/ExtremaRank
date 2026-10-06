"""Lossless deterministic text compression, with original/archive byte hashes."""
import gzip,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
manifest=ROOT/'results/archive_v03.json'
records={r['original_path']:r for r in json.loads(manifest.read_text())['files']} if manifest.exists() else {}
folders=[ROOT/'results/v03',*[ROOT/'data/prepared'/name for name in ('KangPB','CPTAC_runs','Nutrimouse')]]
for folder in folders:
 for p in sorted(folder.rglob('*')):
  if p.suffix not in {'.csv','.tsv'} or p.stat().st_size<100000:continue
  raw=p.read_bytes();packed=gzip.compress(raw,compresslevel=9,mtime=0);assert gzip.decompress(packed)==raw
  dest=p.with_name(p.name+'.gz');dest.write_bytes(packed)
  records[str(p.relative_to(ROOT))]={'original_path':str(p.relative_to(ROOT)),'archived_path':str(dest.relative_to(ROOT)),'original_sha256':hashlib.sha256(raw).hexdigest(),'archive_sha256':hashlib.sha256(packed).hexdigest(),'original_bytes':len(raw),'archive_bytes':len(packed),'lossless_verified':True}
  p.unlink()
  for name in ('audit.json','preparation.json','feature_robustness.json'):
   j=p.parent/name
   if not j.exists():continue
   obj=json.loads(j.read_text())
   if 'outputs' in obj:
    obj['outputs']={key:dest.name if value==p.name else value for key,value in obj['outputs'].items()}
   prep=obj.get('preparation',obj)
   if prep.get('prepared_path')==p.name:prep['prepared_path']=dest.name
   obj['distribution_note']='Large execution tables compressed losslessly after validation; archive manifest binds original and compressed bytes. Reruns write ordinary tables.'
   j.write_text(json.dumps(obj,indent=2)+'\n')
manifest.write_text(json.dumps({'files':list(records.values())},indent=2)+'\n')
for r in records.values():
 packed=(ROOT/r['archived_path']).read_bytes();assert hashlib.sha256(packed).hexdigest()==r['archive_sha256'];assert hashlib.sha256(gzip.decompress(packed)).hexdigest()==r['original_sha256']
print(json.dumps({'lossless_archives_verified':len(records),'hash_mismatches':0}))
