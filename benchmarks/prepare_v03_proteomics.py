"""Frozen complete-positive-run protein preparation; no imputation."""
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
source=ROOT/'data/source/v03/proteinGroups.txt'
out=ROOT/'data/prepared/CPTAC_runs';out.mkdir(parents=True,exist_ok=True)
samples=['A7','A8','A9','B7','B8','B9']
columns=['Intensity 6'+s[0]+'_'+s[1:] for s in samples]
rows=[];excluded={'provider_flags':0,'missing_or_nonpositive':0}
with source.open() as f:
    for row in csv.DictReader(f,delimiter='\t'):
        if any(row[k]=='+' for k in ('Reverse','Potential contaminant','Only identified by site')):
            excluded['provider_flags']+=1;continue
        values=[float(row[c]) if row[c] else math.nan for c in columns]
        if any(not math.isfinite(x) or x<=0 for x in values):
            excluded['missing_or_nonpositive']+=1;continue
        rows.append([row['Protein IDs'],*[format(math.log2(x),'.17g') for x in values]])
assert len({r[0] for r in rows})==len(rows)
with (out/'matrix.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['feature_id',*samples]);w.writerows(rows)
with (out/'metadata.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['sample_id','group','donor_id','unit_type']);w.writerows([[s,s[0],s,'technical_measurement_run'] for s in samples])
manifest={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'retained_proteins':len(rows),'excluded':excluded,
          'source_columns':columns,'transform':'log2 supplied Intensity; complete positive values in all runs; no imputation or additional normalization',
          'independence_claim':'Technical spike-in measurement runs only; no independent patient/donor claim',
          'contract':'benchmarks/v03_data_contract.json'}
(out/'preparation.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest))
