"""Matched exact search timings, frozen separated and difficult inputs."""
import json
from pathlib import Path
import time
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from extremarank.unpaired import audit_welch

records=[]
for n in (20,50,100):
    half=n//2;groups=['T']*half+['R']*half
    rows=[[100+i/10,*[1+((i+j)%5)/10 for j in range(7)]] for i in range(half)]+[[i/10,*[1+((i+j)%5)/10 for j in range(7)]] for i in range(half)]
    for method in ('enumeration','branch-bound'):
        start=time.perf_counter();r=audit_welch(rows,groups,tuple('abcdefgh'),'T','R',1,3,1000000,'up',method,1000000)
        assert r.status=='CERTIFIED'
        records.append({'case':'separated-fixed-input','samples':n,'features':8,'budget':3,'method':method,
                        'status':r.status,'seconds':time.perf_counter()-start,'scenarios_checked':r.subsets_checked,
                        'feasible_deletion_sets':r.feasible_deletion_sets,'nodes':r.nodes,'query_bounds_pruned':r.query_bounds_pruned})
# A deterministic near-boundary case with an intentionally small resource
# limit remains reported as obtained, including UNRESOLVED if encountered.
rows=[[1+(i%3)/10,1+(i%3)/10,0] for i in range(10)]+[[0+(i%3)/10,0+(i%3)/10,0] for i in range(10)]
r=audit_welch(rows,['T']*10+['R']*10,('a','b','c'),'T','R',1,3,5,'up','branch-bound',5)
records.append({'case':'tied-boundary-finite-limit','method':'branch-bound','status':r.status,'scenarios_checked':r.subsets_checked,'nodes':r.nodes})
out={'records':records,'scope':'One-machine matched exact arithmetic; separated stable synthetic inputs demonstrate pruning, not typical or worst-case speed. Difficult outcomes retained. No timing comparison with DESeq2/edgeR/limma.'}
(ROOT/'results/welch_pruning_v03.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
