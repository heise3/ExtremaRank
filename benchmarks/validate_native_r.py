"""Independent Fraction oracle for the native R/C++ kernel (development only).

R is invoked only by this validation script. The installed R package never
starts Python or any other interpreter. Frozen golden files can be replayed
by R alone in R CMD check.
"""
import argparse
from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import random
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tests'))
sys.path.insert(0,str(ROOT/'src'))
from test_unpaired import oracle_rank, oracle_scores, feasible_deletions


def scores(rows, groups, paired, deleted=()):
    if not paired:
        return oracle_scores(rows, groups, deleted)
    kept=[row for i,row in enumerate(rows) if i not in deleted]
    return oracle_scores(kept+[[0]*len(rows[0])]*2,['T']*len(kept)+['R']*2)


def ranks(rows, groups, ids, paired, direction, deleted=()):
    if not paired:
        return oracle_rank(rows,groups,ids,direction,deleted)
    kept=[row for i,row in enumerate(rows) if i not in deleted]
    return oracle_rank(kept+[[0]*len(ids)]*2,['T']*len(kept)+['R']*2,ids,direction)


def compare(a,b):
    sa,qa=a; sb,qb=b
    if sa!=sb:return (sa>sb)-(sa<sb)
    if not sa:return 0
    if qa is None or qb is None:return sa*((qa is None)-(qb is None))
    return sa*((qa>qb)-(qa<qb))


def unpack(z):
    return z['sign'],F(int(z['numerator']),int(z['denominator'])) if int(z['denominator']) else None


def make_cases(n):
    rng=random.Random(40606); out=[]
    for _ in range(n):
        groups=['T']*rng.randint(2,4)+['R']*rng.randint(2,4)
        out.append(dict(rows=[[rng.randint(-20,20)/8 for _ in range(4)] for g in groups],groups=groups))
    tiny=float.fromhex('0x0.0000000000001p-1022')
    for values in ([1e300,-1e300,tiny,-tiny,0.,1e-200,-1e-200,1.], [1e300]*4+[-1e300]*4,
                   [tiny,2*tiny,3*tiny,0.,-tiny,-2*tiny,-3*tiny,0.], [0.]*8):
        out.append(dict(rows=[[v,float(2**40+i),(i-4)*tiny,0.] for i,v in enumerate(values)],groups=['T']*4+['R']*4))
    return out


WORKER = r'''
library(extremarank)
a <- commandArgs(TRUE); source(a[1]); out <- list()
for (c in seq_along(cases)) {
  z <- cases[[c]]; x <- t(z$rows); ids <- c("d","b","a","c"); rownames(x)<-ids; colnames(x)<-paste0("u",seq_len(ncol(x)))
  for (paired in c(TRUE,FALSE)) {
    labels <- if(paired) integer(ncol(x)) else z$labels
    budget<-min(2L,ncol(x)-4L)
    bounds<-list()
    for(r in seq_len(budget)) {
      alloc <- if(paired) list(ncol(x)-r) else Filter(function(n) all(n>=2), lapply(0:r,function(d) c(sum(labels==0)-d,sum(labels==1)-(r-d))))
      for(need in alloc) for(forced in 0:2) {
        inc<-if(forced) 1L else integer(); opt<-setdiff(seq_len(ncol(x)),c(inc,if(forced==2L) 2L else integer()))
        nn<-need-tabulate(labels[inc]+1L,nbins=length(need))
        if(any(nn<0) || any(nn>tabulate(labels[opt]+1L,nbins=length(need)))) next
        b<-extremarank:::cpp_bounds(x,ids,labels,paired,inc,opt,nn)
        bounds[[length(bounds)+1L]]<-list(retained=need,included=inc,optional=opt,bounds=b)
      }
    }
    for(direction in c("up","down","absolute")) {
      audits<-lapply(c(1L,3L,100000L),function(limit) extremarank:::cpp_audit(x,ids,colnames(x),labels,paired,2L,budget,direction,1:4,limit,limit,1L))
      out[[length(out)+1L]]<-list(case=c,paired=paired,direction=direction,budget=budget,
        baseline=extremarank:::cpp_evaluate(x,ids,labels,paired,direction,integer()),audits=audits,
        bounds=if(direction=="up") bounds else list())
    }
  }
}
jsonlite::write_json(out,a[2],auto_unbox=TRUE,digits=NA,null="null",pretty=FALSE)
'''


def main():
    p=argparse.ArgumentParser();p.add_argument('--cases',type=int,default=1000);p.add_argument('--work',type=Path,required=True)
    p.add_argument('--output',type=Path,default=ROOT/'results/v04/native_fraction_validation.json')
    p.add_argument('--golden',action='store_true');a=p.parse_args();a.work.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();cases=make_cases(a.cases)
    # Hexadecimal R literals preserve every bit, including extreme exponents.
    spec='cases <- list('+','.join('list(rows=matrix(c('+','.join(float(v).hex() for row in c['rows'] for v in row)+'),ncol=4,byrow=TRUE),labels=c('+','.join('0L' if g=='T' else '1L' for g in c['groups'])+'))' for c in cases)+')\n'
    (a.work/'cases.R').write_text(spec);(a.work/'worker.R').write_text(WORKER)
    subprocess.run(['Rscript',str(a.work/'worker.R'),str(a.work/'cases.R'),str(a.work/'native.json')],check=True)
    results=json.loads((a.work/'native.json').read_text());checks=0;ids=('d','b','a','c');golden=[]
    for z in results:
        c=cases[z['case']-1];rows,groups=c['rows'],c['groups'];paired=z['paired'];direction=z['direction'];budget=z['budget']
        base=ranks(rows,groups,ids,paired,direction);base_scores=scores(rows,groups,paired)
        assert z['baseline']['order']==[j+1 for j in base]
        for native,truth in zip(z['baseline']['scores'],base_scores):
            assert unpack(native)==truth[:2],(z['case'],native,truth)
            assert F(native['effect'])==truth[2];checks+=2
        deletions=list(itertools.chain.from_iterable(itertools.combinations(range(len(rows)),r) for r in range(1,budget+1))) if paired else list(feasible_deletions(groups,budget))
        minima={(j,p):None for j in range(4) for p in ('membership','direction')}
        truth_cache={d:scores(rows,groups,paired,d) for d in deletions}
        rank_cache={d:ranks(rows,groups,ids,paired,direction,d) for d in deletions}
        original=set(base[:2])
        for d in deletions:
            selected=set(rank_cache[d][:2])
            for j in range(4):
                for prop,changed in (('membership',(j in selected)!=(j in original)),('direction',truth_cache[d][j][0]!=base_scores[j][0])):
                    if changed and minima[j,prop] is None:minima[j,prop]=len(d)
        for limit,report in zip((1,3,100000),z['audits']):
            for j,row in enumerate(report['features']):
                for prop in ('membership','direction'):
                    r=row[prop];truth=minima[j,prop]
                    assert r['minimum_change_lower_bound'] <= (truth or budget+1)
                    if r['status']=='CERTIFIED':assert truth is None
                    if r['minimum_change_upper_bound'] is not None:assert truth is not None and truth<=r['minimum_change_upper_bound']
                    if r['exact_minimum_change'] is not None:assert r['exact_minimum_change']==truth
                    if limit==100000:assert r['status']==('REFUTED' if truth else 'CERTIFIED') and r['exact_minimum_change']==truth
                    if r['witness']:
                        deleted=r['witness']['deleted_indices'];deleted=[deleted] if isinstance(deleted,int) else deleted
                        d=tuple(i-1 for i in deleted)
                        assert d in deletions
                        assert r['witness']['topk']==[ids[j] for j in rank_cache[d][:2]]
                        assert ((j in rank_cache[d][:2])!=(j in original)) if prop=='membership' else truth_cache[d][j][0]!=base_scores[j][0]
                    checks+=1
        for node in z['bounds']:
            inc=node['included'];opt=node['optional'];inc=[inc] if isinstance(inc,int) else inc;opt=[opt] if isinstance(opt,int) else opt
            inc={i-1 for i in inc};allowed=inc|{i-1 for i in opt};retained=node['retained'];retained=[retained] if isinstance(retained,int) else retained
            for d in deletions:
                kept=set(range(len(rows)))-set(d)
                if not inc<=kept<=allowed:continue
                counts=[len(kept)] if paired else [sum(groups[i]==g for i in kept) for g in ('T','R')]
                if counts!=retained:continue
                for b,t in zip(node['bounds'],truth_cache[d]):
                    assert compare(unpack(b['lower']),t[:2])<=0 and compare(t[:2],unpack(b['upper']))<=0,(z['case'],node,t)
                    checks+=1
        if a.golden:
            for j in range(4):golden.append(dict(case=z['case'],paired=paired,direction=direction,budget=budget,feature_id=ids[j],baseline_rank=base.index(j)+1,
                membership_minimum=minima[j,'membership'] or 0,direction_minimum=minima[j,'direction'] or 0,
                sign=base_scores[j][0],numerator=str(base_scores[j][1].numerator) if base_scores[j][1] is not None else '0',denominator=str(base_scores[j][1].denominator) if base_scores[j][1] is not None else '0',effect=str(base_scores[j][2])))
    report=dict(cases=len(cases),design_direction_problems=len(results),independent_checks=checks,failures=0,seed=40606,
        oracle='Fraction means and centered sample variances; complete shared deletion subsets; full ranks, exact scalar scores, membership/sign minima, actual witnesses, capped intervals and root/forced-inclusion enclosures',seconds=time.perf_counter()-start)
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    if a.golden:
        import csv
        folder=ROOT/'r/extremarank/inst/validation';folder.mkdir(parents=True,exist_ok=True)
        (folder/'cases.R').write_text(spec)
        with (folder/'golden.tsv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(golden[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(golden)
        (folder/'README.txt').write_text('Frozen expected results generated by independent Python Fraction centered-variance and full-subset enumeration. Replayed solely in R at package check time. Seed 40606. Regenerate with benchmarks/validate_native_r.py --cases 30 --golden; Python is a development oracle only.\n')


if __name__=='__main__':main()
