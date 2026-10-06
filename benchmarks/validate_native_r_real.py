"""Replay full-feature native R real-data ranks with centered Fraction variance."""
import argparse
import csv
from fractions import Fraction as F
from functools import cmp_to_key
import gzip
import json
from pathlib import Path
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
DIRECTIONS=('up','down','absolute')


def rank(scores,ids,direction):
    def compare(a,b):
        sa,qa=scores[a][:2];sb,qb=scores[b][:2]
        if direction=='absolute':sa,sb=abs(sa),abs(sb)
        if sa!=sb:c=(sa>sb)-(sa<sb)
        elif not sa:c=0
        elif qa is None or qb is None:c=sa*((qa is None)-(qb is None))
        else:c=sa*((qa>qb)-(qa<qb))
        if direction=='down':c=-c
        return -c if c else (ids[a]>ids[b])-(ids[a]<ids[b])
    return sorted(range(len(ids)),key=cmp_to_key(compare))


def centered_scores(features,labels,paired,deleted):
    retained=[i for i in range(len(labels)) if i not in deleted]
    pools=[[i for i in retained if labels[i]==group] for group in ([0] if paired else [0,1])]
    out=[]
    for xs in features:
        means=[];variances=[]
        for pool in pools:
            values=[xs[i] for i in pool];mean=sum(values,F(0))/len(values)
            means.append(mean)
            variances.append(sum(((x-mean)**2 for x in values),F(0))/(len(values)-1)/len(values))
        effect=means[0] if paired else means[0]-means[1]
        variance=sum(variances,F(0));out.append(((effect>0)-(effect<0),effect**2/variance if variance else None,effect))
    return out


WORKER=r'''
library(extremarank)
args<-commandArgs(TRUE); specs<-jsonlite::read_json(args[1],simplifyVector=FALSE); out<-list()
for(z in specs) {
  d<-read.delim(z$path,colClasses="character",check.names=FALSE)
  x<-matrix(as.numeric(as.matrix(d[,-1,drop=FALSE])),nrow=nrow(d),dimnames=list(d[[1]],names(d)[-1]))
  for(deleted in z$deletions) for(direction in c("up","down","absolute")) {
    d<-as.integer(unlist(deleted));r<-extremarank:::cpp_evaluate(x,rownames(x),as.integer(unlist(z$labels)),z$paired,direction,d)
    out[[length(out)+1L]]<-list(dataset=z$dataset,deleted=d,direction=direction,order=r$order)
  }
}
jsonlite::write_json(out,args[2],auto_unbox=TRUE,digits=NA,null="null")
'''


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--work',type=Path,required=True);a=parser.parse_args()
    a.work.mkdir(parents=True,exist_ok=True);start=time.perf_counter();specs=[];data={};reports={}
    for dataset in ('golub','kang-b','cptac','nutrimouse'):
        path=ROOT/'results/v04'/f'{dataset}-up/frozen_hex.tsv'
        if not path.exists():path=path.with_name(path.name+'.gz')
        opener=gzip.open if path.suffix=='.gz' else open
        with opener(path,'rt') as f:
            rows=list(csv.reader(f,delimiter='\t'));units=rows[0][1:];ids=[r[0] for r in rows[1:]]
            features=[[F(float.fromhex(x)) for x in r[1:]] for r in rows[1:]]
        paired=dataset=='kang-b'; labels=[0]*len(units)
        if not paired:
            with (ROOT/'results/v03'/f'{dataset}-up/groups.tsv').open() as f:meta={r['sample_id']:r['group'] for r in csv.DictReader(f,delimiter='\t')}
            target={'golub':'AML','cptac':'B','nutrimouse':'ppar'}[dataset];labels=[0 if meta[u]==target else 1 for u in units]
        deleted_sets={()}
        for direction in DIRECTIONS:
            z=json.loads((ROOT/'results/v04'/f'{dataset}-{direction}/audit.json').read_text());reports[dataset,direction]=z
            for row in z['features']:
                for prop in ('membership','direction'):
                    w=row[prop]['witness']
                    if w:
                        d=w['deleted_indices'];d=[d] if isinstance(d,int) else d
                        deleted_sets.add(tuple(d))
        specs.append(dict(dataset=dataset,path=str(path),paired=paired,labels=labels,deletions=sorted(deleted_sets,key=lambda d:(len(d),d))))
        data[dataset]=(ids,features,labels,paired)
    (a.work/'spec.json').write_text(json.dumps(specs));(a.work/'worker.R').write_text(WORKER)
    subprocess.run(['Rscript',str(a.work/'worker.R'),str(a.work/'spec.json'),str(a.work/'ranks.json')],check=True)
    actual={}
    for z in json.loads((a.work/'ranks.json').read_text()):
        d=z['deleted'];d=[d] if isinstance(d,int) else d
        actual[z['dataset'],tuple(d),z['direction']]=[j-1 for j in z['order']]
    records=[];checks=0
    for spec in specs:
        dataset=spec['dataset'];ids,features,labels,paired=data[dataset];orders={};signs={}
        for deleted in spec['deletions']:
            d=tuple(deleted);truth=centered_scores(features,labels,paired,{i-1 for i in d});signs[d]=[s[0] for s in truth]
            for direction in DIRECTIONS:
                expected=rank(truth,ids,direction);assert actual[dataset,d,direction]==expected
                orders[d,direction]=[ids[j] for j in expected];checks+=1
        witness_checks=0
        for direction in DIRECTIONS:
            z=reports[dataset,direction];k=z['top_k'];assert z['baseline_topk']==orders[(),direction][:k]
            for row in z['features']:
                for prop in ('membership','direction'):
                    w=row[prop]['witness']
                    if not w:continue
                    d=w['deleted_indices'];d=tuple(d if isinstance(d,list) else [d])
                    assert w['topk']==orders[d,direction][:k]
                    j=ids.index(row['feature_id'])
                    assert ((row['feature_id'] in w['topk'])!=row['baseline_topk']) if prop=='membership' else signs[d][j]!=row['baseline_effect_sign']
                    witness_checks+=1
        record=dict(dataset=dataset,features=len(ids),units=len(labels),distinct_baseline_and_witness_sets=len(spec['deletions']),
                    independent_complete_feature_ranks=len(spec['deletions'])*3,witness_checks=witness_checks,mismatches=0)
        records.append(record);print(json.dumps(record),flush=True)
    report=dict(datasets=records,independent_complete_feature_ranks=checks,mismatches=0,seconds=time.perf_counter()-start,
        scope='Every distinct native R reported witness and baseline replayed with Fraction means and centered sample variances; actual full-feature C++ ranks match; real unreported subsets are not exhaustively enumerated')
    (ROOT/'results/v04/real_fraction_validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)


if __name__=='__main__':main()
