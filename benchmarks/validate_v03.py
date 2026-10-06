"""Execute frozen real contracts, check independent sums and replay witnesses."""
import csv, gzip, json, os, sys, time
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'tests'),str(ROOT/'benchmarks')]
from extremarank.study import main as study_main
from extremarank.pseudobulk import aggregate_pseudobulk
from extremarank.preparation import prepare_study, read_table
from extremarank.models import run_refits
from test_unpaired import oracle_scores
from upgrade_real import data_path, rank_fraction
DIRECTIONS=('up','down','absolute')

def table(path):
 h,r,_=read_table(path);return h,[dict(zip(h,x)) for x in r]

def check_sums():
 source=ROOT/'data/prepared/Kang10X';out=ROOT/'data/prepared/KangPB'
 pb=aggregate_pseudobulk(source.relative_to(ROOT),(source/'cell_metadata.tsv.gz').relative_to(ROOT),out.relative_to(ROOT),'10x',10)
 _,index=table(source/'expected_index.tsv');checks=0
 for record,cell in zip(index,pb['cell_types']):
  assert record['cell_type']==cell['cell_type']
  folder=out/cell['directory'];_,meta=table(folder/'metadata.tsv')
  ids={x['sample_id']:x['donor_id']+'|'+x['group'] for x in meta}
  eh,er=table(data_path(source,record['filename']));ah,ar=table(data_path(folder,'matrix.tsv'))
  assert set(eh[1:])==set(ids.values()) and len(er)==len(ar)
  for expected,actual in zip(er,ar):
   assert expected['feature_id']==actual['feature_id']
   for sid,key in ids.items():
    assert int(actual[sid])==int(expected[key]);checks+=1
 print('Independent R pseudobulk sums matched:',checks,flush=True)
 return {'genes':pb['n_genes'],'provider_singlet_cells':pb['n_cells'],'retained_cells':pb['retained_cells'],'pseudobulks':pb['n_pseudobulks'],'exact_integer_sums_checked':checks,'mismatches':0,'cell_types':len(pb['cell_types'])}

def native(name,matrix,meta,target,reference,design,k,transform='none'):
 matrix,meta=Path(matrix).relative_to(ROOT),Path(meta).relative_to(ROOT)
 prepared=prepare_study(matrix,meta,target,reference,design,transform=transform)
 runs=[];witnesses={():None}
 for direction in DIRECTIONS:
  out=ROOT/'results/v03'/f'{name}-{direction}'
  study_main([str(matrix),'--metadata',str(meta),'--target',target,'--reference',reference,'--design',design,'--transform',transform,'--top-k',str(k),'--budget','2','--max-nodes','10000','--max-subsets','10000','--feature-audit','--direction',direction,'--output',str(out)])
  r=json.loads((out/'feature_robustness.json').read_text())
  for row in r['features']:
   for prop in ('membership','direction'):
    if row[prop]['witness']:witnesses[tuple(row[prop]['witness']['deleted_indices'])]=None
  runs.append({'direction':direction,'status':r['status'],'membership':dict(Counter(x['membership']['status'] for x in r['features'])),'effect_sign':dict(Counter(x['direction']['status'] for x in r['features'])),'minimum_topk_change':r['topk_minimum_change'],'scenarios_checked':r['scenarios_checked'],'output':str(out.relative_to(ROOT))})
 # Definition-based independent Fraction means/centered variances. For paired
 # scores compare the effect rows with a constant-zero reference of size two.
 rows=list(prepared.values);groups=['T' if x==target else 'R' for x in prepared.groups]
 if design=='paired':rows += [[0.0]*len(prepared.feature_ids)]*2;groups=['T']*len(prepared.values)+['R']*2
 orders={}; signs={}
 for deleted in witnesses:
  truth=oracle_scores(rows,groups,deleted)
  signs[deleted]={fid:s[0] for fid,s in zip(prepared.feature_ids,truth)}
  for direction in DIRECTIONS:
   order=rank_fraction(truth,prepared.feature_ids,direction)
   orders[deleted,direction]=[prepared.feature_ids[j] for j in order[:k]]
 for run in runs:
  r=json.loads((ROOT/run['output']/'feature_robustness.json').read_text())
  assert r['baseline_topk']==orders[(),run['direction']]
  for row in r['features']:
   for prop in ('membership','direction'):
    w=row[prop]['witness']
    if not w:continue
    d=tuple(w['deleted_indices']);assert w['topk']==orders[d,run['direction']]
    if prop=='membership':assert (row['feature_id'] in w['topk']) != row['baseline_topk']
    else:assert signs[d][row['feature_id']]!=row['baseline_effect_sign']
 print('Full-feature Fraction witness replay complete:',name,len(witnesses)*3,flush=True)
 return {'dataset':name,'units':len(prepared.unit_ids),'features':len(prepared.feature_ids),'native_runs':runs,'independent_full_feature_rank_checks':len(witnesses)*3,'rank_mismatches':0,'certificate_validation':'Exhaustive independent small-problem tests support proofs; real baseline and every distinct reported witness replayed independently. Real unreported subsets not enumerated by this checker.'}

def models(name,matrix,meta,target,reference,model,k,paired=False,categorical=(),positive=0):
 out=ROOT/'results/v03'/name
 r=run_refits(matrix,meta,target,reference,model,out,k,1,'up',categorical,(),paired,min_count_samples=positive)
 _,expected_rows=table(out/'expected_topk.tsv');expected={}
 for row in expected_rows:expected.setdefault((row['scenario_id'],row['direction']),[]).append(row['feature_id'])
 all_scores={};checks=0
 for file in ('baseline.tsv','refits.tsv'):
  if not (out/file).exists():continue
  _,rows=table(out/file)
  for row in rows:all_scores.setdefault(row.get('scenario_id','baseline'),{})[row['feature_id']]=float(row['score'])
 for (sid,direction),topk in expected.items():
  values=all_scores[sid]
  actual=sorted(values,key=lambda fid:(values[fid] if direction=='down' else -abs(values[fid]) if direction=='absolute' else -values[fid],fid))[:k]
  assert actual==topk;checks+=1
 statuses=r.get('fit_status',[])
 record={'dataset_model':name,'status':r['status'],'requested_refits':r['requested_refits'],'valid_refits':sum(x['status']=='OK' and x['scenario_id']!='baseline' for x in statuses),'failed_refits':sum(x['status']!='OK' and x['scenario_id']!='baseline' for x in statuses),'baseline_valid':any(x['scenario_id']=='baseline' and x['status']=='OK' for x in statuses),'matched_R_topk_sets':checks,'R_topk_mismatches':0,'retained_features':r['fixed_feature_filter']['retained'],'output':str(out.relative_to(ROOT))}
 print('R model:',json.dumps(record),flush=True);return record

def main():
 os.chdir(ROOT);start=time.perf_counter();pb=check_sums();records=[]
 g=ROOT/'data/prepared/Golub';b=ROOT/'data/prepared/KangPB/celltype_1';c=ROOT/'data/prepared/CPTAC_runs';l=ROOT/'data/prepared/Nutrimouse'
 for args in [('golub',data_path(g,'matrix.csv'),g/'metadata.tsv','AML','ALL','welch',20),('kang-b',data_path(b,'matrix.tsv'),b/'metadata.tsv','stim','ctrl','paired',20,'cpm-log2'),('cptac',data_path(c,'matrix.tsv'),c/'metadata.tsv','B','A','welch',20),('nutrimouse',data_path(l,'matrix.tsv'),l/'metadata.tsv','ppar','wt','welch',5)]:records.append(native(*args))
 refits=[models('golub-limma',data_path(g,'matrix.csv'),g/'metadata.tsv','AML','ALL','limma',20)]
 for model in ('limma-voom','edgeR','DESeq2'):refits.append(models('kang-b-'+model,data_path(b,'matrix.tsv'),b/'metadata.tsv','stim','ctrl',model,20,True,positive=3))
 refits.append(models('nutrimouse-limma-diet',data_path(l,'matrix.tsv'),l/'metadata.tsv','ppar','wt','limma',5,categorical=('diet',)))
 report={'version':'0.3.0','pseudobulk':pb,'native':records,'automatic_refits':refits,'seconds':time.perf_counter()-start,'claim':'Sample-deletion ranking sensitivity only; technical protein runs and biological lipid/microarray/donor units have different semantics. No biological truth, prediction, FDR or clinical efficacy claim.'}
 (ROOT/'results/validation_v03_real.json').write_text(json.dumps(report,indent=2)+'\n');print('Completed real validation in',report['seconds'],'seconds',flush=True)
if __name__=='__main__':main()
