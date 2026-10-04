"""Execute frozen v0.2 contract and independently replay every delete-one rank."""
import csv
from functools import cmp_to_key
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'tests')]
from extremarank.study import main as study_main
from extremarank.external import main as compare_main
from extremarank.preparation import prepare_study, read_table
from extremarank import PreparedWelch
from extremarank.unpaired import rank_welch
from test_unpaired import oracle_scores


def rank_fraction(scores, ids, direction):
    def compare(a,b):
        sa,qa,_ = scores[a]
        sb,qb,_ = scores[b]
        if direction == 'absolute':
            sa,sb = abs(sa),abs(sb)
        if sa != sb:
            c = (sa>sb)-(sa<sb)
        elif not sa:
            c = 0
        elif qa is None or qb is None:
            c = sa*((qa is None)-(qb is None))
        else:
            c = sa*((qa>qb)-(qa<qb))
        if direction == 'down':
            c = -c
        return -c if c else (ids[a]>ids[b])-(ids[a]<ids[b])
    return tuple(sorted(range(len(ids)),key=cmp_to_key(compare)))


def data_path(directory, name):
    plain = directory/name
    return plain if plain.exists() else directory/(name+'.gz')


def main():
    start = time.perf_counter()
    contract_path = ROOT/'benchmarks/upgrade_contract.json'
    contract = json.loads(contract_path.read_text())
    directory = ROOT/'data/prepared/Golub'
    matrix = data_path(directory,'matrix.csv')
    meta = directory/'metadata.tsv'
    study = prepare_study(matrix,meta,'AML','ALL',design='welch')
    runs = []
    for direction in contract['directions']:
        for budget in contract['budget']:
            out = ROOT/'results/upgrade'/f'welch-{direction}-b{budget}'
            study_main([str(matrix),'--metadata',str(meta),'--target','AML','--reference','ALL',
                        '--design','welch','--direction',direction,'--budget',str(budget),
                        '--top-k','20','--max-subsets','10000','--output',str(out)])
            report = json.loads((out/'audit.json').read_text())
            runs.append(dict(direction=direction,budget=budget,status=report['status'],
                             deleted_units=report['topk_audit']['deleted_units'],
                             witness_topk=report['topk_audit']['witness_topk'],
                             baseline_topk=report['topk_audit']['baseline_topk'],
                             subsets_checked=report['topk_audit']['subsets_checked'],
                             feasible_deletion_sets=report['topk_audit']['feasible_deletion_sets'],
                             delete_one_changed=report['delete_one_summary']['n_changed'],
                             seconds=report['timing_seconds']['study_total_before_writing']))
    print('Public Welch runs complete; checking all delete-one ranks independently.',flush=True)
    normalized_groups = tuple('T' if g=='AML' else 'R' for g in study.groups)
    prepared = PreparedWelch(study.values,study.groups,'AML','ALL')
    removed_sets = {(),*((i,) for i in range(len(study.unit_ids)))}
    for run in runs:
        if run['deleted_units'] is not None:
            removed_sets.add(tuple(study.unit_ids.index(uid) for uid in run['deleted_units']))
    exact_rank_checks = 0
    for deleted in sorted(removed_sets,key=lambda x:(len(x),x)):
        truth = oracle_scores(study.values,normalized_groups,deleted)
        actual_scores = prepared.evaluate(deleted)
        for direction in contract['directions']:
            expected = rank_fraction(truth,study.feature_ids,direction)
            actual = rank_welch(actual_scores,study.feature_ids,direction)
            assert actual == expected
            for run in runs:
                if run['direction'] == direction:
                    if not deleted:
                        assert run['baseline_topk'] == [study.feature_ids[j] for j in expected[:20]]
                    elif run['deleted_units'] == [study.unit_ids[i] for i in deleted]:
                        assert run['witness_topk'] == [study.feature_ids[j] for j in expected[:20]]
                        assert set(run['witness_topk']) != set(run['baseline_topk'])
            exact_rank_checks += 1
    print('Independent complete-feature Welch checks passed; importing limma refits.',flush=True)
    limma_dir = directory/'limma'
    eh,er,_ = read_table(limma_dir/'expected_topk.tsv')
    expected_topk = {}
    for row in er:
        rec = dict(zip(eh,row))
        expected_topk.setdefault((rec['scenario_id'],rec['direction']),[]).append(rec['feature_id'])
    rh,rr,_ = read_table(data_path(limma_dir,'refits.tsv'))
    refits = {}
    for row in rr:
        rec = dict(zip(rh,row))
        refits.setdefault(rec['scenario_id'],{})[rec['feature_id']] = float(rec['score'])
    bh,br,_ = read_table(limma_dir/'baseline.tsv')
    baseline = {dict(zip(bh,row))['feature_id']:float(dict(zip(bh,row))['score']) for row in br}
    external = []
    for direction in contract['directions']:
        out = ROOT/'results/upgrade'/f'limma-{direction}'
        compare_main([str(limma_dir/'baseline.tsv'),str(data_path(limma_dir,'refits.tsv')),
                      '--manifest',str(limma_dir/'manifest.tsv'),'--top-k','20','--direction',direction,
                      '--method','limma refitted moderated t; AML minus ALL','--output',str(out)])
        report = json.loads((out/'comparison.json').read_text())
        def ranking(values):
            return sorted(values,key=lambda fid:(values[fid] if direction=='down' else
                -abs(values[fid]) if direction=='absolute' else -values[fid],fid))[:20]
        assert ranking(baseline)==expected_topk['baseline',direction]==report['baseline_topk']
        for sid,scores in refits.items():
            assert ranking(scores)==expected_topk[sid,direction]
        external.append(dict(direction=direction,status=report['status'],scenarios=report['n_scenarios'],
                             changed=sum(row['topk_changed'] for row in report['scenarios']),
                             matched_R_topk_sets=1+len(refits)))
    output = dict(dataset='Golub ALL/AML, public multtest processed microarray training data',
                  samples=len(study.unit_ids),features=len(study.feature_ids),target_count=11,reference_count=27,
                  contract_sha256=hashlib.sha256(contract_path.read_bytes()).hexdigest(),
                  matrix_sha256=study.manifest['matrix_sha256'],metadata_sha256=study.manifest['metadata_sha256'],
                  welch_runs=runs,independent_full_feature_rank_checks=exact_rank_checks,
                  rank_mismatches=0,limma_runs=external,limma_R_rank_mismatches=0,
                  total_seconds=time.perf_counter()-start,
                  claim='ranking sensitivity only; independent designs verified numerically, no prediction/clinical efficacy')
    (ROOT/'results/upgrade_real.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':
    main()
