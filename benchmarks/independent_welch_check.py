"""Independent Fraction enumeration against exact integer Welch rankings."""
import argparse
import json
from pathlib import Path
import random
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'tests')]
from test_unpaired import oracle_rank, feasible_deletions
from extremarank import PreparedWelch, audit_welch
from extremarank.unpaired import rank_welch


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cases', type=int, default=2000)
    parser.add_argument('--output', type=Path, default=ROOT/'results/independent_welch_check.json')
    args = parser.parse_args()
    rng = random.Random(740125)
    started = time.perf_counter()
    evaluations = audits = 0
    for trial in range(args.cases):
        n, m = rng.randint(2,4), rng.randint(2,4)
        groups = ['T']*n+['R']*m
        rows = [[rng.randint(-20,20)/8 for _ in range(4)] for _ in groups]
        ids = ('10','2','a','z')
        budget = min(2,n+m-4)
        prepared = PreparedWelch(rows,groups,'T','R')
        deletes = list(feasible_deletions(groups,budget))
        assert set(deletes) == set(prepared.deletions(budget))
        for direction in ('up','down','absolute'):
            truth = oracle_rank(rows,groups,ids,direction)
            selected = set(truth[:2])
            changed = False
            for deleted in [(),*deletes]:
                expected = oracle_rank(rows,groups,ids,direction,deleted)
                assert rank_welch(prepared.evaluate(deleted),ids,direction) == expected
                changed |= set(expected[:2]) != selected
                evaluations += 1
            result = audit_welch(rows,groups,ids,'T','R',2,budget,10000,direction)
            assert result.status == ('REFUTED' if changed else 'CERTIFIED')
            if result.deleted_indices is not None:
                assert result.deleted_indices in deletes
                assert tuple(ids[j] for j in oracle_rank(rows,groups,ids,direction,result.deleted_indices)[:2]) == result.witness_topk
            audits += 1
    output = dict(cases=args.cases,rank_evaluations=evaluations,complete_status_checks=audits,
                  failures=0,seed=740125,arithmetic='independent Fraction means and centered sample variances',
                  seconds=time.perf_counter()-started)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output))


if __name__ == '__main__':
    main()
