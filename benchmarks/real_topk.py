"""Practical shared-donor top20 counterexamples, with independent replay."""
from __future__ import annotations
import csv
from fractions import Fraction
from functools import cmp_to_key
import gzip
import hashlib
import json
from pathlib import Path
import time
from extremarank import audit_topk

ROOT = Path(__file__).resolve().parents[1]


def independent_topk(rows, ids, keep, k, direction="up"):
    # Separate Fraction implementation: no Score, extrema, or ranking helper.
    moments = []
    for j in range(len(ids)):
        xs = [Fraction.from_float(rows[i][j]) for i in keep]
        moments.append((sum(xs), sum(x*x for x in xs)))
    def compare(a,b):
        s,q = moments[a]
        t,r = moments[b]
        if direction == "absolute":
            s,t = abs(s),abs(t)
        sa,sb = (s>0)-(s<0),(t>0)-(t<0)
        if sa != sb:
            c = (sa>sb)-(sa<sb)
        elif not sa:
            c = 0
        else:
            l,h = s*s*r,t*t*q
            c = sa*((l>h)-(l<h))
        if direction == "down":
            c = -c
        return -c if c else ((ids[a]>ids[b])-(ids[a]<ids[b]))
    return tuple(ids[j] for j in sorted(range(len(ids)),key=cmp_to_key(compare))[:k])


def main():
    payload = {"target": "full frozen eligible signed paired-t top20 set",
               "max_deletions": 4, "datasets": []}
    for accession in ("GSE87290", "GSE50760"):
        path = ROOT / "data" / "prepared" / accession / "effects.csv.gz"
        with gzip.open(path,"rt",newline="") as f:
            reader = csv.reader(f)
            ids = tuple(next(reader)[1:])
            donors, rows = [], []
            for row in reader:
                donors.append(row[0]); rows.append(tuple(float(x) for x in row[1:]))
        for direction in ("up", "down", "absolute"):
            start = time.perf_counter()
            result = audit_topk(rows, ids, 20, 4, direction=direction)
            seconds = time.perf_counter()-start
            baseline = independent_topk(rows,ids,tuple(range(len(rows))),20,direction)
            assert baseline == result.baseline_topk
            entry = {"accession": accession, "direction": direction,
                     "donors": len(rows), "genes":len(ids),
                     "effects_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                     "audit_seconds": seconds, "audit": result.as_dict(),
                     "baseline_independently_verified": True}
            if result.status == "REFUTED":
                keep = tuple(i for i in range(len(rows)) if i not in result.deleted_indices)
                witness = independent_topk(rows,ids,keep,20,direction)
                assert witness == result.witness_topk and set(witness) != set(baseline)
                entry.update(witness_independently_verified=True,
                             deleted_donor_ids=[donors[i] for i in result.deleted_indices],
                             entered_gene_ids=sorted(set(witness)-set(baseline)),
                             exited_gene_ids=sorted(set(baseline)-set(witness)))
            payload["datasets"].append(entry)
            print(json.dumps({k:v for k,v in entry.items() if k!="audit"}),flush=True)
            (ROOT/"results"/"real_topk.json").write_text(json.dumps(payload,indent=2)+"\n")


if __name__ == "__main__":
    main()
