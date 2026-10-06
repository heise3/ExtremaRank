"""Independent Fraction means/centered variances and complete subset targets."""
import argparse
import json
from pathlib import Path
import random
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"src"),str(ROOT/"tests")]
from test_robustness import verify_problem


def main():
    p=argparse.ArgumentParser();p.add_argument("--cases",type=int,default=1000)
    p.add_argument("--output",type=Path,default=ROOT/"results"/"independent_robustness_v03.json")
    a=p.parse_args();rng=random.Random(60319);start=time.perf_counter();checks=0
    for i in range(a.cases):
        groups=["T"]*rng.randint(2,4)+["R"]*rng.randint(2,4)
        rows=[[rng.randint(-20,20)/8 for j in range(4)] for g in groups]
        for design in ("paired","welch"):
            for direction in ("up","down","absolute"):
                checks+=verify_problem(rows,groups,tuple("abcd"),direction,design,2,min(3,len(rows)-4))
    report={"cases":a.cases,"native_design_direction_problems":a.cases*6,"independent_property_and_scalar_enclosure_checks":checks,
            "failures":0,"seed":60319,"seconds":time.perf_counter()-start,
            "oracle":"Complete deletion subsets, Fraction means and centered sample variances, lexical ranking ties; checks statuses, exact minimum changes, witnesses and root and forced-inclusion scalar enclosures, plus capped-search intervals."}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=="__main__":main()
