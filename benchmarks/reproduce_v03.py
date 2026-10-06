"""Rebuild v0.3 real inputs, independent checks and reports from fixed sources."""
import argparse,gzip,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def run(*args):subprocess.run(args,cwd=ROOT,check=True)
def main():
 p=argparse.ArgumentParser();p.add_argument('--checks-only',action='store_true');a=p.parse_args()
 run(sys.executable,'data/fetch_v03_sources.py','--offline')
 if not a.checks_only:
  run(sys.executable,'benchmarks/prepare_v03_proteomics.py')
  run('Rscript','benchmarks/prepare_v03_lipids.R')
  folder=ROOT/'data/prepared/Kang10X'
  run('Rscript','benchmarks/prepare_v03_single_cell.R','data/source/v03/Kang/EH2259.rds',str(folder))
  for path in sorted(folder.iterdir()):
   if path.suffix in {'.tsv','.mtx'} and path.name!='expected_index.tsv':
    raw=path.read_bytes();path.with_name(path.name+'.gz').write_bytes(gzip.compress(raw,mtime=0));path.unlink()
 run(sys.executable,'benchmarks/independent_robustness_check.py','--cases','1000')
 run(sys.executable,'benchmarks/welch_pruning_v03.py')
 run(sys.executable,'benchmarks/validate_v03.py')
 run(sys.executable,'benchmarks/compact_v03.py')
if __name__=='__main__':main()
