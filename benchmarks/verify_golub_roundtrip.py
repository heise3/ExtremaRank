"""Verify every decimal value against bits independently exported from RData."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import struct

parser = argparse.ArgumentParser()
parser.add_argument('directory', type=Path)
args = parser.parse_args()
directory = args.directory
matrix = directory/'matrix.csv'
binary = directory/'matrix.binary64'
data = matrix.read_bytes() if matrix.exists() else gzip.decompress(matrix.with_suffix('.csv.gz').read_bytes())
bits = binary.read_bytes() if binary.exists() else gzip.decompress(binary.with_suffix('.binary64.gz').read_bytes())
rows = list(csv.reader(data.decode('utf-8').splitlines()))
values = [float(x) for row in rows[1:] for x in row[1:]]
source = struct.unpack('<'+'d'*len(values),bits)
mismatches = sum(a!=b for a,b in zip(values,source))
assert mismatches == 0
result = dict(values=len(values),mismatches=mismatches,
              source_bits_sha256=hashlib.sha256(bits).hexdigest(),
              decimal_csv_sha256=hashlib.sha256(data).hexdigest(),
              note='Python binary64 parser agrees exactly with original RData bits; R read.table can round one ulp differently')
(directory/'roundtrip_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
