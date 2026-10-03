# Contributing

Run the exactness tests before proposing changes:

```bash
python -m unittest discover -s tests -v
python benchmarks/independent_theorem_check.py --cases 10000 --output results/local_check.json
```

Bug reports should include a minimal effects matrix, feature and donor IDs,
ranking direction, deletion budget, node limit, Python version, and the returned
status/witness. Keep raw participant data private when your data require it;
a synthetic minimal counterexample is usually sufficient.

New certificate targets need a stated statistic, tie rule, degeneracy policy,
sound bounds, and an independent exhaustive oracle. A bound overlap is not a
refutation. Changes to preprocessing are separate from deletion of fixed effects.
Performance contributions should compare identical targets and arithmetic,
include sorting/setup cost, and report unresolved cases.

The existing real data are computation checks, not a discovery benchmark.
Additional studies need original metadata, verified biological pairing, exact
source and preparation hashes, and a result-blind preprocessing contract.
