# Release checklist and repository copy

The public repository is [heise3/ExtremaRank](https://github.com/heise3/ExtremaRank).
The owner explicitly authorized uploading this implementation and reproducibility
package. Published versions and downloadable artifacts appear on its
[release page](https://github.com/heise3/ExtremaRank/releases).
No PyPI release, DOI or journal submission has been created by this package.
Use `ExtremaRank` as a working project name; registry-name availability has not
been checked. Replace contributor metadata in CITATION.cff before DOI registration.

Suggested GitHub description:

> Sample-deletion robustness audits for paired and independent-group omics: exact paired extrema, Welch audits, external refit comparison, and offline diagnostic reports.

Suggested topics:

`bioinformatics`, `rna-seq`, `proteomics`, `metabolomics`, `paired-data`,
`robustness`, `exact-optimization`, `reproducible-research`, `python`.

The repository uses the verified authenticated owner's GitHub identity.
No verified DOI is invented in the package.
GitHub Actions contains the 3.10/3.12/3.14 test matrix and executes on pushes.
Local tests on 3.12 and 3.14 were executed; remote run status is available from
the repository's [Actions](https://github.com/heise3/ExtremaRank/actions).
Action versions follow the verified official [checkout](https://github.com/actions/checkout)
and [setup-python](https://github.com/actions/setup-python) documentation at packaging.

Useful project presentation is already included: a short runnable example,
full theorem, documented closest methods, one-command real-data reproduction,
machine-readable witnesses, an offline interactive report, and an editable SVG.
v0.2 adds independent two-group studies, strict matrix/metadata preparation,
delete-one diagnostics and an executed limma refit-table adapter. New local
validation is recorded in `results/VALIDATION_V02.json`; `VALIDATION.json`
retains the historical v0.1 record. Scientific claims remain restricted to the
workflow-specific targets in `docs/STUDY.md`.
