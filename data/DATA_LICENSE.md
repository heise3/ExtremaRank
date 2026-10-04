# Third-party data and attribution

The repository's MIT license covers its original software and documentation.
It does not replace terms or attribution for third-party data and publications.
The bundled data are public NCBI-generated gene counts and GEO series/sample
metadata, with exact original URLs and retrieval hashes in `sources.json`.
Consult those sources and cite the original studies when reusing the data.
The prepared effects are reproducible transformations of those public counts.

- [GSE87290](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87290)
- [GSE50760](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760)
- [NCBI RNA-seq count pipeline](https://www.ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html)

The algorithm is usable on the included examples without fetching private or
restricted participant records. No original article PDF is redistributed.

## Golub public microarray data (v0.2)

`data/source/Golub/` retains the original `golub.RData`, dataset help,
preprocessing script and DESCRIPTION from the public Bioconductor
[multtest 2.68.0 source package](https://bioconductor.org/packages/3.23/bioc/src/contrib/multtest_2.68.0.tar.gz).
The package declares **LGPL**; this third-party material is not relicensed MIT.
Original package metadata and an LGPL license text are retained. The prepared
CSV, annotation and class labels are exported from the original processed
matrix; hashes and source URLs are in `data/source/Golub/provenance.json`.

Cite Golub et al. (1999), *Molecular classification of cancer: class discovery
and class prediction by gene expression monitoring*, Science 286:531–537,
and Dudoit, Fridlyand & Speed (2002), *Comparison of discrimination methods for
the classification of tumors using gene expression data*, JASA 97:77–87,
for the study and preprocessing. The matrix has 3,051 selected genes and 38
training samples, with 27 ALL and 11 AML labels. Generated sample IDs identify
the source column order; they do not add clinical metadata or verify identities.
The limma refit outputs are regenerated from this matrix; cite limma as
specified in its [author guide](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf).
