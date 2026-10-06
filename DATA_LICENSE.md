# Third-party data and attribution

ExtremaRank source code is MIT-licensed. That license **does not relicense**
third-party datasets, provider annotations, package examples or help documents.
Data remain under their original terms. Source files are separated from code,
with pinned versions, original URLs, byte hashes and preparation contracts.
No third-party analysis code is included in the Python wheel.

| Bundled source | Origin and attribution | Source terms and limitations |
|---|---|---|
| GSE87290 / GSE50760 | Public GEO count matrices; generators and retrieval in [DATA.md](docs/DATA.md) | No new per-file MIT/CC license inferred from public availability |
| Golub microarray matrix | Golub et al., *Molecular classification of cancer* (1999), processed `multtest` example | Provider LGPL-2 package license and bundled provenance remain applicable |
| Kang EH2259 | Kang et al., *Multiplexed droplet single-cell RNA-sequencing using natural genetic variation* (2018), GSE96583; Bioconductor `muscData` Kang18_8vs8 | muscData MIT package; original study/provider data attribution retained; raw genotype checks not repeated |
| muscat example_sce | Bioconductor muscat 1.26.0 package example, Crowell et al., *muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data* (2020) | GPL-3 package. Only extracted data/description included; insufficient independent donors for the declared nonzero paired deletion example |
| CPTAC spike-in protein intensities | Public statOmics/msqrob2data `dda/cptacAvsB_lab3`, commit 6338dcfebea8681b20f6f30431b865a37d3ef428 | Repository did not declare a license at retrieval. Measurement tables and sample annotation retained with source attribution; author prose and analysis code omitted. Do not infer that MIT applies |
| Nutrimouse | mixOmics 6.36.0 `nutrimouse`; Martin et al., *Journal of Nutritional Biochemistry* 18, 362–371 (2007) | Package GPL (>=2); extracted dataset and GPL help/description retain those terms. GPL-3 text supplied as a permitted version. These are fatty-acid percentages, not absolute concentrations |

Original metadata and hashes are in `data/source/v03/`. GNU GPL text is bundled
at `data/source/v03/GPL-3.txt`; earlier source-specific license files are retained.
Downloaded full package archives are not bundled in v0.3; selected bytes can be
restored and verified by `data/fetch_v03_sources.py`. The original EH2259 resource
bytes are included for offline reproduction.

Primary source links:

- [muscData provider](https://github.com/HelenaLC/muscData), [EH2259](https://experimenthub.bioconductor.org/fetch/2259), [GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583).
- [CPTAC measurement folder at fixed commit](https://github.com/statOmics/msqrob2data/tree/6338dcfebea8681b20f6f30431b865a37d3ef428/dda/cptacAvsB_lab3).
- [muscat official package](https://bioconductor.org/packages/muscat/).
- [mixOmics official package](https://bioconductor.org/packages/mixOmics/), [official Nutrimouse case study](https://mixomics.org/wp-content/uploads/2025/01/rCCA-Nutrimouse-Case-Study.html).
