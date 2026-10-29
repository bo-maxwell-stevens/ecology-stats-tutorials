# Alpha diversity

Calculates observed richness, Shannon, Pielou evenness, Simpson, inverse Simpson, and read depth. Richness counts taxa; Shannon combines richness and relative abundance; Simpson emphasizes dominant taxa; Pielou is undefined for zero- and one-taxon communities. Empty libraries are excluded and written to a QC CSV. By default each sequence library is retained because `.N#` suffixes are not documented as technical replicates; enable the explicit aggregation setting only when that choice is justified.

Run: `python 01_Diversity/diversity_metrics.py`, `python 01_Diversity/reads_vs_richness.py`, then optionally `Rscript 01_Diversity/inext_coverage.R`.

`iNEXT` reports coverage-standardised richness at 90% sample coverage, only up to twice each sample's observed reads. A missing result means that target was not reached within that limit. Compare only samples that reach the target. Documentation: [iNEXT](https://cran.r-project.org/package=iNEXT).
