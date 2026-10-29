# Alpha diversity

Calculates observed richness, Shannon, Hill diversity q = 1 (exponential Shannon, the effective number of taxa), Pielou evenness, Simpson, inverse Simpson, and read depth. Shannon is entropy; exponentiating it expresses the same diversity on a taxon-count scale. Pielou is undefined for zero- and one-taxon communities. Empty libraries are excluded and written to a QC CSV. By default each sequence library is retained because `.N#` suffixes are not documented as technical replicates; enable the explicit aggregation setting only when that choice is justified.

Run: `python 01_Diversity/diversity_metrics.py`, `python 01_Diversity/reads_vs_richness.py`, then optionally `Rscript 01_Diversity/inext_coverage.R`.

`iNEXT` reports the first richness estimate reaching at least 90% sample coverage, only up to twice each sample's observed reads. It uses zero bootstrap replicates here for a fast descriptive demonstration, so it has no confidence intervals and is not an inferential comparison. A missing result means the target was not reached within that limit. Documentation: [iNEXT](https://cran.r-project.org/package=iNEXT).
