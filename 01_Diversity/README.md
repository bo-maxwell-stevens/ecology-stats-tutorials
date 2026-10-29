# Alpha diversity

Calculates observed richness, Shannon, Pielou evenness, Simpson, inverse Simpson, effective species number, and read depth. Use it for independent community samples. Richness is sensitive to sequencing depth, so inspect `reads_vs_richness.py` before interpreting differences.

Run: `python 01_Diversity/diversity_metrics.py`, `python 01_Diversity/reads_vs_richness.py`, then optionally `Rscript 01_Diversity/inext_coverage.R`.

`iNEXT` reports coverage-standardised richness at 90% sample coverage, only up to twice each sample's observed reads. A missing result means that target was not reached within that limit. Plot its CSV against observed richness if all samples can be standardised. Documentation: [iNEXT](https://cran.r-project.org/package=iNEXT).
