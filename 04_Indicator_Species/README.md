# Indicator species

`indicspecies::multipatt` calculates IndVal.g with permutation tests and Benjamini-Hochberg FDR correction. The synthetic example compares predefined organic (24 fields) and conventional (36 fields) management groups. The former black-dot severity threshold was unsuitable because it created only two high-severity fields. Disease severity is continuous, so categorising it also loses information.

Run `python 04_Indicator_Species/indicator_species.py`. Results include group, IndVal.g statistic, raw and FDR-adjusted p-values, prevalence, and mean abundance. An indicator taxon is associated with a group, not necessarily causal; no significant taxa is a valid result. Documentation: [indicspecies](https://cran.r-project.org/package=indicspecies).
