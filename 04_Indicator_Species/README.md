# Indicator species

`indicspecies::multipatt` calculates IndVal.g with permutation tests and Benjamini-Hochberg FDR correction. Prevalence filtering reduces the burden from extremely rare taxa. An indicator is associated with a group, not proven to cause that group's characteristics.

Run `python 04_Indicator_Species/indicator_species.py`. The example creates a median pH teaching split because the available data have no replicated disease or management group. Replace it with a scientifically defined variable such as high/low disease severity. Documentation: [indicspecies](https://cran.r-project.org/package=indicspecies).
