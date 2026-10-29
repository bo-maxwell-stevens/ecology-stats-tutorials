# Synthetic Potato Example

Every file here is completely fictional and reproducible with `python synthetic_data_generation.py` (seed `20261029`). It contains 60 fictional Danish-style fields, two biological soil samples per field (120 samples), 80 AMF virtual taxa (`VTX_SYN...`), and 200 ITS species hypotheses (`SH_SYN...`). No real DarkDivNet, Apsara, farm, coordinate, or taxon records were used.

`Final_data_with_diversity_prefixed.csv` contains sample and field IDs, replicate number, fictional coordinates/regions, pH, phosphorus (mg/100 g), magnesium (mg/100 g), C/N, organic carbon (%), moisture (%), management, and 0--100 harvest/storage severity plus qPCR outcomes for black dot (`Cc`) and silver scurf (`Hs`).

## Simulation ground truth

This is a teaching simulation, not biology: pH modestly lowers black-dot severity; crop rotation adds a modest black-dot effect; qPCR and severity are positively but imperfectly associated; selected ITS taxa increase with black-dot severity; selected AMF taxa decrease with it. Field effects, replicate noise, and multinomial sequencing variation leave substantial unexplained variation.

Counts arise from field/replicate community probabilities followed by multinomial sequencing. ITS and AMF library depths vary from zero to roughly 15,000 reads; three zero-read libraries occur in each table. Low depth therefore misses rare taxa. The metadata intentionally leave some ITS lifestyles empty or `unspecified`; functional-group analyses must retain those reads in the denominator.

## Glossary

**OTU/VT:** sequence-defined AMF taxon. **SH:** ITS species-hypothesis identifier. **Richness:** observed taxon count. **Evenness:** equality of taxon abundances. **Relative abundance:** reads divided by all reads in a sample. **Sequencing depth:** total reads. **Pseudoreplication:** treating related samples as independent.
