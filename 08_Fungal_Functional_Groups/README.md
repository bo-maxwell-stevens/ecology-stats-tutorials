# Fungal functional groups

Matches ITS `SH` identifiers to `ITS_feature_metadata.tsv`, maps **primary** lifestyles into saprotroph, symbiotroph, and pathotroph/parasite teaching groups, and reports annotation coverage. The denominator is every ITS read, including unmatched and unclassified reads, so the stacked bars sum to 100% without inflating annotated groups. Secondary lifestyles are not double-counted; this discards some ecological complexity.

Run `python 08_Fungal_Functional_Groups/fungal_functional_groups.py`. The correlations use one field-level functional profile and mean field pH, so there are 60 independent fields rather than 120 related soil samples. The broad pathotroph/parasite category is not a plant-pathogen estimate. Correlations are exploratory, use BH correction, and concern compositional read proportions; reads are not biomass.
