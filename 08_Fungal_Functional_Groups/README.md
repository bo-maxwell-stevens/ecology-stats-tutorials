# Fungal functional groups

Matches ITS `SH` identifiers to `ITS_feature_metadata.tsv`, maps **primary** lifestyles into saprotroph, symbiotroph, and putative-pathogen teaching groups, and reports how much of the count table was annotated. Secondary lifestyles are not double-counted; this is explicit but necessarily discards some ecological complexity.

Run `python 08_Fungal_Functional_Groups/fungal_functional_groups.py`. Change `environmental_predictor` to disease severity in a potato study. Sequence reads are not direct organism abundance or biomass. The annotations should be reviewed for the taxonomy/database version used.
