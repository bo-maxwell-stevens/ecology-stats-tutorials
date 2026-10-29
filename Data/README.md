# Data

Raw DarkDivNet tables are deliberately not included: they may be large and contain field coordinates. For permitted local use, run this from the repository root to copy the four supplied files into this directory:

```bash
cp ../Data/AMF_OTU_table_final.tsv ../Data/ITS_OTU_table_final.tsv \
   ../Data/ITS_feature_metadata.tsv ../Data/Final_data_with_diversity_prefixed.csv Data/
```

The source path above matches this workspace, where the repository is alongside the supplied `Data/` folder. Adjust the source path for a different location. Do not commit these files without data-owner permission.

Expected files:

| File | Purpose |
|---|---|
| `AMF_OTU_table_final.tsv` | AMF counts; rows are sequence samples, columns are VTX taxa. |
| `ITS_OTU_table_final.tsv` | ITS counts; rows are sequence samples, columns are UNITE SH taxa. |
| `ITS_feature_metadata.tsv` | Taxonomy and ecological annotations for ITS SH identifiers. |
| `Final_data_with_diversity_prefixed.csv` | Environmental/sample metadata. |

In these data, a sequence ID such as `D018.N1` maps to `site.id == D018` after removing the `.N1` suffix. Some sites have multiple sequencing samples but only one environmental record. This is not repeated environmental sampling.
