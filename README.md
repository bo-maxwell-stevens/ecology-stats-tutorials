# Ecology Statistics Tutorials

A self-contained introduction to microbial-community and potato-disease statistics for beginning graduate students. Every example uses a fictional potato-field dataset: 60 fields, two soil samples per field, AMF and ITS communities, soil chemistry, management, and simulated disease outcomes. Simulated associations are teaching examples, not biological findings.

## Start here

```bash
git clone https://github.com/bo-maxwell-stevens/ecology-stats-tutorials.git
cd ecology-stats-tutorials
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install jupyterlab
jupyter lab
```

On Windows PowerShell activate with `.venv\Scripts\Activate.ps1`. Open `Notebooks/01_Introduction_to_Ecological_Data_Analysis.ipynb`, then run scripts from the repository root, for example `python 01_Diversity/diversity_metrics.py`. R tutorials require `vegan`, `indicspecies`, `iNEXT`, `readr`, and `ggplot2`.

## Learning sequence

1. Notebook: tables, metadata, reads, richness, and Shannon diversity.
2. `01_Diversity`: sequencing depth and alpha diversity.
3. `02_Ordination`: NMDS and PCoA.
4. `03_PERMANOVA`: field-level community associations.
5. `04_Indicator_Species`: high/low black-dot teaching groups.
6. `08_Fungal_Functional_Groups`: annotated ITS proportions.
7. `06_Mixed_Effects_Models`: regression and repeated samples.
8. `05_Spatial_Autocorrelation`: Moran's I.
9. `07_Random_Forest`: grouped prediction and leave-one-field-out validation.

Outputs are regenerated under `Outputs/Synthetic/`. The independent unit for field-level questions is the field, not its two soil samples. Read `Data/Synthetic_Potato_Example/README.md` for the simulation design and glossary.
