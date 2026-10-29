# Ecological Statistics Tutorials

Beginner-friendly Python examples for common microbial-ecology analyses. The examples use the DarkDivNet table layout locally, but every script starts with settings you can change for another study. Educational clarity takes priority over sophisticated modelling.

## Who is this for?

Undergraduate and graduate students, and researchers learning ecological statistics with Python.

## Installation

```bash
git clone YOUR-REPOSITORY-URL ecological-statistics-tutorials
cd ecological-statistics-tutorials
python3 -m venv .venv
source .venv/bin/activate              # Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

Copy permitted input tables as described in [`Data/README.md`](Data/README.md). Tutorials 03 and 04 use R because `vegan::adonis2` and `indicspecies::multipatt` are the established implementations:

```r
install.packages(c("vegan", "indicspecies", "iNEXT"))
```

## Getting started

Run `01_Diversity/diversity_metrics.py` first. It creates diversity values from count data, which are used by later examples. Run scripts from the repository root, for example `python 02_Ordination/nmds_pcoa.py`.

| Tutorial | Statistical method | Research question |
|---|---|---|
| 01 | Alpha diversity | How diverse are microbial communities? |
| 02 | Ordination | How do communities differ among sites? |
| 03 | PERMANOVA | Which predictors are associated with composition? |
| 04 | Indicator species | Which taxa are associated with groups? |
| 05 | Moran's I | Are observations spatially autocorrelated? |
| 06 | Mixed-effects models | How do predictors relate to a response with grouped observations? |
| 07 | Random Forest | Can we predict ecological responses? |
| 08 | Functional groups | How do fungal ecological groups vary? |

## DarkDivNet example variables

The available environmental table has `pH_KCl`, `P_Mehlich3_mg_100g`, `K_Mehlich3_mg_100g`, `N_pct`, `C_pct`, climate variables, `lat`, `lon`, and a unique `site.id` for each environmental record. Its `canonical` value retains a sequence replicate suffix, while `site.id` is the environmental matching key. It has no potato disease severity, cultivar, crop rotation, management category, or repeated environmental observations within site. The scripts therefore use pH and Mehlich-3 phosphorus as executable examples. They explain where to replace them with your disease or management variables.

## Adapting to your study

At the top of each script, change the `USER SETTINGS`: filenames, `response_variable`, `predictor_columns`, and grouping column. Community tables require samples as rows and taxa as columns. Metadata requires a sample ID column. For potato data, replace `pH_KCl` or `P_Mehlich3_mg_100g` with your disease-severity and management columns only after checking their type, missingness, and independent sampling unit.

## Statistical limitations

Correlation is not causation. Avoid pseudoreplication: samples from one field are not automatically independent fields. Correct many taxon-level tests for multiple testing. Sequencing depth and compositional counts affect diversity and community analyses. Small samples make estimates unstable and can overfit random forests. Assess spatial autocorrelation and validate predictions without leaking held-out sites or samples into preprocessing. Correlated predictors make individual effects difficult to interpret.

## Reproducibility

All scripts set a random seed and save CSV results plus 300 dpi PNG and PDF figures under `Outputs/`. Raw data and outputs are ignored by Git. Package versions are bounded in `requirements.txt`; record the exact environment with `pip freeze > requirements-lock.txt` for a publication.

## License

The code is MIT licensed. This license does not grant rights to third-party datasets.
