"""
TOPIC: Indicator species analysis

WHAT DOES THIS ANALYSIS DO? Finds taxa statistically associated with predefined groups.
WHEN SHOULD I USE IT? When groups have independent replicated samples.
INPUT: Community counts and a grouping variable. OUTPUT: IndVal.g results and a taxon plot.
MAIN STATISTIC: IndVal.g combines a taxon's specificity and fidelity to a group.
"""
from pathlib import Path
import subprocess
import pandas as pd

# =====================================================
# USER SETTINGS
# =====================================================
data_folder = Path("Data")
output_folder = Path("Outputs/Indicator_Species")
community_type = "ITS"
group_source_column = "pH_KCl"
minimum_prevalence = 0.10
permutations = 999

# =====================================================
# STEP 1: Load and match samples
# =====================================================
output_folder.mkdir(parents=True, exist_ok=True)
filename = {"ITS": "ITS_OTU_table_final.tsv", "AMF": "AMF_OTU_table_final.tsv"}[community_type]
counts = pd.read_csv(data_folder / filename, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
counts["canonical"] = counts.index.str.replace(r"\.N\d+$", "", regex=True)
counts = counts.groupby("canonical").sum()
metadata = pd.read_csv(data_folder / "Final_data_with_diversity_prefixed.csv").set_index("site.id")
metadata = metadata.loc[counts.index.intersection(metadata.index)].dropna(subset=[group_source_column])
counts = counts.loc[metadata.index]

# =====================================================
# STEP 2: Define an example group and filter rare taxa
# =====================================================
# DarkDivNet has no replicated management/disease group. This median pH split is a teaching example,
# not a claim that pH is a natural treatment. Replace with your disease/management grouping.
metadata["analysis_group"] = pd.Series("low_pH", index=metadata.index)
metadata.loc[metadata[group_source_column] > metadata[group_source_column].median(), "analysis_group"] = "high_pH"
counts = counts.loc[:, (counts > 0).mean(axis=0) >= minimum_prevalence]
counts.to_csv(output_folder / "community_for_R.csv")
metadata[["analysis_group"]].to_csv(output_folder / "groups_for_R.csv")

# =====================================================
# STEP 3: Run established IndVal.g implementation in R
# =====================================================
try:
    subprocess.run(["Rscript", "04_Indicator_Species/indicator_species.R", str(output_folder), str(permutations)], check=True)
except FileNotFoundError:
    print("Rscript is unavailable. Install R and indicspecies; prepared inputs were saved.")
except subprocess.CalledProcessError:
    print("Indicator analysis failed. Inspect R output and group replication.")
