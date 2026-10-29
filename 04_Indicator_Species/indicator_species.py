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
data_folder = Path("Data/Synthetic_Potato_Example")
output_folder = Path("Outputs/Synthetic/Indicator_Species")
community_type = "ITS"
group_column = "management"  # 24 organic and 36 conventional independent fields.
minimum_prevalence = 0.10
permutations = 999
random_seed = 42

# =====================================================
# STEP 1: Load and match samples
# =====================================================
output_folder.mkdir(parents=True, exist_ok=True)
filename = {"ITS": "ITS_OTU_table_final.tsv", "AMF": "AMF_OTU_table_final.tsv"}[community_type]
counts = pd.read_csv(data_folder / filename, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
counts["canonical"] = counts.index.str.replace(r"\.N\d+$", "", regex=True)
counts = counts.groupby("canonical").sum()
metadata = pd.read_csv(data_folder / "Final_data_with_diversity_prefixed.csv").groupby("site.id", as_index=True).first()
metadata = metadata.loc[counts.index.intersection(metadata.index)].dropna(subset=[group_column])
counts = counts.loc[metadata.index]
zero_reads = counts.sum(axis=1) == 0
counts, metadata = counts.loc[~zero_reads], metadata.loc[~zero_reads]
print(f"{community_type}: {len(counts)} sites used; {zero_reads.sum()} zero-read sites excluded.")

# =====================================================
# STEP 2: Define an example group and filter rare taxa
# =====================================================
# Management is predefined. Black-dot severity is continuous, so the previous
# threshold discarded information and left an unusably small high-severity group.
metadata["analysis_group"] = metadata[group_column]
if metadata["analysis_group"].nunique() < 2 or metadata["analysis_group"].value_counts().min() < 2:
    raise ValueError("Indicator analysis needs at least two groups with two or more independent sites each.")
counts = counts.loc[:, (counts > 0).mean(axis=0) >= minimum_prevalence]
counts.to_csv(output_folder / "community_for_R.csv")
metadata[["analysis_group"]].to_csv(output_folder / "groups_for_R.csv")

# =====================================================
# STEP 3: Run established IndVal.g implementation in R
# =====================================================
try:
    subprocess.run(["Rscript", "04_Indicator_Species/indicator_species.R", str(output_folder), str(permutations), str(random_seed)], check=True)
except FileNotFoundError:
    raise RuntimeError("Rscript is unavailable. Install R and indicspecies; no indicator result was created.")
except subprocess.CalledProcessError:
    raise RuntimeError("Indicator analysis failed; inspect the R error and group replication.")
