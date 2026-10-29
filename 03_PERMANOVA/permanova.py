"""
TOPIC: Multivariable PERMANOVA

WHAT DOES THIS ANALYSIS DO? Tests whether predictors are associated with community composition.
WHEN SHOULD I USE IT? For independent samples with a community table and predictors.
INPUT: AMF/ITS counts and metadata. OUTPUT: adonis2 and dispersion tables.
MAIN STATISTIC: Pseudo-F compares among-predictor to residual multivariate variation.
"""
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd

# =====================================================
# USER SETTINGS
# =====================================================
data_folder = Path("Data")
output_folder = Path("Outputs/PERMANOVA")
predictor_columns = ["pH_KCl", "P_Mehlich3_mg_100g"]
categorical_predictor = None  # Example for potato data: "Crop_rotation". Must have replication in every group.
permutations = 999
random_seed = 42
permutation_block_column = None  # Set only for a documented repeated-sample block.

# =====================================================
# STEP 1: Load data and match independent sampling units
# =====================================================
output_folder.mkdir(parents=True, exist_ok=True)
metadata = pd.read_csv(data_folder / "Final_data_with_diversity_prefixed.csv").set_index("site.id")
for community_name, filename in {"AMF": "AMF_OTU_table_final.tsv", "ITS": "ITS_OTU_table_final.tsv"}.items():
    counts = pd.read_csv(data_folder / filename, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
    counts["canonical"] = counts.index.str.replace(r"\.N\d+$", "", regex=True)
    counts = counts.groupby("canonical").sum()
    model_columns = predictor_columns + ([categorical_predictor] if categorical_predictor else [])
    selected_metadata = metadata.loc[counts.index.intersection(metadata.index), model_columns].dropna()
    selected_counts = counts.loc[selected_metadata.index]
    zero_reads = selected_counts.sum(axis=1) == 0
    excluded = selected_counts.index[zero_reads]
    selected_counts, selected_metadata = selected_counts.loc[~zero_reads], selected_metadata.loc[~zero_reads]
    pd.DataFrame({"excluded_zero_read_site": excluded}).to_csv(output_folder / f"{community_name}_excluded_zero_reads.csv", index=False)
    print(f"{community_name}: {len(selected_counts)} sites used; {len(excluded)} zero-read sites excluded.")

    # =====================================================
    # STEP 2: Write aligned inputs for vegan::adonis2
    # =====================================================
    selected_counts.to_csv(output_folder / f"{community_name}_community_for_R.csv")
    selected_metadata.to_csv(output_folder / f"{community_name}_metadata_for_R.csv")

    # =====================================================
    # STEP 3: Run multivariable PERMANOVA in R
    # =====================================================
    # adonis2 is used because it performs one multivariable model, unlike separate one-factor tests.
    blocks = "NONE"
    if permutation_block_column:
        metadata.loc[selected_counts.index, [permutation_block_column]].to_csv(output_folder / f"{community_name}_permutation_blocks.csv")
        blocks = str(output_folder / f"{community_name}_permutation_blocks.csv")
    command = ["Rscript", "03_PERMANOVA/permanova.R", community_name, str(output_folder), ",".join(model_columns), str(permutations), categorical_predictor or "NONE", blocks, str(random_seed)]
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as error:
        raise RuntimeError("Rscript is unavailable. Install R and vegan before running PERMANOVA.") from error
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f"vegan failed for {community_name}; no valid result was created.") from error
print("Null hypothesis: communities are unrelated to the listed predictors after the other model terms.")
