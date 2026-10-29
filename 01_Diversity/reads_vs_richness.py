"""
TOPIC: Sequencing depth and richness

WHAT DOES THIS ANALYSIS DO? Tests whether samples with more reads have more observed taxa.
WHEN SHOULD I USE IT? Before comparing observed richness among samples.
INPUT: Community count tables. OUTPUT: read-depth plots and Spearman correlations.
MAIN STATISTIC: Spearman rho measures monotonic association without assuming normality.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

# =====================================================
# USER SETTINGS
# =====================================================
dataset = "synthetic"
data_folder = Path("Data/Synthetic_Potato_Example") if dataset == "synthetic" else Path("Data")
output_folder = Path("Outputs/Diversity")
community_files = {"AMF": "AMF_OTU_table_final.tsv", "ITS": "ITS_OTU_table_final.tsv"}

# =====================================================
# STEP 1: Load and summarize the data
# =====================================================
output_folder.mkdir(parents=True, exist_ok=True)
sns.set_theme(style="whitegrid", context="talk")
results = []
for community_name, filename in community_files.items():
    counts = pd.read_csv(data_folder / filename, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
    summary = pd.DataFrame({"total_reads": counts.sum(axis=1), "observed_richness": (counts > 0).sum(axis=1)})

    # =====================================================
    # STEP 2: Test the read-depth relationship
    # =====================================================
    rho, p_value = spearmanr(summary["total_reads"], summary["observed_richness"])
    results.append({"community": community_name, "n_samples": len(summary), "spearman_rho": rho, "p_value": p_value})

    # =====================================================
    # STEP 3: Plot and save the results
    # =====================================================
    figure, axis = plt.subplots(figsize=(7, 5))
    sns.regplot(data=summary, x="total_reads", y="observed_richness", lowess=True, scatter_kws={"alpha": .7}, ax=axis)
    axis.set(xlabel="Total sequencing reads", ylabel="Observed richness", title=f"{community_name}: reads versus richness\nSpearman rho = {rho:.2f}, p = {p_value:.3g}")
    figure.tight_layout()
    figure.savefig(output_folder / f"{community_name}_reads_vs_richness.png", dpi=300)
    figure.savefig(output_folder / f"{community_name}_reads_vs_richness.pdf")
    plt.close(figure)
pd.DataFrame(results).to_csv(output_folder / "reads_vs_richness_spearman.csv", index=False)
print("Observed richness can increase with read depth; use iNEXT coverage standardisation below when appropriate.")
